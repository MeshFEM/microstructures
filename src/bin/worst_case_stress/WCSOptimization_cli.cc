////////////////////////////////////////////////////////////////////////////////
// WCSOptimization_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evolves the microstructure parameters to bring the structure's
//      homogenized elasticity tensor closer to a target while also minimizing
//      worst-case stress.
//
//      Example invocation using BoundaryPerturbationInflator:
//      ./WCSOptimization_cli -p demo.msh
//              -m $MICRO_DIR/docs/materials/B9Creator.material ./test_2D_job.opt
//              --alpha 0 --inflator boundary_perturbation -o iterate -d1
//              --step 1e-21 -P4 -n 80
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/12/2014 01:15:28
////////////////////////////////////////////////////////////////////////////////
#include "WCSObjectiveTerm.hh"

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/Materials.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/Future.hh>

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <functional>

#include <json.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <inflators/Inflator.hh>
#include <inflators/MakeInflator.hh>

#include <pattern_optimization/PatternOptimizationJob.hh>
#include <pattern_optimization/PatternOptimizationConfig.hh>
#include <pattern_optimization/IterateFactory.hh>
#include <pattern_optimization/IterateManager.hh>
#include <pattern_optimization/PatternOptimizationIterate.hh>

#include <optimizers/BoundConstraints.hh>
#include <optimizers/wrappers/ceres.hh>
#include <optimizers/wrappers/dlib.hh>
#include <optimizers/wrappers/gradient_descent.hh>
#include <optimizers/wrappers/nlopt.hh>
#include <optimizers/wrappers/knitro.hh>

#include <optimizers/OptimizerConfig.hh>
#include <pattern_optimization/PatternOptimizationConfig.hh>
#include <pattern_optimization/objective_terms/TensorFit.hh>
#include <pattern_optimization/objective_terms/IsotropicFit.hh>
#include <pattern_optimization/objective_terms/IsotropicFitRel.hh>
#include <pattern_optimization/objective_terms/ProximityRegularization.hh>
#include <pattern_optimization/objective_terms/TargetVolume.hh>
#include <pattern_optimization/objective_terms/PeriodicSmoothingRegularization.hh>

#include <pattern_optimization/constraints/TensorFit.hh>
#include <pattern_optimization/constraints/Printability.hh>

#include <MeshFEM/Parallelism.hh>

namespace po = boost::program_options;
namespace PO = PatternOptimization;
using json = nlohmann::json;
using namespace std;

[[ noreturn ]] void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: WCSOptimization_cli [options] job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

using OptimizerMap =
    map<string, std::function<void(ScalarField<Real> &, const PO::BoundConstraints &,
            PO::IterateManagerBase &, const PO::OptimizerConfig &, const string &)>>;

OptimizerMap optimizers = {
    {"levenberg_marquardt",  optimize_ceres_lm},
    {"dogleg",               optimize_ceres_dogleg},
    {"bfgs",                 optimize_dlib_bfgs},
    {"custom_bfgs",          optimize_dlib_custom_bfgs},
    {"slsqp",                optimize_nlopt_slsqp},
    {"active_set",           optimize_knitro_active_set},
    {"gradient_descent",     optimize_gd},
    {"custom_gradient_descent", optimize_gd_smartstep}
};

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job", po::value<string>(), "job configuration file")
        ;
    po::positional_options_description p;
    p.add("job", 1);

    po::options_description patternOptions;
    patternOptions.add_options()
        ("pattern,p",    po::value<string>(),                              "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
        ("inflator,i",   po::value<string>()->default_value("isosurface"), "Which inflator to use: Isosurface (default), LpHole, BoundaryPerturbation, Luigi, James")
        ("symmetry",     po::value<string>()->default_value("orthotropic"),"Symmetries to enforce (orthotropic (default), cubic, square, triply_periodic, doubly_periodic)")
        ("vertexThickness,V",                                              "Use vertex thickness instead of edge thickness (3D only)")
        ("cell_size,c",  po::value<double>(),                              "Inflation cell size (3D only)")
        ("params",       po::value<string>(),                              "Initial params (overrides those specified in job file).")
        ("metaParams",   po::value<string>(),                              "Meta params.")
        ("deformedCell", po::value<string>(),                              "Specify the Jacobian of a deformation linearly warping the pattern after meshing (scanline order; in 3D: xx xy xz yx yy yz zx zy zz)")
        ;

    po::options_description meshingOptions;
    meshingOptions.add_options()
        ("meshingOptions,M",      po::value<string>(), "Meshing options configuration file")
        ("max_volume,v",          po::value<double>(), "Maximum element area for remeshing (overrides meshing options)")
        ("hole_segments",         po::value<size_t>(), "Number of segments in hole boundary for LpHoleInflator (default: 64)")
        ("subdivide,S",           po::value<size_t>(), "Number of subdivisions to run for James' inflator (default: 0)")
        ("sub_algorithm,A",       po::value<string>(), "Subdivision algorithm for James' inflator (simple or loop, default: simple)")
        ("inflation_dump_path,D", po::value<string>(), "Dump the inflated geometry immediately after meshing.")
        ("ortho_cell,O",                               "Only mesh and optimize the orthotropic base cell (for orthotropic patterns only")
        ("fullCellInflator",                           "Use the full period cell inflator instead of the reflection-based one")
        ("inflation_graph_radius",po::value<size_t>(), "Number of edges to traverse outward from the symmetry cell when building the inflation graph (defaults to 2)")
        ;

    po::options_description optimizerOptions;
    optimizerOptions.add_options()
        ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
        ("tensor_fit_tolerance", po::value<double>(),                            "tolerance for tensor fitting (stops the moment tensor is reached regardless of other objective terms, works only for slsqp)")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, slsqp, bfgs, lbfgs")
        ;

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
        ("ignoreShear",                                                   "Ignore the shear components in the isotropic tensor fitting")
        ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp global worst case stress measure")
        ("usePthRoot,R",                                                  "Use the true Lp norm for global worst case stress measure (applying pth root)")
        ("WCSWeight",    po::value<double>()->default_value(1.0),         "Weight for the WCS term of the objective")
        ("WCSTarget",    po::value<double>()->default_value(0.0),         "Target for max stress allowed for the structure")
        ("WCSMeasure,w", po::value<string>()->default_value("frobenius"), "Which worst-case stress measure to analyze ('frobenius', 'vonMises')")
        ("FixedMacroStress", po::value<string>(),                         "Use a fixed macroscopic stress tensor instead of the worst-case")
        ("JSWeight",     po::value<double>(),                             "Use the NLLS tensor fitting term with specified weight.")
        ("JIsoWeight",   po::value<double>(),                             "Use the NLLS isotropy fitting term with specified weight (shrinkage side effect).")
        ("JIsoRelWeight",po::value<double>(),                             "Use the relative isotropy fitting term with specified weight.")
        ("proximityRegularizationWeight", po::value<double>(),            "Use a quadratic proximity regularization term with the specified weight.")
        ("proximityRegularizationTarget", po::value<string>(),            "The target parameter values for the proximity regularization term (defaults to initial parameters.)")
        ("proximityRegularizationZeroTarget",                             "Use 0 vector as target parameter of proximity regularization cost function term.")
        ("smoothingRegularizationWeight", po::value<double>(),            "Use a smoothing regularization term for the boundary of the mesh with the specified weight (when using the BoundaryPerturbation inflator).")
        ("LaplacianRegWeight,r", po::value<double>()->default_value(0.0), "Weight for the boundary Laplacian regularization term")
        ("JIsoFixedTarget",                                               "Make JIso just fit to the closest isotropic tensor to the *original* tensor.")
        ("targetVolWeight", po::value<double>()->default_value(0.0),      "Weight for the target volume term of the objective")
        ("targetVol", po::value<double>(),                                "Define target volume")
        ;

    po::options_description constraintOptions;
    constraintOptions.add_options()
        ("TensorFitConstraint,C",  "Enforce homogenized tensor fitting as a nonlinear equality constraint (for optimizers that support this)")
        ("PrintabilityConstraint", "Enforce self-supporting printability constraints as inequality constraints (for optimizers that support this)")
        ;

    po::options_description elasticityOptions;
    elasticityOptions.add_options()
        ("material,m",   po::value<string>(),                    "Base material")
        ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
        ;

    po::options_description gvOptions;
    gvOptions.add_options()
        ("validateGradientComponent", po::value<size_t>(),                   "Run gradient component validation instead of optimization")
        ("nsamples",                  po::value<size_t>()->default_value(5), "Number of gradient component validation samples")
        ("range",                     po::value<string>(),                   "Absolute sweep range (lower:upper)")
        ("rangeRelative",             po::value<double>(),                   "Relative sweep range: current +/- arg * paramBound(compIdx)")
        ("singleIteration",           po::value<size_t>(),                   "Only run the ith iteration of the validation")
        ;

    po::options_description generalOptions;
    generalOptions.add_options()
        ("help,h",                                    "Produce this help message")
        ("output,o",             po::value<string>(), "Output mesh and fields at each iteration")
        ("dumpShapeDerivatives", po::value<string>(), "Dump shape derivative fields for JVol, JS, and WCS")
        ("numProcs",             po::value<size_t>(), "Number of threads to use for TBB parallelism (CGAL mesher, etc.)")
        ("dumpJson",             po::value<string>(), "Dump some information into the specified json file")
        ("hideGradientInformation",                          "Hide gradient informtion (useful when there are too many parameters)")
        ;

    po::options_description visibleOptions;
    visibleOptions.add(patternOptions).add(meshingOptions).add(optimizerOptions)
                  .add(objectiveOptions).add(constraintOptions)
                  .add(elasticityOptions).add(gvOptions).add(generalOptions);

    po::options_description cli_opts;
    cli_opts.add(visibleOptions).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visibleOptions);
    }

    bool fail = false;
    if (vm.count("job") == 0) {
        cout << "Error: must specify input job.opt file" << endl;
        fail = true;
    }

    if (vm.count("pattern") == 0) {
        cout << "Error: must specify pattern mesh" << endl;
        fail = true;
    }

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    if (optimizers.count(vm["solver"].as<string>()) == 0) {
        cout << "Illegal solver specified" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visibleOptions);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
typedef ScalarField<Real> SField;

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, PO::Job<_N> *job)
{
    auto infl_ptr = make_inflator<_N>(args["inflator"].as<string>(),
                                     filterInflatorOptions(args),
                                     job->parameterConstraints);
    auto &inflator = *infl_ptr;

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();
    bool gradientValidationMode = args.count("validateGradientComponent");
    if (!gradientValidationMode) {
        cout << "Target moduli:\t";
        targetC.printOrthotropic(cout);
        cout << endl;
    }

    // cout << "target tensor: " << targetC << endl;

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    // Design a microstructure to achieve the desired material property in a
    // linearly deformed tiling.
    if (args.count("deformedCell")) {
        // Parse Jacobian of the linear deformation.
        Eigen::Matrix<Real, _N, _N> jacobian;

        vector<string> jacobianComponents;
        string jacobianString = args["deformedCell"].as<string>();
        boost::trim(jacobianString);
        boost::split(jacobianComponents, jacobianString, boost::is_any_of("\t "),
                boost::token_compress_on);
        if (jacobianComponents.size() != _N * _N) {
            throw runtime_error("Invalid deformation jacobian");
        }
        for (size_t i = 0; i < _N; ++i) {
            for (size_t j = 0; j < _N; ++j) {
                jacobian(i, j) = stod(jacobianComponents[_N * i + j]);
            }
        }

        Eigen::Matrix<Real, 3, 3> meshing_jacobian = inflator.meshingOptions().jacobian;
        std::cerr << meshing_jacobian << std::endl << std::endl;
        if (!meshing_jacobian.topLeftCorner<_N, _N>().isApprox(jacobian)) {
            std::cerr << "-- Warning: different Jacobian specified in the meshing options and the --deformedCell arguments." << std::endl << std::endl;
        }

        // Apply the appropriate transformation to analyze the deformed cell's
        // properties on the undeformed microstructure geometry.
        // Material properties must be transformed with the inverse Jacobian:
        Eigen::Matrix<Real, _N, _N> jacobianInv(jacobian.inverse());
        mat.setTensor(mat.getTensor().transform(jacobianInv));

        // We compute the deformed tiling's homogenized tensor by homogenizing
        // the undeformed geometry with this transformed material to obtain
        // "Eh" and then computing:
        //      EhDefo = Eh.transform(jacobian);
        // So, Eh = EhDefo.transform(jacobianInv), and we must transform our
        // target tensor accordingly:
        targetC = targetC.transform(jacobianInv);

        cout << "Transformed target tensor:" << endl << targetC << endl << endl;
        cout << "Transformed target moduli:\t";
        targetC.printOrthotropic(cout);
        cout << endl;
        targetS = targetC.inverse();
    }

    auto parseParams = [](string pstring) -> vector<Real> {
        boost::trim(pstring);
        vector<string> tokens;
        boost::split(tokens, pstring, boost::is_any_of("\t "),
                     boost::token_compress_on);
        vector<Real> pvals;
        for (string &s : tokens) pvals.push_back(std::stod(s));
        return pvals;
    };

    auto parseMetaParams = [](string pstring) -> vector<string> {
        boost::trim(pstring);
        vector<string> tokens;
        boost::split(tokens, pstring, boost::is_any_of("\t "),
                     boost::token_compress_on);
        return tokens;
    };

    // If requested, override the initial parameters set in the job file
    if (args.count("params"))
        job->initialParams = parseParams(args["params"].as<string>());

    if (args.count("metaParams"))
        job->metaParams = parseMetaParams(args["metaParams"].as<string>());

    SField params = job->validatedInitialParams(inflator);

    PO::BoundConstraints bdcs(inflator, job->radiusBounds, job->translationBounds, job->blendingBounds, job->metaBounds,
                              job->custom1Bounds, job->custom2Bounds, job->custom3Bounds, job->custom4Bounds,
                              job->custom5Bounds, job->custom6Bounds, job->custom7Bounds, job->custom8Bounds,
                              job->varLowerBounds, job->varUpperBounds);

    // TODO: Laplacian regularization term (probably only needed for boundary
    // perturbation version.

    using WCSTermConfig          = PO::ObjectiveTerms::IFConfigWorstCaseStress<Simulator>;
    using TensorFitTermConfig    = PO::ObjectiveTerms::IFConfigTensorFit<Simulator>;
    using IsotropyFitConfig      = PO::ObjectiveTerms::IFConfigIsotropyFit<Simulator>;
    using IsoFitRelConfig        = PO::ObjectiveTerms::IFConfigIsotropyFitRel<Simulator>;
    using PRegTermConfig         = PO::ObjectiveTerms::IFConfigProximityRegularization;
    using SRegTermConfig         = PO::ObjectiveTerms::IFConfigPeriodicSmoothingRegularization<Simulator>;
    using TargetVolumeTermConfig = PO::ObjectiveTerms::IFConfigTargetVolume<Simulator>;
    using TFConstraintConfig  = PO::   Constraints::IFConfigTensorFit<Simulator>;
    using  PConstraintConfig  = PO::   Constraints::IFConfigPrintability<Simulator>;

    bool outputGradientInformation = !args.count("hideGradientInformation");

    auto ifactory = PO::make_iterate_factory<PO::Iterate<Simulator>,
         WCSTermConfig,
         TensorFitTermConfig,
         IsotropyFitConfig,
         IsoFitRelConfig,
         PRegTermConfig,
         SRegTermConfig,
         TFConstraintConfig,
         TargetVolumeTermConfig,
         PConstraintConfig>(inflator, bdcs, outputGradientInformation);

    ////////////////////////////////////////////////////////////////////////////
    // Configure the objective terms
    ////////////////////////////////////////////////////////////////////////////
    ifactory->WCSTermConfig           ::enabled = args["WCSWeight"].as<double>() != 0;
    ifactory->TensorFitTermConfig     ::enabled = args.count("JSWeight");
    ifactory->IsotropyFitConfig       ::enabled = args.count("JIsoWeight");
    ifactory->IsoFitRelConfig         ::enabled = args.count("JIsoRelWeight");
    ifactory->PRegTermConfig          ::enabled = args.count("proximityRegularizationWeight");
    ifactory->SRegTermConfig          ::enabled = args.count("smoothingRegularizationWeight");
    ifactory->TFConstraintConfig      ::enabled = false;
    ifactory->PConstraintConfig       ::enabled = false;
    ifactory->TargetVolumeTermConfig  ::enabled = args["targetVolWeight"].as<double>() > 0.0;

    // Configure WCS Objective
    // By default, an "Lp norm" objective is really the p^th power of the Lp norm.
    // To use the true "Lp norm", globalObjectiveRoot must be set to
    // 2.0 * globalObjectivePNorm (since pointwise WCS is already squared (e.g. Frobenius) norm)
    // Set target as squared value (since pointwise WCS is already squared (e.g. Frobenius) norm)
    ifactory->WCSTermConfig::weight = args["WCSWeight"].as<double>();
    double target = args["WCSTarget"].as<double>();
    ifactory->WCSTermConfig::target = target * target;
    Real pnorm = args["pnorm"].as<double>();
    ifactory->WCSTermConfig::globalObjectivePNorm = pnorm;
    ifactory->WCSTermConfig::globalObjectiveRoot  = args.count("usePthRoot") ? 2.0 * pnorm : 1.0;

    // Configure the pointwise worst-case stress measure
    auto measure = args["WCSMeasure"].as<string>();
    boost::algorithm::to_lower(measure);
    ifactory->WCSTermConfig::measure = measure;

    if (args.count("FixedMacroStress")) {
        auto macroStressStr = args["FixedMacroStress"].as<string>();
        vector<string> stressComponents;
        boost::trim(macroStressStr);
        boost::split(stressComponents, macroStressStr, boost::is_any_of("\t ,"),
                     boost::token_compress_on);
        if (stressComponents.size() != flatLen(_N))
            throw runtime_error("Invalid FixedMacroStress tensor");
        ifactory->WCSTermConfig::macroLoad =
            Future::make_unique<typename Simulator::SMatrix>();
        for (size_t i = 0; i < flatLen(_N); ++i)
            (*ifactory->WCSTermConfig::macroLoad)[i] = std::stod(stressComponents[i]);
    }

    if (args.count("JSWeight")) {
        ifactory->TensorFitTermConfig::weight  = args["JSWeight"].as<double>();
        ifactory->TensorFitTermConfig::targetS = targetS;
        ifactory->TensorFitTermConfig::ignoreShear = args.count("ignoreShear");
        if (ifactory->TensorFitTermConfig::ignoreShear) cout << "Ignoring shear components" << endl;
    }

    if (args.count("JIsoRelWeight")) {
        ifactory->IsoFitRelConfig::weight = args["JIsoRelWeight"].as<double>();
    }

    if (args.count("JIsoWeight")) {
        ifactory->IsotropyFitConfig::weight = args["JIsoWeight"].as<double>();
        ifactory->IsotropyFitConfig::useFixedTarget = args.count("JIsoFixedTarget");
    }

    if (args.count("TensorFitConstraint")) {
        ifactory->TFConstraintConfig::enabled             = true;
        ifactory->TFConstraintConfig::targetS             = targetS;
        ifactory->TFConstraintConfig::ignoreShear         = args.count("ignoreShear");
        ifactory->TFConstraintConfig::orthotropicSymmetry = inflator.hasOrthotropicSymmetry();
    }

    if (args.count("PrintabilityConstraint")) {
        ifactory->PConstraintConfig::enabled = true;
    }

    if (args.count("proximityRegularizationWeight")) {
        ifactory->PRegTermConfig::enabled = true;
        ifactory->PRegTermConfig::weight = args["proximityRegularizationWeight"].as<double>();

        // if parameter of provided, then use them. Otherwise, assume we should set initial parameters as target
        if (job->targetParams.size())
            ifactory->PRegTermConfig::targetParams = job->targetParams;
        else
            ifactory->PRegTermConfig::targetParams = job->validatedInitialParams(inflator);


        if (args.count("proximityRegularizationTarget")) {
            ifactory->PRegTermConfig::targetParams = parseParams(args["proximityRegularizationTarget"].as<string>());
            if (ifactory->PRegTermConfig::targetParams.size() != params.domainSize())
                throw runtime_error("Invalid proximity regularization target parameter count");
        }
        if (args.count("proximityRegularizationZeroTarget")) {
            vector<Real> zeroTargetParams(job->numParams(), 0.0);
            ifactory->PRegTermConfig::IFConfigProximityRegularization::targetParams = zeroTargetParams;
        }
    }

    if (args.count("smoothingRegularizationWeight")) {
        ifactory->SRegTermConfig::enabled = true;
        ifactory->SRegTermConfig::weight = args["smoothingRegularizationWeight"].as<double>();
    }

    if (args["targetVolWeight"].as<double>() > 0.0) {
        if (args.count("targetVol") > 0.0) {
            ifactory->TargetVolumeTermConfig::enabled = true;
            ifactory->TargetVolumeTermConfig::weight = args["targetVolWeight"].as<double>();
            ifactory->TargetVolumeTermConfig::targetVolume = args["targetVol"].as<double>();
        }
        else if (job->targetVolume) {
            ifactory->TargetVolumeTermConfig::enabled = true;
            ifactory->TargetVolumeTermConfig::weight = args["targetVolWeight"].as<double>();
            ifactory->TargetVolumeTermConfig::targetVolume = *(job->targetVolume);
        }
    }

    auto imanager = PO::make_iterate_manager(std::move(ifactory));

    ////////////////////////////////////////////////////////////////////////////
    // Gradient component validation, if requested, bypasses optimization
    ////////////////////////////////////////////////////////////////////////////
    if (gradientValidationMode) {
        size_t compIdx = args["validateGradientComponent"].as<size_t>();
        if (compIdx >= params.domainSize()) throw runtime_error("Gradient component index out of bounds");
        if (args.count("range") == args.count("rangeRelative"))
            throw runtime_error("Either range or rangeRelative must be specified (not both)");

        if (!bdcs.hasLowerBound.at(compIdx) || !bdcs.hasUpperBound.at(compIdx))
            throw runtime_error("Swept parameters must be bounded");

        Real prlb = bdcs.lowerBound[compIdx], prub = bdcs.upperBound[compIdx];
        Real lb, ub;

        if (args.count("range")) {
            auto rangeStr = args["range"].as<string>();
            vector<string> rangeComponents;
            boost::trim(rangeStr), boost::split(rangeComponents, rangeStr, boost::is_any_of(":"));
            if (rangeComponents.size() != 2) throw runtime_error("Invalid range; expected lower:upper");
            lb = stod(rangeComponents[0]), ub = stod(rangeComponents[1]);
        }
        else {
            Real rr = args["rangeRelative"].as<double>();
            Real prSize = prub - prlb;
            lb = params[compIdx] - rr * prSize, ub = params[compIdx] + rr * prSize;
        }

        if ((lb < prlb) || (ub > prub)) {
            std::cerr << "WARNING: Specified sweep range of " << lb << ":" << ub
                << " outside parameter range of " << prlb << ":" << prub << std::endl;
        }

        params[compIdx] = lb; // make dummy iterate (actually 0th iterate, will be reused)
        inflator.meshingOptions().debugSVelPath = "svels.msh";
        cout << "it\tparam\tJFull\tgradp JFull";
        {
            const auto &it = imanager->get(params.size(), params.data());
            for (const auto &etermptr : it.evaluatedObjectiveTerms())
                cout << "\t" << etermptr->name << "\tgradp " << etermptr->name;
        }
        cout << endl;

        const size_t nsamples = args["nsamples"].as<size_t>();
        for (size_t i = 0; i < nsamples; ++i) {
            if (args.count("singleIteration")) i = args["singleIteration"].as<size_t>();
            params[compIdx] = lb + ((nsamples == 1) ? 0.0 : (ub - lb) * (double(i) / (nsamples - 1)));
            auto &it = imanager->get(params.size(), params.data());
            cout << i << "\t" << params[compIdx] << "\t" << it.evaluate() << "\t" << it.gradp()[compIdx];
            for (const auto &etermptr : it.evaluatedObjectiveTerms())
                cout << "\t" << etermptr->value() << "\t" << etermptr->gradp[compIdx];
            cout << endl;

            if (args.count("output")) it.writeMeshAndFields(args["output"].as<string>() + "_" + std::to_string(i) + ".msh");
            if (args.count("singleIteration")) break;
        }

        // BENCHMARK_REPORT();
        return;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Run the optimizer
    ////////////////////////////////////////////////////////////////////////////
    string solver = args["solver"].as<string>(), output;
    if (args.count("output")) output = args["output"].as<string>();

    PO::OptimizerConfig oconfig;
    if (args.count("nIters")) oconfig.niters = args["nIters"].as<size_t>();
    oconfig.gd_step = args["step"].as<double>();

    if (args.count("tensor_fit_tolerance"))
        oconfig.tensor_fit_tolerance = args["tensor_fit_tolerance"].as<double>();

    if (solver == "lbfgs") oconfig.lbfgs_memory = 10;
    optimizers.at(solver)(params, bdcs, *imanager, oconfig, output);

    ////////////////////////////////////////////////////////////////////////////
    // Extract and process the result.
    ////////////////////////////////////////////////////////////////////////////
    std::vector<Real> result(params.domainSize());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = params[i];

    json output_json;
    if (inflator.isParametric()) {
        output_json["final_p"] = json::array();
        cout << "Final p:";
        for (size_t i = 0; i < result.size(); ++i) {
            cout << "\t" << result[i];
            output_json["final_p"].push_back(result[i]);
        }
        cout << endl;
    }

    if (args.count("dumpJson")) {
        std::ofstream out(args["dumpJson"].as<string>());
        out << output_json;
    }

    BENCHMARK_REPORT_NO_MESSAGES();
}

int main(int argc, const char *argv[]) {
    po::variables_map args = parseCmdLine(argc, argv);

#if MICRO_WITH_TBB
    size_t np = tbb::task_scheduler_init::default_num_threads();
    if (args.count("numProcs")) {
        size_t manualNP = args["numProcs"].as<size_t>();
        if (manualNP > np)
            std::cerr << "WARNING: specifying more than the default number of TBB threads." << std::endl;
        np = manualNP;
    }
    tbb::task_scheduler_init init(np);
#else
    if (args.count("numProcs"))
        std::cerr << "WARNING: parallelism disabled; numProcs argument ignored." << std::endl;
#endif

    cout << setprecision(16);
    auto job = PO::parseJobFile(args["job"].as<string>());

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<PO::Job<2> *>(job.get())) {
        if (deg == 1) execute<2, 1>(args, job2D);
        if (deg == 2) execute<2, 2>(args, job2D);
    }
    else if (auto job3D = dynamic_cast<PO::Job<3> *>(job.get())) {
        if (deg == 1) execute<3, 1>(args, job3D);
        if (deg == 2) execute<3, 2>(args, job3D);
    }
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
