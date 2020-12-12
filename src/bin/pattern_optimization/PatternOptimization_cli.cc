////////////////////////////////////////////////////////////////////////////////
// PatternOptimization_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evolves the microstructure parameters to bring the structure's
//      homogenized elasticity tensor closer to a target tensor.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/12/2014 01:15:28
////////////////////////////////////////////////////////////////////////////////
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/Materials.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/GlobalBenchmark.hh>

#include <inflators/Inflator.hh>
#include <inflators/MakeInflator.hh>

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <map>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>

#include <optimizers/wrappers/ceres.hh>
#include <optimizers/wrappers/dlib.hh>
#include <optimizers/wrappers/gradient_descent.hh>
#include <optimizers/wrappers/nlopt.hh>
#include <optimizers/wrappers/knitro.hh>

#include <optimizers/BoundConstraints.hh>
#include <pattern_optimization/PatternOptimizationIterate.hh>
#include <pattern_optimization/PatternOptimizationJob.hh>
#include <pattern_optimization/IterateFactory.hh>
#include <pattern_optimization/IterateManager.hh>

#include <optimizers/OptimizerConfig.hh>
#include <pattern_optimization/PatternOptimizationConfig.hh>
#include <pattern_optimization/objective_terms/TensorFit.hh>
#include <pattern_optimization/objective_terms/ProximityRegularization.hh>
#include <pattern_optimization/constraints/TensorFit.hh>

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

using OptimizerMap =
    map<string, std::function<void(ScalarField<Real> &, const BoundConstraints &,
            IterateManagerBase &, const OptimizerConfig &, const string &)>>;
OptimizerMap optimizers = {
    {"levenberg_marquardt",  optimize_ceres_lm},
    {"dogleg",               optimize_ceres_dogleg},
    {"bfgs",                 optimize_dlib_bfgs},
    {"lbfgs",                optimize_dlib_bfgs},
    {"slsqp",                optimize_nlopt_slsqp},
    {"active_set",           optimize_knitro_active_set},
    {"gradient_descent",     optimize_gd}
};

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: PatternOptimization_cli [options] job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job", po::value<string>(), "job configuration file")
        ;
    po::positional_options_description p;
    p.add("job", 1);

    po::options_description misc_opts;
    misc_opts.add_options()("help",        "Produce this help message")
        ("inflator,i",   po::value<string>()->default_value("Isosurface"),       "inflator to use (defaults to Isosurface)")
        ("symmetry",     po::value<string>()->default_value("orthotropic"),      "Symmetries to enforce (orthotropic (default), cubic, square, triply_periodic, doubly_periodic)")
        ("pattern,p",    po::value<string>(),                                    "Pattern wire mesh (.obj|wire)")
        ("material,m",   po::value<string>(),                                    "base material")
        ("degree,d",     po::value<size_t>()->default_value(2),                  "FEM Degree")
        ("output,o",     po::value<string>()->default_value(""),                 "output .js mesh + fields at each iteration")
        ("cell_size,c",  po::value<double>(),                                    "Inflation cell size (James' inflator only. Default: 5mm)")
        ("vertexThickness,V",                                                    "Use vertex thickness instead of edge thickness (3D only)")
        ("proximityRegularizationWeight", po::value<double>(),                   "Use a quadratic proximity regularization term with the specified weight.")
        ("proximityRegularizationZeroTarget",                                    "Use 0 vector as target parameter of proximity regularization cost function term.")
        ;

    po::options_description constraintOptions;
    constraintOptions.add_options()
        ("TensorFitConstraint,C", "Enforce homogenized tensor fitting as a nonlinear equality constraint (for optimizers that support this)")
        ;

    po::options_description optimizer_options;
    optimizer_options.add_options()
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs, slsqp, levenberg_marquardt")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient descent step size")
        ("nIters,n",     po::value<size_t>(),                                    "number of iterations (infinite by default)")
        ("tensor_fit_tolerance", po::value<double>(),                            "tolerance for tensor fitting (stops the moment tensor is reached regardless of other objective terms, works only for slsqp)")
        ("ignoreShear",                                                          "Ignore the shear components in the isotropic tensor fitting")
        ;

    po::options_description meshingOptions;
    meshingOptions.add_options()
        ("meshingOptions,M", po::value<string>(),  "Meshing options configuration file")
        ("max_volume,v",     po::value<double>(),  "Maximum element area for remeshing (overrides meshing options)")
        ("subdivide,S",  po::value<size_t>(),      "Number of subdivisions to run for James' inflator (default: 0)")
        ("sub_algorithm,A", po::value<string>(),   "Subdivision algorithm for James' inflator (simple or loop, default: simple)")
        ("ortho_cell,O",                           "Only mesh and optimize the orthotropic base cell (for orthotropic patterns only")
        ("fullCellInflator",                       "Use the full period cell inflator instead of the reflection-based one")
        ("inflation_graph_radius",po::value<size_t>(), "Number of edges to traverse outward from the symmetry cell when building the inflation graph (defaults to 2)")
        ;

    po::options_description visible_opts;
    visible_opts.add(misc_opts).add(constraintOptions).add(meshingOptions).add(optimizer_options);

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
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
        usage(fail, visible_opts);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
typedef ScalarField<Real> SField;

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, const Job<_N> *job)
{
    auto infl_ptr = make_inflator<_N>(args["inflator"].as<string>(),
                                     filterInflatorOptions(args),
                                     job->parameterConstraints);

    Inflator<_N> &inflator = *infl_ptr;

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    cout << "target tensor: " << targetC << endl;

    SField params = job->validatedInitialParams(inflator);

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    BoundConstraints bdcs(inflator, job->radiusBounds, job->translationBounds, job->blendingBounds, job->metaBounds,
                          job->custom1Bounds, job->custom2Bounds, job->custom3Bounds, job->custom4Bounds,
                          job->custom5Bounds, job->custom6Bounds, job->custom7Bounds, job->custom8Bounds,
                          job->varLowerBounds, job->varUpperBounds);

    using TFConstraintConfig  = Constraints::IFConfigTensorFit<Simulator>;
    using TensorFitTermConfig = ObjectiveTerms::IFConfigTensorFit<Simulator>;
    using Iterate = Iterate<Simulator>;
    auto ifactory = make_iterate_factory<Iterate,
                                         TensorFitTermConfig,
                                         ObjectiveTerms::IFConfigProximityRegularization,
                                         TFConstraintConfig>(inflator, bdcs, true);

    bool ignoreShear = args.count("ignoreShear");
    ifactory->TensorFitTermConfig::ignoreShear = ignoreShear;
    if (ignoreShear) cout << "Ignoring shear components" << endl;
    ifactory->ObjectiveTerms::IFConfigProximityRegularization::enabled = false;
    if (args.count("proximityRegularizationWeight")) {
        ifactory->ObjectiveTerms::IFConfigProximityRegularization::enabled      = true;
        ifactory->ObjectiveTerms::IFConfigProximityRegularization::targetParams = job->validatedInitialParams(inflator);
        if (args.count("proximityRegularizationZeroTarget")) {
            // Split up params.
            vector<Real> targetParams(job->numParams(), 0.0);
            ifactory->ObjectiveTerms::IFConfigProximityRegularization::targetParams = targetParams;
        }
        else {
        }

        ifactory->ObjectiveTerms::IFConfigProximityRegularization::weight       = args["proximityRegularizationWeight"].as<double>();
    }

    ifactory->TensorFitTermConfig::targetS = targetS;

    ifactory->TFConstraintConfig ::enabled = false;
    if (args.count("TensorFitConstraint")) {
        ifactory->TFConstraintConfig::enabled             = true;
        ifactory->TFConstraintConfig::targetS             = targetS;
        ifactory->TFConstraintConfig::ignoreShear         = args.count("ignoreShear");
        ifactory->TFConstraintConfig::orthotropicSymmetry = inflator.hasOrthotropicSymmetry();
    }

    auto imanager = make_iterate_manager(std::move(ifactory));

    string solver = args["solver"].as<string>(),
           output = args["output"].as<string>();
    OptimizerConfig oconfig;
    if (args.count("nIters")) oconfig.niters = args["nIters"].as<size_t>();
    oconfig.gd_step = args["step"].as<double>();

    if (args.count("tensor_fit_tolerance"))
        oconfig.tensor_fit_tolerance = args["tensor_fit_tolerance"].as<double>();

    if (solver == "lbfgs") oconfig.lbfgs_memory = 10;
    optimizers.at(solver)(params, bdcs, *imanager, oconfig, output);

    if (inflator.isParametric()) {
        cout << "Final p:";
        vector<Real> result(params.domainSize());
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] = params[i];
            cout << "\t" << params[i];
        }
        cout << endl;
    }

    BENCHMARK_REPORT_NO_MESSAGES();
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    cout << setprecision(16);

    po::variables_map args = parseCmdLine(argc, argv);

    auto job = parseJobFile(args["job"].as<string>());

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<Job<2> *>(job.get())) {
        if (deg == 1) execute<2, 1>(args, job2D);
        if (deg == 2) execute<2, 2>(args, job2D);
    }
    else if (auto job3D = dynamic_cast<Job<3> *>(job.get())) {
        if (deg == 1) execute<3, 1>(args, job3D);
        if (deg == 2) execute<3, 2>(args, job3D);
    }
    else throw runtime_error("Invalid job file.");

    return 0;
}
