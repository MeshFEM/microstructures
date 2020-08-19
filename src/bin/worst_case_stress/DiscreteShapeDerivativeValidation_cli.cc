////////////////////////////////////////////////////////////////////////////////
// DiscreteShapeDerivativeValidation_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validate the discrete shape derivative of periodic homogenization and
//      worst case stress quantities.
//
//      Boundary vertices are offset in the normal direction to create a
//      perturbed mesh, and post- and pre-perturbation quantities are
//      subtracted to compute forward/centered difference (material)
//      derivatives.
//
//      The BoundaryPerturbationInflator is used to ease the creation of a
//      perturbed mesh.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  03/16/2016 18:14:07
////////////////////////////////////////////////////////////////////////////////

#include "WorstCaseStress.hh"

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/Materials.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/Laplacian.hh>
#include <inflators/wrappers/BoundaryPerturbationInflator.hh>
#include <isosurface_inflator/ShapeVelocityInterpolator.hh>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iomanip>

namespace po = boost::program_options;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: DiscreteShapeDerivativeValidation_cli [options] mesh.msh" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("mesh", po::value<string>(), "input mesh")
        ;
    po::positional_options_description p;
    p.add("mesh", 1);

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
        ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp global worst case stress measure")
        ("usePthRoot,R",                                                  "Use the true Lp norm for global worst case stress measure (applying pth root)")
        ;

    po::options_description elasticityOptions;
    elasticityOptions.add_options()
        ("material,m",   po::value<string>(),                    "Base material")
        ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
        ;

    po::options_description generalOptions;
    generalOptions.add_options()
        ("help,h",                                               "Produce this help message")
        ("output,o",     po::value<string>(),                    "Output the Lagrangian derivatives computed by forward difference and the discrete shape derivative.")
        ("fullDegreeFieldOutput,D",                              "Output full-degree nodal fields (don't do piecewise linear subsample)")
        ("perturbationAmplitude,a", po::value<double>()->default_value(0.01), "Amplitude of boundary perturbation")
        ("perturbationFrequency,f", po::value<double>()->default_value(1.0),  "Frequency of boundary perturbation")
        ;

    po::options_description visibleOptions;
    visibleOptions.add(objectiveOptions).add(elasticityOptions).add(generalOptions);

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

    if (vm.count("mesh") == 0) {
        cout << "Error: must specify input mesh" << endl;
        fail = true;
    }

    if (vm.count("output") == 0) {
        cout << "Error: must specify output mesh" << endl;
        fail = true;
    }

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visibleOptions);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args,
             const std::vector<MeshIO::IOVertex>  &vertices,
             const std::vector<MeshIO::IOElement> &elements)
{
    using Mesh      = typename LinearElasticity::Mesh<_N, _FEMDegree, HMG>;
    using Simulator = typename LinearElasticity::Simulator<Mesh>;
    using ETensor   = typename Simulator::ETensor;
    using VField    = typename Simulator::VField;

    // Original mesh and simulator
    BoundaryPerturbationInflator<_N> bpi(vertices, elements);
    vector<Real> perturbParams(bpi.numParameters());
    bpi.inflate(perturbParams);
    Simulator sim(bpi.elements(), bpi.vertices());

    // Perturb mesh
    {
        auto &mesh = sim.mesh();
        VField perturbation(sim.mesh().numBoundaryVertices());
        Real A = args["perturbationAmplitude"].as<Real>();
        Real f = args["perturbationFrequency"].as<Real>();
        auto normals = bpi.boundaryVertexNormals();
        for (auto bv : mesh.boundaryVertices()) {
            auto pt = bv.node().volumeNode()->p;
            Real a = A * cos(M_PI * f * pt[0]) * cos(M_PI * f * pt[1]);
            perturbation(bv.index()) = a * normals(bv.index());
        }
        bpi.paramsFromBoundaryVField(perturbation).getFlattened(perturbParams);
    }

    // Perturbed meshes and simulators
    bpi.inflate(perturbParams);
    Simulator perturbed_sim(bpi.elements(), bpi.vertices());
    for (Real &p : perturbParams) p *= -1.0;
    bpi.inflate(perturbParams);
    Simulator neg_perturbed_sim(bpi.elements(), bpi.vertices());

    // Determine change in each vertex's position.
    VField delta_p(sim.mesh().numVertices());
    for (auto v : sim.mesh().vertices()) {
        delta_p(v.index()) = perturbed_sim.mesh().vertex(v.index()).node()->p;
        delta_p(v.index()) -= v.node()->p;
    }

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    std::vector<VField> w;
    PeriodicHomogenization::solveCellProblems(w, sim);
    auto delta_w = PeriodicHomogenization::deltaFluctuationDisplacements(sim, w, delta_p);

    std::vector<VField> perturbed_w, neg_perturbed_w;
    PeriodicHomogenization::solveCellProblems(perturbed_w, perturbed_sim);
    PeriodicHomogenization::solveCellProblems(neg_perturbed_w, neg_perturbed_sim);

    std::vector<VField> delta_w_forward_diff = perturbed_w;
    std::vector<VField> delta_w_centered_diff = perturbed_w;
    for (size_t ij = 0; ij < w.size(); ++ij) {
        delta_w_forward_diff[ij] -= w[ij];
        delta_w_centered_diff[ij] -= neg_perturbed_w[ij];
        delta_w_centered_diff[ij] *= 0.5;
    }

    string output = args["output"].as<string>();
    bool linearSubsampleFields = args.count("fullDegreeFieldOutput") == 0;
    MSHFieldWriter writer(output, sim.mesh(), linearSubsampleFields);

    for (size_t ij = 0; ij < w.size(); ++ij) {
        writer.addField("w " + std::to_string(ij), w[ij]);
        writer.addField("delta w " + std::to_string(ij), delta_w[ij]);
        writer.addField("forward difference delta w "  + std::to_string(ij), delta_w_forward_diff[ij]);
        writer.addField("centered difference delta w " + std::to_string(ij), delta_w_centered_diff[ij]);
    }

    for (size_t ij = 0; ij < w.size(); ++ij) {
        auto origStrain = sim.averageStrainField(w[ij]);
        auto forwardDiffStrain   = perturbed_sim.averageStrainField(perturbed_w[ij]);
        auto centeredDiffStrain  = forwardDiffStrain;
        forwardDiffStrain  -= origStrain;
        centeredDiffStrain -= neg_perturbed_sim.averageStrainField(neg_perturbed_w[ij]);
        centeredDiffStrain *= 0.5;

        writer.addField("strain w "       + std::to_string(ij), origStrain);
        writer.addField("delta strain w " + std::to_string(ij), sim.deltaAverageStrainField(w[ij], delta_w[ij], delta_p));
        writer.addField("forward difference delta strain w "  + std::to_string(ij),  forwardDiffStrain);
        writer.addField("centered difference delta strain w " + std::to_string(ij), centeredDiffStrain);
    }

    MSHFieldWriter perturbed_writer(output + ".perturbed.msh", perturbed_sim.mesh(), linearSubsampleFields);
    for (size_t ij = 0; ij < w.size(); ++ij) {
        perturbed_writer.addField("w "+ std::to_string(ij), perturbed_w[ij]);
        perturbed_writer.addField("strain w " + std::to_string(ij), perturbed_sim.averageStrainField(perturbed_w[ij]));
    }

    MSHFieldWriter neg_perturbed_writer(output + ".neg_perturbed.msh", neg_perturbed_sim.mesh(), linearSubsampleFields);
    for (size_t ij = 0; ij < w.size(); ++ij) {
        neg_perturbed_writer.addField("w "+ std::to_string(ij), neg_perturbed_w[ij]);
        neg_perturbed_writer.addField("strain w " + std::to_string(ij), perturbed_sim.averageStrainField(neg_perturbed_w[ij]));
    }

    // Validate homogenized elasticity tensor shape derivative
    auto deltaCh = PeriodicHomogenization::deltaHomogenizedElasticityTensor(sim, w, delta_p);
    auto Ch = PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(w, sim);
    auto perturbed_Ch = PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(perturbed_w, perturbed_sim);

    auto fdDeltaCh = perturbed_Ch - Ch;
    cout << "deltaCh: " << endl << deltaCh << endl;
    cout << "fdDeltaCh: " << endl << fdDeltaCh << endl;
    cout << "Relative error: " << sqrt((deltaCh - fdDeltaCh).frobeniusNormSq() / fdDeltaCh.frobeniusNormSq()) << endl;
    cout << endl;

    auto neg_perturbed_Ch = PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(neg_perturbed_w, neg_perturbed_sim);
    auto cdDeltaCh = perturbed_Ch - neg_perturbed_Ch;
    cdDeltaCh *= 0.5;
    cout << "cdDeltaCh: " << endl << cdDeltaCh << endl;
    cout << "Relative error: " << sqrt((deltaCh - cdDeltaCh).frobeniusNormSq() / cdDeltaCh.frobeniusNormSq()) << endl;
    cout << endl;


    auto buildWCSObjective = [&](const Simulator &sim_, const ETensor &Ch_, vector<VField> w_) {
        PthRootObjective<IntegratedWorstCaseObjective<_N, WCStressIntegrandLp>> objective;

        objective.integrand.p = args["pnorm"].as<double>();
        objective.p = args.count("usePthRoot") ? 2.0 * args["pnorm"].as<double>() : 1.0;

        objective.setPointwiseWCS(sim_.mesh(),
            worstCaseFrobeniusStress(mat.getTensor(), Ch_.inverse(),
                PeriodicHomogenization::macroStrainToMicroStrainTensors(w_, sim_)));
        return objective;
    };

    auto origWCSObjective          = buildWCSObjective(sim, Ch, w);
    auto perturbedWCSObjective     = buildWCSObjective(perturbed_sim, perturbed_Ch, perturbed_w);
    auto neg_perturbedWCSObjective = buildWCSObjective(neg_perturbed_sim, neg_perturbed_Ch, neg_perturbed_w);

    cout << "WCS:\t" << origWCSObjective.evaluate() << endl;
    cout << "Perturbed WCS:\t" << perturbedWCSObjective.evaluate() << endl;
    cout << "Neg Perturbed WCS:\t" << neg_perturbedWCSObjective.evaluate() << endl;

    cout << "Forward  difference WCS:\t" << perturbedWCSObjective.evaluate() - origWCSObjective.evaluate() << endl;
    cout << "Centered difference WCS:\t" << 0.5 * (perturbedWCSObjective.evaluate() - neg_perturbedWCSObjective.evaluate()) << endl;
    cout << "Discrete shape derivative WCS:\t"    << origWCSObjective.deltaJ(sim, w, delta_p) << endl;

#if 0 // worst-case stress code no longer implements the inaccurate continuous shape derivative formula.
    // Compute shape derivative using continuous formulation
    using NSVI = Interpolant<Real, _N - 1, 1>;
    std::vector<NSVI> delta_p_nsv; delta_p_nsv.reserve(sim.mesh().numBoundaryElements());
    NSVI nsv;
    for (auto be : sim.mesh().boundaryElements()) {
        // Compute boundary element's linear normal velocity under delta_p
        for (size_t i = 0; i < be.numVertices(); ++i)
            nsv[i] = be->normal().dot(delta_p(be.vertex(i).volumeVertex().index()));
        delta_p_nsv.push_back(nsv);
    }
    cout << "Continuous shape derivative WCS:\t" << origWCSObjective.directDerivative(sim, w, delta_p_nsv) << endl;
#endif

    // Compute shape derivative ignoring the motion of internal vertices
    VField delta_p_bdryonly(sim.mesh().numVertices());
    delta_p_bdryonly.clear();
    for (auto bv : sim.mesh().boundaryVertices()) {
        size_t vi = bv.volumeVertex().index();
        delta_p_bdryonly(vi) = delta_p(vi);
    }

    cout << "Boundary-perturbation-only discrete shape derivative WCS:\t" << origWCSObjective.deltaJ(sim, w, delta_p_bdryonly) << endl;

    // Compute shape derivative assuming smooth motion of internal vertices
    // (Solve for interior vertices using deg 1 Laplace equation)
    const auto &mesh = sim.mesh();
    VField laplacian_delta_p(delta_p.domainSize());
    auto L = Laplacian::construct<1>(mesh);
    std::vector<size_t> bdryVertices;
    std::vector<Real> bdry_delta_p;
    bdryVertices.reserve(mesh.numBoundaryVertices());
    bdry_delta_p.reserve(mesh.numBoundaryVertices());
    for (auto bv : mesh.boundaryVertices())
        bdryVertices.push_back(bv.volumeVertex().index());
    for (size_t c = 0; c < _N; ++c) {
        // TODO: avoid refactorization of SPSDSystem when fixed var *values*
        // are changed. (keep track of variables contributing to RHS)
        SPSDSystem<Real> Lsys(L);
        bdry_delta_p.clear();
        for (auto bv : mesh.boundaryVertices())
            bdry_delta_p.push_back(delta_p(bv.volumeVertex().index())[c]);
        Lsys.fixVariables(bdryVertices, bdry_delta_p);
        std::vector<Real> x;
        Lsys.solve(std::vector<Real>(mesh.numVertices()), x);
        for (size_t vi = 0; vi < mesh.numVertices(); ++vi)
            laplacian_delta_p(vi)[c] = x.at(vi);
    }
    writer.addField("delta_p", delta_p);
    writer.addField("laplacian_delta_p", laplacian_delta_p);

    cout << "Laplacian delta p discrete shape derivative WCS:\t" << origWCSObjective.deltaJ(sim, w, laplacian_delta_p) << endl;

    VField bdry_svel(mesh.numBoundaryVertices());
    for (auto bv : mesh.boundaryVertices())
        bdry_svel(bv.index()) = delta_p(bv.volumeVertex().index());

    ShapeVelocityInterpolator interpolator(sim);
    cout << "Periodic laplacian delta p discrete shape derivative WCS:\t" << origWCSObjective.deltaJ(sim, w, interpolator.interpolate(sim, bdry_svel)) << endl;

    auto cellOps = constructBaseCellOps(BaseCellType::TriplyPeriodic, sim);
    OneForm<Real, _N> dJ = origWCSObjective.adjointDeltaJ(*cellOps);
    cout << "Adjoint discrete shape derivative WCS (volume):\t" << dJ[interpolator.interpolate(sim, bdry_svel)] << endl;
    OneForm<Real, _N> dJbdry = interpolator.adjoint(sim, dJ);
    cout << "Adjoint discrete shape derivative WCS (boundary):\t" << dJbdry[bdry_svel] << endl;
}

int main(int argc, const char *argv[])
{
    po::variables_map args = parseCmdLine(argc, argv);

    vector<MeshIO::IOVertex>  inVertices;
    vector<MeshIO::IOElement> inElements;
    string meshPath = args["mesh"].as<string>();
    auto type = load(meshPath, inVertices, inElements, MeshIO::FMT_GUESS,
                     MeshIO::MESH_GUESS);

    // Infer dimension from mesh type.
    size_t dim;
    if      (type == MeshIO::MESH_TET) dim = 3;
    else if (type == MeshIO::MESH_TRI) dim = 2;
    else    throw std::runtime_error("Mesh must be triangle or tet.");

    // Look up and run appropriate homogenizer instantiation.
    int deg = args["degree"].as<size_t>();
    auto exec = (dim == 3) ? ((deg == 2) ? execute<3, 2> : execute<3, 1>)
                           : ((deg == 2) ? execute<2, 2> : execute<2, 1>);

    cout << setprecision(19) << endl;

    exec(args, inVertices, inElements);

    return 0;
}
