////////////////////////////////////////////////////////////////////////////////

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/MSHFieldParser.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/Materials.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/util.h>
#include <CLI/CLI.hpp>
#include <json.hpp>
#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

using json = nlohmann::json;

////////////////////////////////////////////////////////////////////////////////

struct Args {
    std::string mesh;
    std::string material;
    std::string matFieldName;
    std::string boundaryConditions;
    std::string outputMSH;
    std::string dumpMatrix;
    int degree = 2;
    bool fullDegreeFieldOutput = false;
};

////////////////////////////////////////////////////////////////////////////////S

template<size_t _N, size_t _Deg>
void execute(const Args &args,
             const std::vector<MeshIO::IOVertex> &inVertices,
             const std::vector<MeshIO::IOElement> &inElements)
{
    size_t numElements = inElements.size();
    typedef LinearElasticity::Mesh<_N, _Deg> Mesh;
    using Simulator = LinearElasticity::Simulator<Mesh>;
    Simulator sim(inElements, inVertices);

    typedef ScalarField<Real> SField;

    if (fileExtension(args.material) == ".msh") {
        MSHFieldParser<_N> fieldParser(args.material);
        // Read heterogenous material from .msh file.
        // Guess isotropic or orhotropic based on fields present
        // Isotropic names: E nu
        // Orthotropic names: E_x E_y [E_z] nu_yx [nu_zx nu_zy] [mu_yz mu_zx] mu[_xy]
        auto domainSizeChecker = [=](const std::vector<SField> &fs) -> bool {
            return all_of(fs.begin(), fs.end(),
               [=](const SField &f) { return f.domainSize() == numElements; } );
        };
        std::runtime_error sizeErr("Material parameter fields of incorrect size.");
        std::runtime_error notFound("No complete material parameter field was found.");

        std::vector<SField> paramFields;
        std::vector<std::string> isotropicNames = { "E", "nu" };
        for (std::string name : isotropicNames) {
            name = args.matFieldName + name;
            try { paramFields.push_back(fieldParser.scalarField(name,
                        DomainType::PER_ELEMENT)); }
            catch (...) { /* Don't complain yet--try orthotropic */ }
        }
        if (paramFields.size() == 2) {
            if (!domainSizeChecker(paramFields)) throw sizeErr;
            // Valid isotropic material field--load it into simulator.
            LinearElasticity::ETensorStoreGetter<_N> store;
            for (size_t i = 0; i < sim.mesh().numElements(); ++i) {
                store().setIsotropic(paramFields[0][i], paramFields[1][i]);
                sim.mesh().element(i)->configure(store);
            }
            std::cout << "Loaded " << _N << "D isotropic material" << std::endl;
        }
        else {
            // If isotropic field wasn't found, try orthotropic.
            paramFields.clear();
            std::vector<std::vector<std::string> > orthotropicNames =
                { { "E_x", "E_y", "nu_yx", "mu" },
                  { "E_x", "E_y", "E_z", "nu_yx", "nu_zx", "nu_zy", "mu_yz", "mu_zx", "mu_xy" } };
            for (std::string name : orthotropicNames.at(_N - 2)) {
                name = args.matFieldName + name;
                try { paramFields.push_back(fieldParser.scalarField(name,
                            DomainType::PER_ELEMENT)); }
                catch (...) { throw notFound; }
            }
            if (!domainSizeChecker(paramFields)) throw sizeErr;
            // Valid orthotropic material field--load it into simulator.
            LinearElasticity::ETensorStoreGetter<_N> store;
            for (size_t i = 0; i < sim.mesh().numElements(); ++i) {
                if (_N == 2) {
                    store().setOrthotropic2D(
                        paramFields[0][i], paramFields[1][i],
                        paramFields[2][i], paramFields[3][i]);
                }
                else {
                    store().setOrthotropic3D(
                        paramFields[0][i], paramFields[1][i], paramFields[2][i],
                        paramFields[3][i], paramFields[4][i], paramFields[5][i],
                        paramFields[6][i], paramFields[7][i], paramFields[8][i]);
                }
                sim.mesh().element(i)->configure(store);
            }
            std::cout << "Loaded " << _N << "D Orthotropic material" << std::endl;
        }
    }
    else {
        // Read homogenous material from .material file (or use default material
        // if no file is given).
        Materials::Constant<_N> mat;
        if (args.material != "")
            mat.setFromFile(args.material);
        LinearElasticity::ETensorStoreGetter<_N> store(mat.getTensor());
        for (size_t i = 0; i < sim.mesh().numElements(); ++i)
            sim.mesh().element(i)->configure(store);
    }

    // Check if we're just dumping the stiffness matrix without simulating
    if ((args.dumpMatrix != "") && (args.boundaryConditions == "")) {
        typename Simulator::TMatrix K;
        sim.m_assembleStiffnessMatrix(K);
        K.sumRepeated();
        K.dumpBinary(args.dumpMatrix);
        return;
    }

    bool noRigidMotion;
    std::vector<PeriodicPairDirichletCondition<_N>> pps;
    ComponentMask pinTranslationComponents;
    auto bconds = readBoundaryConditions<_N>(args.boundaryConditions, sim.mesh().boundingBox(), noRigidMotion, pps, pinTranslationComponents);
    sim.applyTranslationPins(pinTranslationComponents);
    sim.applyBoundaryConditions(bconds);
    sim.applyPeriodicPairDirichletConditions(pps);
    if (noRigidMotion) sim.applyNoRigidMotionConstraint();

    if (args.dumpMatrix != "") sim.dumpSystem(args.dumpMatrix);


    BENCHMARK_START_TIMER_SECTION("Simulation");
    auto u = sim.solve();
    auto e = sim.averageStrainField(u);
    auto s = sim.averageStressField(u);
    auto f = sim.dofToNodeField(sim.neumannLoad());
    BENCHMARK_STOP_TIMER_SECTION("Simulation");

    bool linearSubsampleFields = !args.fullDegreeFieldOutput;

    MSHFieldWriter writer(args.outputMSH, sim.mesh(), linearSubsampleFields);
    writer.addField("u",      u, DomainType::PER_NODE);
    writer.addField("load",   f, DomainType::PER_NODE);
    if ((Simulator::Strain::Deg == 0) || linearSubsampleFields) {
        // Output constant (average) strain/stress for piecewise linear u
        writer.addField("strain", e, DomainType::PER_ELEMENT);
        writer.addField("stress", s, DomainType::PER_ELEMENT);
    } else {
        // Output full-degree per-element strain. (Wasteful since
        // strain fields are of degree - 1, but Gmsh/MSHFieldWriter
        // only supports full-degree ElementNodeData).
        auto linearField = sim.strainField(u);
        using Upsampled = SymmetricMatrixInterpolant<typename Simulator::SMatrix, _N, _Deg>;
        std::vector<Upsampled> upsampledField;
        upsampledField.reserve(linearField.size());
        for (const auto s: linearField) upsampledField.emplace_back(s);
        writer.addField("strain", upsampledField, DomainType::PER_ELEMENT);

        linearField = sim.stressField(u);
        upsampledField.clear();
        for (const auto s: linearField) upsampledField.emplace_back(s);
        writer.addField("stress", upsampledField, DomainType::PER_ELEMENT);
    }

    // // Write mat parameter fields
    // SField Ex(numElements), Ey(numElements), nuYX(numElements), mu(numElements);
    // for (size_t i = 0; i < sim.mesh().numElements(); ++i)
    //     sim.mesh().element(i)->E().getOrthotropic2D(Ex[i], Ey[i], nuYX[i], mu[i]);
    // writer.addField("E_x",    Ex,    DomainType::PER_ELEMENT);
    // writer.addField("E_y",    Ey,    DomainType::PER_ELEMENT);
    // writer.addField("nu_yx",  nuYX,  DomainType::PER_ELEMENT);
    // writer.addField("mu",     mu,    DomainType::PER_ELEMENT);

    sim.reportRegionSurfaceForces(u);
    writer.addField("Ku", sim.applyStiffnessMatrix(u), DomainType::PER_NODE);

    BENCHMARK_REPORT();
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    Args args;

    // Parse arguments
    CLI::App app{"simulate"};

    app.add_option("mesh", args.mesh, "input mesh")->required()->check(CLI::ExistingFile);
    app.add_option("-m,--material", args.material , "simulation material material");
    app.add_option("-f,--matFieldName", args.matFieldName , "name of material field to load from .msh passed as --material");
    app.add_option("-b,--boundaryConditions", args.boundaryConditions , "boundary conditions");
    app.add_option("-o,--outputMSH", args.outputMSH , "output mesh");
    app.add_option("--dumpMatrix", args.dumpMatrix , "dump system matrix in triplet format");
    app.add_option("-d,--degree", args.degree , "FEM degree (1 or 2)");
    app.add_flag("-D,--fullDegreeFieldOutput", args.fullDegreeFieldOutput, "Output full-degree nodal fields (don't do piecewise linear subsample)");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    if (app.count("dumpMatrix") == 0) {
        if (app.count("outputMSH") == 0) {
            std::string e = "Error: must specify output msh file (unless dumping a stiffness matrix)";
            return app.exit(CLI::ValidationError(e));
        }
    }

    if (app.count("outputMSH") && (app.count("boundaryConditions") == 0)) {
        return app.exit(CLI::ValidationError("Error: must specify boundary conditions to run a simulation"));
    }

    std::vector<MeshIO::IOVertex>  inVertices;
    std::vector<MeshIO::IOElement> inElements;

    auto type = load(args.mesh, inVertices, inElements, MeshIO::FMT_GUESS, MeshIO::MESH_GUESS);

    // Infer dimension from mesh type.
    size_t dim;
    if (type == MeshIO::MESH_TET) dim = 3;
    else if (type == MeshIO::MESH_TRI) dim = 2;
    else throw std::runtime_error("Mesh must be pure triangle or tet.");

    // Look up and run appropriate simulation instantiation.
    auto exec = (dim == 3) ? ((args.degree == 2) ? execute<3, 2> : execute<3, 1>)
                           : ((args.degree == 2) ? execute<2, 2> : execute<2, 1>);
    exec(args, inVertices, inElements);

    return 0;
}
