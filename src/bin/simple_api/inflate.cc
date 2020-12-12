////////////////////////////////////////////////////////////////////////////////

#include <isosurface_inflator/IsosurfaceInflator.hh>
#include <isosurface_inflator/MeshingOptions.hh>
#include <isosurface_inflator/IsosurfaceInflatorConfig.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/Parallelism.hh>
#include <CLI/CLI.hpp>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

struct Args {
    std::string mesher;
    std::string wire;
    std::string outMSH;
    std::string dumpInflationGraph;
    std::string dumpReplicatedGraph;
    std::string mopts;
    std::string params;
    std::string paramsFile;
    size_t inflation_graph_radius = 2;
    std::string dumpShapeVelocities;
    std::string loadMesh;
    std::string rasterize;
    std::string rasterResolution = "20x20x20";
    std::string dumpBaseUnitGraph;
    bool nonReflectiveInflator = false;
    bool disablePostprocessing = false;
    bool cheapPostprocessing = false;
    bool ortho_cell = false;
    bool assertPlanarNormals = false;
    #if MICRO_WITH_TBB
    int numProcs = tbb::task_scheduler_init::default_num_threads();
    #endif
};

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    Args args;

    // Parse arguments
    CLI::App app{"inflate"};

    app.add_option("mesher", args.mesher, "name of mesher to use")->required();
    app.add_option("wire",   args.wire, "input wire file")->required()->check(CLI::ExistingFile);
    app.add_option("outMSH", args.outMSH, "output msh file");
    app.add_option("-D,--dumpInflationGraph",  args.dumpInflationGraph,     "Output the inflation graph to the specified path");
    app.add_option("-R,--dumpReplicatedGraph", args.dumpReplicatedGraph,    "Output the replicated pattern graph to the specified path");
    app.add_option("-m,--mopts",               args.mopts,                  "Meshing options file")->check(CLI::ExistingFile);
    app.add_option("-p,--params",              args.params,                 "Pattern parameters");
    app.add_option("--paramsFile",             args.paramsFile,             "Pattern parameters file")->check(CLI::ExistingFile);
    app.add_option("--inflation_graph_radius", args.inflation_graph_radius, "Number of edges to traverse outward from the symmetry cell when building the inflation graph");
    app.add_option("-S,--dumpShapeVelocities", args.dumpShapeVelocities,    "Dump the shape velocities for debugging");
    app.add_option("-M,--loadMesh",            args.loadMesh,               "Skip meshing process, loading existing mesh instead (for debugging)");
    app.add_option("-r,--rasterize",           args.rasterize,              "Rasterize the pattern to an indicator field on a regular grid");
    app.add_option("--rasterResolution",       args.rasterResolution,       "Size of the rasterization grid (2D or 3D)");
    app.add_option("-B,--dumpBaseUnitGraph",   args.dumpBaseUnitGraph,      "Output the base unit inflation graph to the specified path");
    app.add_flag("--nonReflectiveInflator",    args.nonReflectiveInflator,  "use non-reflective inflator (reflective by default)");
    app.add_flag("-d,--disablePostprocessing", args.disablePostprocessing,  "Disable post-processing of mesher output");
    app.add_flag("--cheapPostprocessing",      args.cheapPostprocessing,    "Cheap post-processing of mesher output");
    app.add_flag("-O,--ortho_cell",            args.ortho_cell,             "Generate the ortho cell only (for ortho-cell meshers)");
    app.add_flag("--assertPlanarNormals",      args.assertPlanarNormals,    "Verify that normals have a zero z component (relevant in 2D)");
    #if MICRO_WITH_TBB
    app.add_option("--numProcs", args.numProcs, "Number of threads to use for TBB parallelism (CGAL mesher, etc.)");
    #endif

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

#if MICRO_WITH_TBB
    if (args.numProcs > tbb::task_scheduler_init::default_num_threads()) {
        std::cerr << "WARNING: specifying more than the default number of TBB threads." << std::endl;
    }
    tbb::task_scheduler_init init(args.numProcs);
#endif

    IsosurfaceInflator inflator(args.mesher, true, args.wire, args.inflation_graph_radius);

    std::vector<Real> params(inflator.defaultParameters());
    if (args.params.size()) {
        // Split up params.
        std::istringstream ss(args.params);
        params.clear();
        Real x;
        while (ss >> x) {
            params.push_back(x);
        }
        std::cout << params.size() << std::endl;
    } else if (!args.paramsFile.empty()) {
        std::ifstream paramsFile(args.paramsFile, std::ifstream::in);
        if (paramsFile) {
            std::string line;
            while (getline(paramsFile, line)) {
                std::istringstream ss(line);
                params.clear();
                Real x;
                while (ss >> x) {
                    params.push_back(x);
                }
            }
        }
    } else {
        std::cout << "Inflating default parameters: " << std::endl;
        for (Real p : params) {
            std::cout << p << "\t";
        }
        std::cout << std::endl;
    }

    auto &config = IsosurfaceInflatorConfig::get();
    if (args.dumpInflationGraph.size()) {
        // dump inflation graph directly without running inflation
        inflator.dumpInflationGraph(args.dumpInflationGraph, params);
    }
    if (args.dumpReplicatedGraph.size()) {
        config.replicatedGraphPath = args.dumpReplicatedGraph;
    }
    if (args.dumpBaseUnitGraph.size()) {
        config.baseUnitGraphPath = args.dumpBaseUnitGraph;
    }

    if (args.mopts.size()) {
        inflator.meshingOptions().load(args.mopts);
    }
    if (args.disablePostprocessing) {
        inflator.disablePostprocess();
    }
    if (args.cheapPostprocessing) {
        inflator.enableCheapPostprocess();
    }
    inflator.setReflectiveInflator(!args.nonReflectiveInflator);

    if (args.dumpShapeVelocities.size()) {
        inflator.meshingOptions().debugSVelPath = args.dumpShapeVelocities;
    }
    if (args.loadMesh.size()) {
        inflator.meshingOptions().debugLoadMeshPath = args.loadMesh;
    }

    if (args.outMSH.size()) {
        inflator.setGenerateFullPeriodCell(!args.ortho_cell);
        inflator.inflate(params);

        MeshIO::save(args.outMSH, inflator.vertices(), inflator.elements());

        if (args.assertPlanarNormals) {
            const auto &n = inflator.vertexNormals();
            double maxZMag = 0;
            size_t maxZMagVtx = 0;
            for (size_t vi = 0; vi < n.size(); ++vi) {
                Real zmag = std::abs(n[vi][2]);
                if (zmag > maxZMag) {
                    maxZMag = zmag;
                    maxZMagVtx = vi;
                }
            }

            if (maxZMag > 0.1) {
                std::cerr << "Large normal z component: " << maxZMag << std::endl;
                auto n = inflator.trackSignedDistanceGradient(inflator.vertices().at(maxZMagVtx).point);
                std::cout << "normal: " << n.transpose() << std::endl;
            }
        }
    }

    if (args.rasterize.size()) {
        inflator.rasterize(params, args.rasterResolution, args.rasterize);
    }

    BENCHMARK_REPORT();
}
