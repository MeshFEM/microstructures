////////////////////////////////////////////////////////////////////////////////
// Example run:
// ./isosurface_inflator/stitch_cells_cli $MICRO_DIR/isosurface_inflator/tests/patch.json -o out.obj
////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/StitchedWireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>
////////////////////////////////////////////////////////////////////////////////

using json = nlohmann::json;
using WireMeshBasePtr = std::shared_ptr<WireMeshBase>;

////////////////////////////////////////////////////////////////////////////////

std::string lowercase(std::string data) {
    std::transform(data.begin(), data.end(), data.begin(), ::tolower);
    return data;
}

#define TRY_SYMMETRY(s, x, p)                                  \
    if (lowercase(x) == lowercase(#s))                         \
    {                                                          \
        return std::make_shared<WireMesh<Symmetry::s<>>>((p)); \
    }

#define TRY_KEY_VAL(s, a, x, p)                                \
    if (lowercase(x) == lowercase(#a))                         \
    {                                                          \
        return std::make_shared<WireMesh<Symmetry::s<>>>((p)); \
    }

WireMeshBasePtr load_wire_mesh(const std::string &sym, const std::string &path) {
    TRY_SYMMETRY(Square, sym, path);
    TRY_SYMMETRY(Cubic, sym, path);
    TRY_SYMMETRY(Orthotropic, sym, path);
    TRY_SYMMETRY(Diagonal, sym, path);
    TRY_KEY_VAL(DoublyPeriodic, Doubly_Periodic, sym, path);
    TRY_KEY_VAL(TriplyPeriodic, Triply_Periodic, sym, path);
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

template<size_t N>
void execute(const std::string &patchFilename, const std::string &meshingOptions, const std::string &outname) {
    // Load patch config
    json patch;
    std::ifstream patchFile(patchFilename);
    try {
        patchFile >> patch;
    } catch (...) {
        std::cerr << "Error parsing the json file" << std::endl;
        return;
    }

    // Create mesher and load meshing options
    std::unique_ptr<MesherBase> mesher;
    if (N == 2) { mesher = std::make_unique<MidplaneMesher>(); }
    if (N == 3) { mesher = std::make_unique<IGLSurfaceMesherMC>(); }
    if (!meshingOptions.empty()) { mesher->meshingOptions.load(meshingOptions); }

    // Assign topologies and parameters to each cell
    NDCubeArray<WireMeshBasePtr, N, 3> topologyGrid;
    NDCubeArray<std::vector<double>, N, 3> parameterGrid;

    for (auto entry : patch) {
        NDArrayIndex<N> index;
        std::copy_n(entry["index"].begin(), N, index.idxs.begin());
        topologyGrid(index) = load_wire_mesh(entry["symmetry"], entry["pattern"]);
        parameterGrid(index) = entry["params"].get<std::vector<double>>();
    }

    auto swm = make_stitched_wire_mesh(topologyGrid);

    auto params = swm.paramsFromParamGrid(parameterGrid);

    PatternSignedDistance<double, StitchedWireMesh<N>> sdf(swm);
    sdf.setParameters(params, mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);

    std::vector<MeshIO::IOVertex> vertices;
    std::vector<MeshIO::IOElement> elements;

    mesher->meshInterfaceConsistently = true;
    mesher->mesh(sdf, vertices, elements);

    MeshIO::save(outname, vertices, elements);
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
    // Default arguments
    struct {
        std::string patch_config;
        std::string meshing_options = "";
        std::string output = "out.msh";
        int dimension = 2;
    } args;

    // Parse arguments
    CLI::App app{"stitch_cells_cli"};
    app.add_option("patch,-p,--patch", args.patch_config, "3x3 patch description (json file).")->required();
    app.add_option("output,-o,--output", args.output, "Output triangle mesh.");
    app.add_option("-m,--meshing", args.meshing_options, "Meshing options (json file).");
    app.add_option("-d,--dimension", args.dimension, "Dimension of the problem (2 or 3).")->check(CLI::Range(2, 3));
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    if (args.dimension == 2) { execute<2>(args.patch_config, args.meshing_options, args.output); }
    else if (args.dimension == 3) { execute<3>(args.patch_config, args.meshing_options, args.output); }

    return 0;
}
