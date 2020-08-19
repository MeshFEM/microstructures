////////////////////////////////////////////////////////////////////////////////
// Inflate and stitch a pattern on a given quad-mesh.
// Given a quad-mesh .obj as input, and a json describing the per-quad patterns,
// builds an inflation graph and SDF accordingly, and inflate the desired cell.
////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <isosurface_inflator/WireQuadMesh.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>
////////////////////////////////////////////////////////////////////////////////

using json = nlohmann::json;
using WireMeshBasePtr = std::shared_ptr<WireMeshBase>;

////////////////////////////////////////////////////////////////////////////////

template<size_t N>
void execute(const std::string &inputMesh,
    const std::string &meshingOptions,
    const std::string &patternParameters,
    const std::string &outputMesh,
    int index)
{
    // Load pattern parameters
    json pattern;
    {
        std::ifstream in(patternParameters);
        try {
            in >> pattern;
        } catch (...) {
            std::cerr << "Error parsing the json file." << std::endl;
            return;
        }
    }

    // Create mesher and load meshing options
    std::unique_ptr<MesherBase> mesher;
    if (N == 2) { mesher = std::make_unique<MidplaneMesher>(); }
    if (N == 3) { mesher = std::make_unique<IGLSurfaceMesherMC>(); }
    if (!meshingOptions.empty()) { mesher->meshingOptions.load(meshingOptions); }

    // Load input quad mesh
    std::vector<MeshIO::IOVertex> VI;
    std::vector<MeshIO::IOElement> FI;
    auto mesh_type = MeshIO::load(inputMesh, VI, FI);
    if (mesh_type != MeshIO::MESH_QUAD) {
        std::cerr << "Invalid element type, expected quads." << std::endl;
        return;
    }

    // Setup SDF function
    WireQuadMesh wm(VI, FI, pattern);
    wm.setActiveQuad(index);

    PatternSignedDistance<double, WireQuadMesh, WireQuadMesh::MapToBaseUnit> sdf(wm);
    sdf.setParameters(wm.params(), mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);
    sdf.setMapFunctor(wm.mapFunctor());
    sdf.setBoundingBox(wm.boundingBox());

    // Mesh active quad and save result
    {
        std::vector<MeshIO::IOVertex> vertices;
        std::vector<MeshIO::IOElement> elements;

        mesher->meshInterfaceConsistently = true;
        mesher->mesh(sdf, vertices, elements);

        MeshIO::save(outputMesh, vertices, elements);
    }
}

////////////////////////////////////////////////////////////////////////////////

// Input json for the quad mesh:
// [
//     {
//         "params": [...],
//         "jacobian": [...],
//         "symmetry": "diagonal",
//         "pattern": "0905.obj"
//     },
//     {
//         "params": [...],
//         "jacobian": [...],
//         "symmetry": "diagonal",
//         "pattern": "0905.obj"
//     },
//     ...
// ]


int main(int argc, char * argv[]) {
    // Default arguments
    struct {
        std::string input;
        std::string output = "out.msh";
        std::string meshing_options = "";
        std::string params = "";
        int active_quad = 0;
    } args;

    // Parse arguments
    CLI::App app{"stitch_cells_cli"};
    app.add_option("input,-i,--input", args.input, "Input quad mesh.")->required();
    app.add_option("output,-o,--output", args.output, "Output triangle mesh.");
    app.add_option("-x,--active_quad", args.active_quad, "Index of the quad to mesh.");
    app.add_option("-m,--meshing", args.meshing_options, "Meshing options (json file).")->check(CLI::ExistingFile);
    app.add_option("-p,--params", args.params, "Pattern parameters per quad (json file).")->check(CLI::ExistingFile);
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    execute<2>(args.input, args.meshing_options, args.params, args.output, args.active_quad);

    return 0;
}
