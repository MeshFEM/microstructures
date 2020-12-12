#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#include "test_common.hh"
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IsosurfaceInflator.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <isosurface_inflator/WireQuadMesh.hh>
#include <MeshFEM/StringUtils.hh>
#include <catch2/catch.hpp>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

const nlohmann::json default_meshing_options = R"({
    "domainErrorBound"        : 1e-5,
    "facetAngle"              : 30.0,
    "facetSize"               : 0.02,
    "facetDistance"           : 1e-3,
    "cellSize"                : 0.50,
    "edgeSize"                : 0.015,
    "cellRadiusEdgeRatio"     : 2.0,

    "marchingSquaresGridSize" : 128,
    "marchingSquaresCoarsening": 1,
    "marchingCubesGridSize"   : 128,

    "forceConsistentInterfaceMesh": false,
    "forceMSGridSize": true,

    "maxArea"                 : 0.001,
    "featureAngleThreshold"   : 0.7853981633974483
})"_json;

////////////////////////////////////////////////////////////////////////////////

template<int N, typename SymmetryType>
void test_quad_mesh(const std::string &mesh, const std::string &topology) {
    // Create mesher and load meshing options
    std::unique_ptr<MesherBase> mesher;
    if (N == 2) { mesher = std::make_unique<MidplaneMesher>(); }
    if (N == 3) { mesher = std::make_unique<IGLSurfaceMesherMC>(); }

    // Load input quad mesh
    std::vector<MeshIO::IOVertex> vertices_in;
    std::vector<MeshIO::IOElement> elements_in;
    auto mesh_type = MeshIO::load(mesh, vertices_in, elements_in);
    if (mesh_type != MeshIO::MESH_QUAD) {
        std::cerr << "Invalid element type, expected quads." << std::endl;
        return;
    }

    // Build dummy json input
    json data = json::array();
    {
        std::vector<MeshIO::IOVertex> VI;
        std::vector<MeshIO::IOElement> FI;
        MeshIO::load(topology, VI, FI);
        WireMesh<SymmetryType> wm(VI, FI);
        for (int i = 0; i < (int) elements_in.size(); ++i) {
            json entry;
            entry["params"] = wm.defaultParameters();
            entry["symmetry"] = SymmetryTraits<SymmetryType>::value;
            entry["pattern"] = topology;
            data.push_back(entry);
        }
    }

    // Wire mesh embedded into a quad mesh
    WireQuadMesh wm(vertices_in, elements_in, data);
    Eigen::VectorXd A = wm.areas();

    for (int index : {-1, 0}) {
        // Set meshing options
        json opt = default_meshing_options;
        double factor = (index < 0 ? A.sum() / A.mean() : 1.0);
        opt["maxArea"] = opt["maxArea"].get<double>() / factor;
        opt["marchingSquaresGridSize"] = std::round(opt["marchingSquaresGridSize"].get<double>() * std::sqrt(factor));
        mesher->meshingOptions.load(opt);

        // Set SDF
        wm.setActiveQuad(index);
        PatternSignedDistance<double, WireQuadMesh, WireQuadMesh::MapToBaseUnit> sdf(wm);
        sdf.setParameters(wm.params(), mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);
        sdf.setMapFunctor(wm.mapFunctor());
        sdf.setBoundingBox(wm.boundingBox());

        // Mesh pattern and save result
        {
            std::vector<MeshIO::IOVertex> vertices_out;
            std::vector<MeshIO::IOElement> elements_out;

            mesher->meshInterfaceConsistently = true;
            mesher->mesh(sdf, vertices_out, elements_out);

            #ifdef DUMP_OUTPUT
            std::cout << "Saving stuff" << std::endl;
            std::string pattern_name = MeshFEM::replace_ext(MeshFEM::split(topology, "/").back(), "");
            std::string background_name = MeshFEM::replace_ext(MeshFEM::split(mesh, "/").back(), "");
            std::string suffix = (index < 0 ? "whole" : "q" + std::to_string(index));
            MeshIO::save(pattern_name + "_" + background_name + "_" + SymmetryTraits<SymmetryType>::value + "_" + suffix + ".obj", vertices_out, elements_out);
            #endif
        }
    }
}

// -----------------------------------------------------------------------------

template<typename SymmetryType>
void test_quad_meshes(const std::vector<std::string> &meshes, const std::string &pattern) {
    for (auto filename : meshes) {
        test_quad_mesh<2, SymmetryType>(filename, pattern);
    }
}

////////////////////////////////////////////////////////////////////////////////

#ifndef SECTION
#define SECTION(x)
#endif

// int main(void) {
TEST_CASE("inflate_and_stitch", "[isosurface_inflation]") {

    SECTION("pattern_0001") {
        std::string pattern_2d = DATA_DIR "patterns/2D/topologies/0001.obj";

        std::vector<std::string> meshes = {
            DATA_DIR "tests/quad_grid_orient_perfect.obj",
            DATA_DIR "tests/quad_grid_orient_fuzzy.obj",
            DATA_DIR "tests/quad_irregular_orient_fuzzy.obj",
            DATA_DIR "tests/quad_self_repeating.obj",
            DATA_DIR "tests/quad_rhombi_regular.obj",
        };

        SECTION("2d_square")          { test_quad_meshes<Symmetry::Square<>>(meshes, pattern_2d);          }
        SECTION("2d_diagonal")        { test_quad_meshes<Symmetry::Diagonal<>>(meshes, pattern_2d);        }
        SECTION("2d_orthotropic")     { test_quad_meshes<Symmetry::Orthotropic<>>(meshes, pattern_2d);     }
        SECTION("2d_doubly_periodic") { test_quad_meshes<Symmetry::DoublyPeriodic<>>(meshes, pattern_2d);  }
    }

    SECTION("pattern_0905") {
        std::string pattern_2d = DATA_DIR "patterns/2D/topologies/0905.obj";

        std::vector<std::string> meshes = {
            DATA_DIR "tests/quad_grid_orient_perfect.obj",
            DATA_DIR "tests/quad_grid_orient_fuzzy.obj",
            DATA_DIR "tests/quad_irregular_orient_fuzzy.obj",
            DATA_DIR "tests/quad_self_repeating.obj",
            // DATA_DIR "tests/quad_rhombi_regular.obj", // TODO: DEBUG SEGFAULT HERE
        };

        SECTION("2d_diagonal")        { test_quad_meshes<Symmetry::Diagonal<>>(meshes, pattern_2d); }
    }
}

#endif
