////////////////////////////////////////////////////////////////////////////////
#include "test_common.hh"
#include <isosurface_inflator/IsosurfaceInflator.hh>
#include <isosurface_inflator/MeshingOptions.hh>
#include <isosurface_inflator/IsosurfaceInflatorConfig.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <MeshFEM/StringUtils.hh>
#include <catch2/catch.hpp>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

namespace {

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

// -----------------------------------------------------------------------------

// From MeshIO to Eigen (triangle mesh)
void meshfem_to_eigen(
    const std::vector<MeshIO::IOVertex> &VI,
    const std::vector<MeshIO::IOElement> &FI,
    Eigen::MatrixXd &VO,
    Eigen::MatrixXi &FO)
{
    VO.resize(VI.size(), 3);
    FO.resize(FI.size(), 3);
    for (int v = 0; v < (int) VI.size(); ++v) {
        for (int c = 0; c < 3; ++c) {
            VO(v, c) = VI[v][c];
        }
    }
    for (int f = 0; f < (int) FI.size(); ++f) {
        assert(FI[f].size() == 3);
        for (int lv = 0; lv < 3; ++lv) {
            FO(f, lv) = FI[f][lv];
        }
    }
}

// -----------------------------------------------------------------------------

template<int N, typename SymmetryType>
void test_aabb_meshing(const std::string &topology, bool skewed = false) {
    // Create mesher and load meshing options
    std::unique_ptr<MesherBase> mesher;
    if (N == 2) { mesher = std::make_unique<MidplaneMesher>(); }
    if (N == 3) { mesher = std::make_unique<IGLSurfaceMesherMC>(); }

    nlohmann::json opt = default_meshing_options;
    if (skewed) {
        opt["jacobian"] = {
            0.537285, -0.537285, 0,
            0.930605,  0.930605, 0,
            0       ,  0       , 1
        };
    }
    mesher->meshingOptions.load(opt);

    // Load input wire mesh
    std::vector<MeshIO::IOVertex> vertices_in;
    std::vector<MeshIO::IOElement> elements_in;
    MeshIO::load(topology, vertices_in, elements_in);

    // Setup SDF function
    WireMesh<SymmetryType> wm(vertices_in, elements_in);

    // Change the pattern's meshing domain if we're forcing meshing of the
    // full TriplyPeriodic base cell.
    // if (generateFullPeriodCell && !reflectiveInflator)
    //     sdf.setBoundingBox(Symmetry::TriplyPeriodic<>::representativeMeshCell<Real>());

    // Mesh pattern with and without AABB
    Eigen::MatrixXd V1, V2;
    Eigen::MatrixXi F1, F2;
    {
        PatternSignedDistance<double, WireMesh<SymmetryType>> sdf(wm);
        sdf.setParameters(wm.defaultParameters(), mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);
        std::vector<MeshIO::IOVertex> vertices_out;
        std::vector<MeshIO::IOElement> elements_out;
        mesher->meshInterfaceConsistently = true;
        mesher->mesh(sdf, vertices_out, elements_out);
        meshfem_to_eigen(vertices_out, elements_out, V1, F1);
    }

    {
        PatternSignedDistance<double, WireMesh<SymmetryType>> sdf(wm);
        sdf.setUseAabbTree(true);
        sdf.setParameters(wm.defaultParameters(), mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);
        std::vector<MeshIO::IOVertex> vertices_out;
        std::vector<MeshIO::IOElement> elements_out;
        mesher->meshInterfaceConsistently = true;
        mesher->mesh(sdf, vertices_out, elements_out);
        meshfem_to_eigen(vertices_out, elements_out, V2, F2);
    }

    REQUIRE(V1.rows() == V2.rows());
    REQUIRE(F1.rows() == F2.rows());
    std::cout << (V1 - V2).norm() << std::endl;
    REQUIRE((V1 - V2).norm() <= 1e-10);
    REQUIRE((F1 - F2).norm() == 0);
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("inflate_aabb", "[isosurface_inflation]") {
    SECTION("pattern_0001") {
        std::string pattern_2d = DATA_DIR "patterns/2D/topologies/0001.obj";
        std::string pattern_3d = DATA_DIR "patterns/3D/reference_wires/pattern0000.wire";

        SECTION("2d_orthotropic")     { test_aabb_meshing<2, Symmetry::Orthotropic<>>(pattern_2d);    }
        SECTION("2d_diagonal")        { test_aabb_meshing<2, Symmetry::Diagonal<>>(pattern_2d);       }
        SECTION("2d_doubly_periodic") { test_aabb_meshing<2, Symmetry::DoublyPeriodic<>>(pattern_2d); }
        SECTION("2d_square")          { test_aabb_meshing<2, Symmetry::Square<>>(pattern_2d);         }
    }

    SECTION("patter_0905") {
        std::string pattern_2d = DATA_DIR "patterns/2D/topologies/0905.obj";

        SECTION("2d_diagonal")        { test_aabb_meshing<2, Symmetry::Diagonal<>>(pattern_2d, true); }
    }
}
