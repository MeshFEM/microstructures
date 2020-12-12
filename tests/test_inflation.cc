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

void test_inflation(const std::string &mesher, const std::string &pattern) {
    size_t inflation_graph_radius = 2;

    IsosurfaceInflator inflator(mesher, true, pattern, inflation_graph_radius);
    std::vector<Real> params(inflator.defaultParameters());

    inflator.meshingOptions().load(default_meshing_options);

    // inflator.disablePostprocess();
    // inflator.enableCheapPostprocess();
    // inflator.setReflectiveInflator(!args.nonReflectiveInflator);
    // inflator.setGenerateFullPeriodCell(!args.ortho_cell);

    std::cout << "Dilating " << mesher << ":\n\t" << pattern << std::endl;
    inflator.inflate(params);

    #ifdef DUMP_OUTPUT
    if (MeshFEM::startswith(mesher, "2d_")) {
        MeshIO::save(mesher + ".obj", inflator.vertices(), inflator.elements());
    } else {
        std::cout << "Saving " << mesher << ".mesh" << std::endl;
        MeshIO::save(mesher + ".mesh", inflator.vertices(), inflator.elements());
    }
    #endif

    // Verify that normals have a zero z component in 2D
    if (MeshFEM::startswith(mesher, "2d_")) {
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
            auto normal = inflator.trackSignedDistanceGradient(inflator.vertices().at(maxZMagVtx).point);
            std::cout << "normal: " << normal.transpose() << std::endl;
        }
        REQUIRE(maxZMag <= 0.1);
    }
}

// -----------------------------------------------------------------------------

template<int N, typename SymmetryType>
void test_mesher(const std::string &topology, bool skewed = false) {
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

    PatternSignedDistance<double, WireMesh<SymmetryType>> sdf(wm);
    sdf.setParameters(wm.defaultParameters(), mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);

    // Change the pattern's meshing domain if we're forcing meshing of the
    // full TriplyPeriodic base cell.
    // if (generateFullPeriodCell && !reflectiveInflator)
    //     sdf.setBoundingBox(Symmetry::TriplyPeriodic<>::representativeMeshCell<Real>());

    // Mesh pattern and save result
    {
        std::vector<MeshIO::IOVertex> vertices_out;
        std::vector<MeshIO::IOElement> elements_out;

        mesher->meshInterfaceConsistently = true;
        mesher->mesh(sdf, vertices_out, elements_out);

        #ifdef DUMP_OUTPUT
        std::string basename = MeshFEM::replace_ext(MeshFEM::split(topology, "/").back(), "");
        std::string symmetry_name = SymmetryTraits<SymmetryType>::value;
        std::string postfix = (skewed ? "_skewed" : "");
        MeshIO::save(basename + "_" + symmetry_name + postfix + ".obj", vertices_out, elements_out);
        #endif
    }
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("inflate_default_pattern", "[isosurface_inflation]") {
    auto pattern_2d = DATA_DIR "patterns/2D/topologies/0001.obj";
    auto pattern_3d = DATA_DIR "patterns/3D/reference_wires/pattern0000.wire";

    SECTION("2d_orthotropic")     { test_inflation("2d_orthotropic", pattern_2d);     }
    SECTION("2d_diagonal")        { test_inflation("2d_diagonal", pattern_2d);        }
    SECTION("2d_non_periodic")    { test_inflation("2d_non_periodic", pattern_2d);    }
    SECTION("2d_doubly_periodic") { test_inflation("2d_doubly_periodic", pattern_2d); }
    SECTION("2d_square")          { test_inflation("2d_square", pattern_2d);          }
    #ifdef MICRO_WITH_CGAL
    SECTION("3d_cubic")           { test_inflation("cubic", pattern_3d);              }
    SECTION("3d_orthotropic")     { test_inflation("orthotropic", pattern_3d);        }
    // SECTION("3d_non_periodic")    { test_inflation("non_periodic", pattern_3d);       } // Not working
    // SECTION("3d_triply_periodic") { test_inflation("triply_periodic", pattern_3d);    } // Needs CGAL triply periodic meshing (4.13)
    SECTION("3d_square")          { test_inflation("square", pattern_3d);             }
    #endif
}

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("inflate_mesher_explicit", "[isosurface_inflation]") {
    SECTION("pattern_0001") {
        std::string pattern_2d = DATA_DIR "patterns/2D/topologies/0001.obj";
        std::string pattern_3d = DATA_DIR "patterns/3D/reference_wires/pattern0000.wire";

        SECTION("2d_orthotropic")     { test_mesher<2, Symmetry::Orthotropic<>>(pattern_2d);    }
        SECTION("2d_diagonal")        { test_mesher<2, Symmetry::Diagonal<>>(pattern_2d);       }
        SECTION("2d_doubly_periodic") { test_mesher<2, Symmetry::DoublyPeriodic<>>(pattern_2d); }
        SECTION("2d_square")          { test_mesher<2, Symmetry::Square<>>(pattern_2d);         }
    }

    SECTION("patter_0905") {
        std::string pattern_2d = DATA_DIR "patterns/2D/topologies/0905.obj";

        SECTION("2d_diagonal")        { test_mesher<2, Symmetry::Diagonal<>>(pattern_2d, true); }
    }
}
