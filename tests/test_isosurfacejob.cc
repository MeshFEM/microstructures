////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/WireMesh.hh>
#include <pattern_optimization/PatternOptimizationJob.hh>
#include <nlohmann/json.hpp>
#include <catch2/catch.hpp>
////////////////////////////////////////////////////////////////////////////////

template<typename SymmetryType>
void test_pattern(const std::string &filename) {
    // Infer dimension from pattern
    std::vector<MeshIO::IOVertex> vertices;
    std::vector<MeshIO::IOElement> elements;
    MeshIO::load(filename, vertices, elements);

    double zmag = 0;
    for (const auto &v : vertices) {
        zmag = std::max(zmag, std::abs(v[2]));
    }
    size_t dim = (zmag < 1e-2) ? 2 : 3;

    if ( (dim == 3 && std::is_same<SymmetryType, Symmetry::Diagonal<>>::value)
        || (dim == 3 && std::is_same<SymmetryType, Symmetry::DoublyPeriodic<>>::value))
    {
        // * Diagonal symmetry not implemented in 3D
        // * Doubly periodic is for 2D only
        return;
    }

    // Configuration
    nlohmann::json config = R"({
        "offsetBounds": [],
        "translationBounds": [],
        "defaultThickness": 0.07,
        "radiusBounds": [0.04, 0.2],
        "blendingBounds": [0.005, 0.2],
        "elasticityTensor": [1, 0],
        "initialParams": [],
        "parameterConstraints": "",
        "limitedOffset": false
    })"_json;
    config["dim"] = dim;

    WireMesh<SymmetryType> wm(vertices, elements);
    auto job = PatternOptimization::jobFromWireMesh(wm, config);
    std::cout << job->getJson().dump(4) << std::endl;
}

template<typename SymmetryType>
void test_all_patterns(const std::vector<std::string> &patterns) {
    for (auto filename : patterns) {
        test_pattern<SymmetryType>(filename);
    }
}

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("generate_isosurface_job", "[isosurface_job]") {

    const std::vector<std::string> patterns = {
        DATA_DIR "patterns/2D/topologies/0001.obj",
        DATA_DIR "patterns/3D/reference_wires/pattern0000.wire",
    };

    SECTION("cubic")           { test_all_patterns<Symmetry::Cubic<>>(patterns);          }
    SECTION("orthotropic")     { test_all_patterns<Symmetry::Orthotropic<>>(patterns);    }
    SECTION("diagonal")        { test_all_patterns<Symmetry::Diagonal<>>(patterns);       }
    SECTION("square")          { test_all_patterns<Symmetry::Square<>>(patterns);         }
    SECTION("triply_periodic") { test_all_patterns<Symmetry::TriplyPeriodic<>>(patterns); }
    SECTION("doubly_periodic") { test_all_patterns<Symmetry::DoublyPeriodic<>>(patterns); }
    SECTION("non_periodic")    { test_all_patterns<Symmetry::NonPeriodic<>>(patterns);    }
}
