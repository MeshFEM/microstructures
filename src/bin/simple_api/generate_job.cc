////////////////////////////////////////////////////////////////////////////////
// generate_job.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Create a pattern optimization job intended to be used with the
//      isosurface inflator. This is particularly useful for generating sane
//      initial parameters (from the defaults) and from converting from offset
//      bounds to translation bounds:
//
//          The new isosurface inflator uses translations as variables because these
//          allow the printability constraints to be formulated more transparently.
//
//          However, often it is desirable to bound the offset of vertices from
//          known good positions (e.g. positions initialized from topology
//          enumeration). This utility will express such an offset bound as a
//          separate bound on each translation variable.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/04/2016 17:47:39
////////////////////////////////////////////////////////////////////////////////

#include <MeshFEM/Future.hh>
#include <isosurface_inflator/WireMesh.hh>
#include <pattern_optimization/PatternOptimizationJob.hh>
#include <nlohmann/json.hpp>
#include <CLI/CLI.hpp>

#include <iostream>
#include <cstdlib>
#include <string>
#include <cassert>
#include <stdexcept>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////

struct Args {
    std::string topology;
    std::vector<double> offsetBounds;
    std::vector<double> translationBounds;
    double defaultThickness = 0.07;
    std::vector<double> radiusBounds = { 0.04, 0.2 };
    std::vector<double> blendingBounds = { 0.005, 0.2};
    std::vector<double> elasticityTensor = { 1, 0 };
    std::vector<double> initialParams;
    std::string parameterConstraints;
    std::string symmetry = "orthotropic";
    bool limitedOffset = false;
};

// -----------------------------------------------------------------------------

int main(int argc, const char *argv[]) {
    Args args;

    // Parse arguments
    CLI::App app{"homogenize"};

    app.add_option("topology", args.topology, "Topology (line mesh)")->required()->check(CLI::ExistingFile);
    auto offset_opt = app.add_option("-o,--offsetBounds", args.offsetBounds, "offset bounds specifier (lower,upper)")->expected(2);
    auto translation_opt = app.add_option("-t,--translationBounds", args.translationBounds, "translation bounds specifier (lower,upper)")->expected(2);
    app.add_option("--defaultThickness", args.defaultThickness, "default thickness");
    app.add_option("-r,--radiusBounds", args.radiusBounds, "radius bounds specifier (lower,upper)")->expected(2);
    app.add_option("-b,--blendingBounds", args.blendingBounds, "blending bounds specifier (lower,upper)")->expected(2);
    app.add_option("-e,--elasticityTensor", args.elasticityTensor, "target tensor specifier (Young,Poisson)")->expected(2);
    app.add_option("-p,--initialParams", args.initialParams, "initial parameters (optional)");
    app.add_option("-c,--parameterConstraints", args.parameterConstraints, "parameter constraint expressions (semicolon-separated, optional)");
    app.add_set("-s,--symmetry", args.symmetry,
        { "cubic", "orthotropic", "diagonal", "square", "triply_periodic", "doubly_periodic", "non_periodic", },
        "symmetries to enforce (orthotropic (default), cubic, square, diagonal, triply_periodic, doubly_periodic)");
    app.add_flag("-L,--limitedOffset", args.limitedOffset, "Limit offset of nodes to within the base unit (0, 1)");

    translation_opt->excludes(offset_opt);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    // Infer dimension from pattern
    std::vector<MeshIO::IOVertex> vertices;
    std::vector<MeshIO::IOElement> elements;
    MeshIO::load(args.topology, vertices, elements);

    double zmag = 0;
    for (const auto &v : vertices) {
        zmag = std::max(zmag, std::abs(v[2]));
    }
    size_t dim = (zmag < 1e-2) ? 2 : 3;

    auto writeJob = [&args, dim] (auto wm) {
        nlohmann::json config = {
            {"dim", dim},
            {"offsetBounds", args.offsetBounds},
            {"translationBounds", args.translationBounds},
            {"defaultThickness", args.defaultThickness},
            {"radiusBounds", args.radiusBounds},
            {"blendingBounds", args.blendingBounds},
            {"elasticityTensor", args.elasticityTensor},
            {"initialParams", args.initialParams},
            {"parameterConstraints", args.parameterConstraints},
            {"limitedOffset", false},
        };

        auto job = PatternOptimization::jobFromWireMesh(wm, config);
        job->writeJobFile(std::cout);
    };

    const auto sym = args.symmetry;
    if      (sym == "cubic"          ) { writeJob(WireMesh<Symmetry::Cubic<>>         (vertices, elements)); }
    else if (sym == "orthotropic"    ) { writeJob(WireMesh<Symmetry::Orthotropic<>>   (vertices, elements)); }
    else if (sym == "diagonal"       ) { writeJob(WireMesh<Symmetry::Diagonal<>>      (vertices, elements)); }
    else if (sym == "square"         ) { writeJob(WireMesh<Symmetry::Square<>>        (vertices, elements)); }
    else if (sym == "triply_periodic") { writeJob(WireMesh<Symmetry::TriplyPeriodic<>>(vertices, elements)); }
    else if (sym == "doubly_periodic") { writeJob(WireMesh<Symmetry::DoublyPeriodic<>>(vertices, elements)); }
    else if (sym == "non_periodic"   ) { writeJob(WireMesh<Symmetry::NonPeriodic<>>(vertices, elements)); }
    else throw std::runtime_error("Unknown symmetry type: " + sym);
    return 0;
}
