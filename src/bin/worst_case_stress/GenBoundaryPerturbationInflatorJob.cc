////////////////////////////////////////////////////////////////////////////////
// GenBoundaryPerturbationInflatorJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Create a job for the boundary perturbation inflator.
*/
////////////////////////////////////////////////////////////////////////////////

#include <MeshFEM/Future.hh>
#include <inflators/wrappers/BoundaryPerturbationInflator.hh>
#include <pattern_optimization/PatternOptimizationJob.hh>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <iostream>
#include <cstdlib>
#include <string>
#include <cassert>
#include <stdexcept>
#include <cmath>

namespace po = boost::program_options;
using namespace std;

[[ noreturn ]] void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: GenBoundaryPerturbationInflatorJob [options] shape.msh" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("shape", po::value<string>(), "Shape (mesh)")
        ;
    po::positional_options_description p;
    p.add("shape", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",                                        "Produce this help message")
        ("maxOffset,o",      po::value<double>(),                             "maximum offset")
        ("elasticityTensor,e",  po::value<string>()->default_value("1,0"),    "target tensor specifier (Young,Poisson)")
        ("targetVolume,v",  po::value<double>()->default_value(0.0),          "target volume")
        ("targetInitialVolume",                                               "set target volume as initial one")
        ;

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;

    if (vm.count("shape") == 0) {
        cout << "Error: must specify shape mesh" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

vector<Real> parseVecArg(const string &arg, size_t N) {
    vector<string> components;
    boost::split(components, arg, boost::is_any_of(","));
    if (components.size() != N) throw runtime_error("Invalid number of components in argument " + arg);

    vector<Real> result;
    for (const string &c : components)
        result.push_back(stod(c));
    return result;
}

int main(int argc, const char *argv[])
{
    auto args = parseCmdLine(argc, argv);

    // Infer dimension from pattern
    vector<MeshIO::IOVertex> vertices;
    vector<MeshIO::IOElement> elements;
    MeshIO::load(args.at("shape").as<string>(), vertices, elements);

    Real zmag = 0;
    for (const auto &v : vertices)
        zmag = max(zmag, std::abs(v[2]));
    size_t dim = (zmag < 1e-2) ? 2 : 3;

    unique_ptr<InflatorBase> infl;

    unique_ptr<BoundaryPerturbationInflator<2>> inflator = Future::make_unique<BoundaryPerturbationInflator<2>>(vertices, elements, 1e-10);

    vector<Real> targetModuli = parseVecArg(args["elasticityTensor"].as<string>(), 2);

    unique_ptr<PatternOptimization::JobBase> job;
    auto job2D = Future::make_unique<PatternOptimization::Job<2>>();
    job2D->targetMaterial.setIsotropic(targetModuli[0], targetModuli[1]);
    job = move(job2D);

    if (args.count("maxOffset")) {
        std::vector<std::vector<Real>> boundsPerParam = inflator->customBounds(args["maxOffset"].as<double>(), 1e-3);

        for (size_t p = 0; p < boundsPerParam.size(); ++p) {
            job->varLowerBounds.emplace(p, boundsPerParam[p][0]);
            job->varUpperBounds.emplace(p, boundsPerParam[p][1]);
        }
    }

    if (args.count("targetInitialVolume")) {
        job->targetVolume = inflator->volume();
    }
    else {
        job->targetVolume = args["targetVolume"].as<double>();
    }

    // Fake value for other bounds
    job->translationBounds = {0.1, 0.8};
    job->radiusBounds   = {0.1, 0.2};
    job->blendingBounds = {0.0, 0.1};

    job->initialParams = std::vector<Real>(inflator->numParameters(), 0.0);

    job->writeJobFile(cout);

    return 0;
}
