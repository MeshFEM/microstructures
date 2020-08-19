////////////////////////////////////////////////////////////////////////////////
// PatternOptimizationJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Reads in JSON files specifying an optimization job. These job
//      specifications include a target material and the initial pattern
//      parameters from which to start optimization.
//      {
//          "dim": 2,
//          "target": {
//              <material_spec>
//          },
//          "initial_params": [ ... ],
//          "radiusBounds": [ low, high ],
//          "translationBounds": [ low, high ],
//          "paramConstraints": [ "p1 = p2 + 5", ... ],
//          "bounds": [{ "var": 0, "lower": 1, "upper": 2 }]
//      }
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/04/2014 17:15:02
////////////////////////////////////////////////////////////////////////////////

#include "PatternOptimizationJob.hh"
#include <isosurface_inflator/WireMesh.hh>
#include <MeshFEM/Future.hh>
#include <MeshFEM/StringUtils.hh>
#include <stdexcept>
#include <fstream>
#include <string>
#include <memory>

using namespace std;

namespace PatternOptimization {

// -----------------------------------------------------------------------------

namespace {

#if 0
void parseVector(const ptree &pt, vector<Real> &v) {
    v.clear();
    for (const auto &val : pt) {
        if (!val.first.empty()) throw runtime_error("Failed to parse vector");
        v.push_back(val.second.get_value<Real>());
    }
}

void parseVector(const ptree &pt, vector<bool> &v) {
    v.clear();
    for (const auto &val : pt) {
        if (!val.first.empty()) throw runtime_error("Failed to parse vector");

        if (val.second.get_value<int>() == 1) {
            v.push_back(true);
        }
        else {
            v.push_back(false);
        }
    }
}
#endif

template<typename T>
void read_vector_if_present(const nlohmann::json &entry, const std::string &key, std::vector<T> &x) {
    if (entry.count(key)) {
        x = entry[key].get<std::vector<T>>();
    }
}

} // anonymous namespace

// -----------------------------------------------------------------------------

size_t JobBase::numParams() const {
    if (paramsMask.empty()) {
        return initialParams.size();
    }
    else {
        assert(paramsMask.size() == initialParams.size());
        return std::count(paramsMask.begin(), paramsMask.begin()+paramsMask.size(), false);
    }
}

// -----------------------------------------------------------------------------

// Verifies that the correct number of parameters were specified in the job
// (must match inflator). For non-parametric inflators (like the
// BoundaryPerturbationInflator) an incorrect number of parameters is
// tolerated with a warning, in which case the initial parameters are taken
// to be all zero.
// Returns the (possibly modified) initial parameters
std::vector<Real> JobBase::validatedInitialParams(const InflatorBase &inflator) const {
    std::vector<Real> params(initialParams);
    // Set params to the default if they're omitted in the job file
    size_t np;
    if (numParams() == 0) {
        params = inflator.defaultParameters();
        np = params.size();
    }
    else {
        np = numParams();
    }
    //const size_t np = params.size();
    // Allow non-parametric inflator to ignore initialParams (if size mismatch)
    if (!inflator.isParametric()) {
        if (np != inflator.numParameters()) {
            std::cerr << "WARNING: ignoring incorrectly-sized initial parameters for non-parametric inflator." << std::endl;
            params.assign(inflator.numParameters(), 0.0);
        }
        return params;
    }
    if (np != inflator.numParameters()) {
        for (size_t i = 0; i < inflator.numParameters(); ++i)
            std::cerr << "param " << i
                      << " role: " << parameterTypeString(inflator.parameterType(i))
                      << std::endl;
        std::cerr <<  "Inflator was expecting " << inflator.numParameters() << " parameters, but input contained only "
                  << np << " values" << std::endl;
        throw std::runtime_error("Invalid number of parameters.");
    }

    for (const auto &boundEntry : varLowerBounds) {
        if (boundEntry.first > params.size())
            std::cerr << "WARNING: bound on nonexistent variable" << std::endl;
    }

    for (size_t p = 0; p < params.size(); ++p) {
        if (varLowerBounds.count(p)) {
            if (params[p] < varLowerBounds.at(p)) {
                std::cerr << "WARNING: param " << p << " clamped to lower bound." << std::endl;
                params[p] = varLowerBounds.at(p);
            }
            if (params[p] > varUpperBounds.at(p)) {
                std::cerr << "WARNING: param " << p << " clamped to upper bound." << std::endl;
                params[p] = varUpperBounds.at(p);
            }
        }
    }

    return params;
}

// -----------------------------------------------------------------------------

void JobBase::writeJobFile(const std::string &jobFile) const {
    std::ofstream os(jobFile);
    if (!os.is_open()) {
        throw std::runtime_error("Couldn't open output job file " + jobFile);
    }
    writeJobFile(os);
}

////////////////////////////////////////////////////////////////////////////////

template<size_t _N>
void Job<_N>::writeJobFile(std::ostream &os) const {
    os << std::setprecision(19);
    os << getJson().dump(4) << std::endl;
    return;

#if 0
    os << "{" << std::endl
       << "\t\"dim\": " << _N << "," << std::endl
       << "\t\"target\": " << targetMaterial << "," << std::endl;
    if (targetVolume)
        os << "\t\"targetVolume\": " << *targetVolume << "," << std::endl;
    if (initialParams.size()) {
        os << "\t\"initial_params\": [";
        for (size_t i = 0; i < initialParams.size(); ++i)
            os << (i ? ", " : "") << std::setprecision(10) << initialParams[i];
        os << "]," << std::endl;
    }
    if (paramsMask.size()) {
        os << "\t\"paramsMask\": [";
        for (size_t i = 0; i < paramsMask.size(); ++i)
            os << (i ? ", " : "") << paramsMask[i];
        os << "]," << std::endl;
    }

    if (trueParams.size() == initialParams.size()) {
        os << "\t\"# true_params\": [";
        for (size_t i = 0; i < trueParams.size(); ++i)
            os << (i ? ", " : "") << trueParams[i];
        os << "]," << std::endl;
    }

    if (parameterConstraints.size()) {
        os << "\t\"paramConstraints\": [";
        bool first = true;
        for (const std::string &pc : parameterConstraints) {
            if (!first) os << ",";
            first = false;
            os << std::endl << "\t\t\"" << pc << "\"";
        }
        os << "\t]," << std::endl;
    }

    if (varLowerBounds.size() + varUpperBounds.size()) {
        os << "\t\"bounds\": [";
        bool first = true;
        for (size_t p = 0; p < initialParams.size(); ++p) {
            if (varLowerBounds.count(p) + varUpperBounds.count(p) == 0)
                continue;
            if (!first) os << ",";
            first = false;
            os << std::endl << "\t\t{\"var\": " << p;
            if (varLowerBounds.count(p)) os << ", \"lower\": " << varLowerBounds.at(p);
            if (varUpperBounds.count(p)) os << ", \"upper\": " << varUpperBounds.at(p);
            os << "}";
        }
        os << std::endl << "\t]," << std::endl;
    }

    if (numberCustomTypes > 0)
        os << "\t\"custom1Bounds\": [" << custom1Bounds[0] << ", " << custom1Bounds[1] << "]," << std::endl;
    if (numberCustomTypes > 1)
        os << "\t\"custom2Bounds\": [" << custom2Bounds[0] << ", " << custom2Bounds[1] << "]," << std::endl;
    if (numberCustomTypes > 2)
        os << "\t\"custom3Bounds\": [" << custom3Bounds[0] << ", " << custom3Bounds[1] << "]," << std::endl;
    if (numberCustomTypes > 3)
        os << "\t\"custom4Bounds\": [" << custom4Bounds[0] << ", " << custom4Bounds[1] << "]," << std::endl;
    if (numberCustomTypes > 4)
        os << "\t\"custom5Bounds\": [" << custom5Bounds[0] << ", " << custom5Bounds[1] << "]," << std::endl;
    if (numberCustomTypes > 5)
        os << "\t\"custom6Bounds\": [" << custom6Bounds[0] << ", " << custom6Bounds[1] << "]," << std::endl;
    if (numberCustomTypes > 6)
        os << "\t\"custom7Bounds\": [" << custom7Bounds[0] << ", " << custom7Bounds[1] << "]," << std::endl;
    if (numberCustomTypes > 7)
        os << "\t\"custom8Bounds\": [" << custom8Bounds[0] << ", " << custom8Bounds[1] << "]," << std::endl;

    os << "\t\"radiusBounds\": [" << radiusBounds[0] << ", " << radiusBounds[1] << "]," << std::endl
       << "\t\"translationBounds\": [" << translationBounds[0] << ", " << translationBounds[1] << "]," << std::endl
       << "\t\"blendingBounds\": [" << blendingBounds[0] << ", " << blendingBounds[1] << "]" << std::endl;

    os << "}" << std::endl;
#endif
}

template<size_t _N>
nlohmann::json Job<_N>::getJson() const {
    using json = nlohmann::json;

    json job;

    job["dim"] = _N;
    job["target"] = targetMaterial.getJson();
    if (targetVolume) {
        job["targetVolume"] = *targetVolume;
    }
    if (initialParams.size()) {
        job["initial_params"] = initialParams;
    }
    if (targetParams.size()) {
        job["target_params"] = targetParams;
    }
    if (paramsMask.size()) {
        job["paramsMask"] = paramsMask;
    }
    if (trueParams.size() == initialParams.size()) {
        job["true_params"] = trueParams;
    }
    if (parameterConstraints.size()) {
        job["paramConstraints"] = parameterConstraints;
    }
    if (varLowerBounds.size() + varUpperBounds.size()) {
        job["bounds"] = json::array();
        for (size_t p = 0; p < initialParams.size(); ++p) {
            if (varLowerBounds.count(p) + varUpperBounds.count(p) == 0) {
                continue;
            }
            json entry = {
                {"var", p},
                {"lower", varLowerBounds.at(p)},
                {"upper", varUpperBounds.at(p)}
            };
            job["bounds"].push_back(entry);
        }
    }

    if (numberCustomTypes > 0) { job["custom1Bounds"] = custom1Bounds; }
    if (numberCustomTypes > 1) { job["custom2Bounds"] = custom2Bounds; }
    if (numberCustomTypes > 2) { job["custom3Bounds"] = custom3Bounds; }
    if (numberCustomTypes > 3) { job["custom4Bounds"] = custom4Bounds; }
    if (numberCustomTypes > 4) { job["custom5Bounds"] = custom5Bounds; }
    if (numberCustomTypes > 5) { job["custom6Bounds"] = custom6Bounds; }
    if (numberCustomTypes > 6) { job["custom7Bounds"] = custom7Bounds; }
    if (numberCustomTypes > 7) { job["custom8Bounds"] = custom8Bounds; }

    job["radiusBounds"] = radiusBounds;
    job["translationBounds"] = translationBounds;
    job["blendingBounds"] = blendingBounds;

    return job;
}

////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<JobBase> parseJobFile(const string &jobFile) {
    ifstream is(jobFile);
    if (!is.is_open()) {
        throw runtime_error("Couldn't open job file " + jobFile);
    }

    nlohmann::json config;
    is >> config;
    return jobFromJson(config);

#if 0
    ptree pt;
    read_json(is, pt);

    size_t dim = pt.get<size_t>("dim");
    auto materialSpec = pt.get_child("target");
    auto radiusBounds = pt.get_child("radiusBounds");
    auto translationBounds = pt.get_child("translationBounds");

    std::unique_ptr<JobBase> job;
    if (dim == 2) {
        auto job2D = Future::make_unique<Job<2>>();
        job2D->targetMaterial.setFromPTree(materialSpec);
        job = std::move(job2D);
    }
    else if (dim == 3) {
        auto job3D = Future::make_unique<Job<3>>();
        job3D->targetMaterial.setFromPTree(materialSpec);
        job = std::move(job3D);
    }
    else throw runtime_error("Invalid dimension.");

    if (pt.count("initial_params")) {
        auto paramVals = pt.get_child("initial_params");
        parseVector(paramVals, job->initialParams);
    }

    parseVector(radiusBounds, job->radiusBounds);
    parseVector(translationBounds, job->translationBounds);

    if (pt.count("paramsMask"))
        parseVector(pt.get_child("paramsMask"), job->paramsMask);

    if (pt.count("blendingBounds"))
        parseVector(pt.get_child("blendingBounds"), job->blendingBounds);

    if (pt.count("metaBounds"))
        parseVector(pt.get_child("metaBounds"), job->metaBounds);

    if (pt.count("custom1Bounds"))
        parseVector(pt.get_child("custom1Bounds"), job->custom1Bounds);

    if (pt.count("custom2Bounds"))
        parseVector(pt.get_child("custom2Bounds"), job->custom2Bounds);

    if (pt.count("custom3Bounds"))
        parseVector(pt.get_child("custom3Bounds"), job->custom3Bounds);

    if (pt.count("custom4Bounds"))
        parseVector(pt.get_child("custom4Bounds"), job->custom4Bounds);

    if (pt.count("custom5Bounds"))
        parseVector(pt.get_child("custom5Bounds"), job->custom5Bounds);

    if (pt.count("custom6Bounds"))
        parseVector(pt.get_child("custom6Bounds"), job->custom6Bounds);

    if (pt.count("custom7Bounds"))
        parseVector(pt.get_child("custom7Bounds"), job->custom7Bounds);

    if (pt.count("custom8Bounds"))
        parseVector(pt.get_child("custom8Bounds"), job->custom8Bounds);

    if (pt.count("paramConstraints") != 0) {
        auto constraints = pt.get_child("paramConstraints");
        for (const auto &val : constraints) {
            if (!val.first.empty()) throw runtime_error("Failed to read constraints");
            job->parameterConstraints.emplace_back(val.second.get_value<string>());
        }
    }

    // individual variable overriding radius/translation bounds
    if (pt.count("bounds") !=  0) {
        try {
            auto bounds = pt.get_child("bounds");
            for (const auto &bound : bounds) {
                size_t var = bound.second.get<size_t>("var");
                job->varLowerBounds[var] = bound.second.get<Real>("lower");
                job->varUpperBounds[var] = bound.second.get<Real>("upper");
            }
        }
        catch (...) {
            throw std::runtime_error("Couldn't parse variable bounds.");
        }
    }

    job->targetVolume = pt.get_optional<double>("targetVolume");

    return job;
#endif
}

std::unique_ptr<JobBase> jobFromJson(const nlohmann::json &entry) {
    size_t dim = entry["dim"];
    auto materialSpec = entry["target"];

    std::unique_ptr<JobBase> job;
    if (dim == 2) {
        auto job2D = Future::make_unique<Job<2>>();
        job2D->targetMaterial.setFromJson(materialSpec);
        job = std::move(job2D);
    } else if (dim == 3) {
        auto job3D = Future::make_unique<Job<3>>();
        job3D->targetMaterial.setFromJson(materialSpec);
        job = std::move(job3D);
    } else {
        throw runtime_error("Invalid dimension.");
    }

    if (entry.count("initial_params")) {
        job->initialParams = entry["initial_params"].get<std::vector<Real>>();
    }

    if (entry.count("target_params")) {
        job->targetParams = entry["target_params"].get<std::vector<Real>>();
    }

    job->radiusBounds = entry["radiusBounds"].get<std::vector<Real>>();;
    job->translationBounds = entry["translationBounds"].get<std::vector<Real>>();;

    read_vector_if_present(entry, "paramsMask", job->paramsMask);
    read_vector_if_present(entry, "blendingBounds", job->blendingBounds);
    read_vector_if_present(entry, "metaBounds", job->metaBounds);
    read_vector_if_present(entry, "custom1Bounds", job->custom1Bounds);
    read_vector_if_present(entry, "custom2Bounds", job->custom2Bounds);
    read_vector_if_present(entry, "custom3Bounds", job->custom3Bounds);
    read_vector_if_present(entry, "custom4Bounds", job->custom4Bounds);
    read_vector_if_present(entry, "custom5Bounds", job->custom5Bounds);
    read_vector_if_present(entry, "custom6Bounds", job->custom6Bounds);
    read_vector_if_present(entry, "custom7Bounds", job->custom7Bounds);
    read_vector_if_present(entry, "custom8Bounds", job->custom8Bounds);
    read_vector_if_present(entry, "paramConstraints", job->parameterConstraints);

    // individual variable overriding radius/translation bounds
    if (entry.count("bounds") !=  0) {
        try {
            for (const auto &bound : entry["bounds"]) {
                size_t var = bound["var"];
                job->varLowerBounds[var] = bound["lower"];
                job->varUpperBounds[var] = bound["upper"];
            }
        } catch (...) {
            throw std::runtime_error("Couldn't parse variable bounds.");
        }
    }

    if (entry.count("targetVolume")) {
        job->targetVolume = entry["targetVolume"].get<Real>();
    }

    return job;
}

////////////////////////////////////////////////////////////////////////////////

template<typename WireMesh>
std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh &wm, const nlohmann::json &new_args) {
    // Default args
    nlohmann::json args = R"({
        "dim": 2,
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

    // Override with user-specified args
    args.update(new_args);

    auto targetModuli = args["elasticityTensor"];

    std::unique_ptr<JobBase> job;
    if (args["dim"] == 2) {
        auto job2D = Future::make_unique<PatternOptimization::Job<2>>();
        job2D->targetMaterial.setIsotropic(targetModuli[0], targetModuli[1]);
        job = move(job2D);
    } else if (args["dim"] == 3) {
        auto job3D = Future::make_unique<PatternOptimization::Job<3>>();
        job3D->targetMaterial.setIsotropic(targetModuli[0], targetModuli[1]);
        job = move(job3D);
    } else {
        throw std::runtime_error("Invalid dimension.");
    }

    if (!args["offsetBounds"].empty()) {
        std::vector<double> offsetBds = args["offsetBounds"];
        auto defaultPositions = wm.defaultPositionParams();

        for (size_t p = 0; p < defaultPositions.size(); ++p) {
            // Position parameters should be first in the isosurface inflator
            assert(wm.isPositionParam(p));
            double lowerBound;
            double upperBound;
            if (args["limitedOffset"]) {
                lowerBound = (defaultPositions[p] + offsetBds[0]) > 0 ? (defaultPositions[p] + offsetBds[0]) : 0;
                upperBound = (defaultPositions[p] + offsetBds[1]) < 1 ? (defaultPositions[p] + offsetBds[1]) : 1;
            } else {
                lowerBound = defaultPositions[p] + offsetBds[0];
                upperBound = defaultPositions[p] + offsetBds[1];
            }

            job->varLowerBounds.emplace(p, lowerBound);
            job->varUpperBounds.emplace(p, upperBound);
        }
    }

    // Set translation bounds, which will be ignored by the optimizer if
    // offsetBounds introduced per-variable bounds above
    job->translationBounds = {0.1, 0.8};
    if (!args["translationBounds"].empty()) {
        job->translationBounds = args["translationBounds"].get<std::vector<double>>();
    }

    job->radiusBounds = args["radiusBounds"].get<std::vector<double>>();
    job->blendingBounds = args["blendingBounds"].get<std::vector<double>>();

    if (!args["parameterConstraints"].empty()) {
        job->parameterConstraints = MeshFEM::split(args["parameterConstraints"], ";");
    }

    job->initialParams = wm.defaultParameters(args["defaultThickness"]);
    if (!args["initialParams"].empty()) {
        job->initialParams = args["initialParams"].get<std::vector<double>>();
    }

    return job;
}

////////////////////////////////////////////////////////////////////////////////
// Explicit template specialization
////////////////////////////////////////////////////////////////////////////////

template class Job<2>;
template class Job<3>;

template std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh<Symmetry::Cubic<>> &wm, const nlohmann::json &args);
template std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh<Symmetry::Orthotropic<>> &wm, const nlohmann::json &args);
template std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh<Symmetry::Diagonal<>> &wm, const nlohmann::json &args);
template std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh<Symmetry::Square<>> &wm, const nlohmann::json &args);
template std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh<Symmetry::TriplyPeriodic<>> &wm, const nlohmann::json &args);
template std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh<Symmetry::DoublyPeriodic<>> &wm, const nlohmann::json &args);
template std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh<Symmetry::NonPeriodic<>> &wm, const nlohmann::json &args);

////////////////////////////////////////////////////////////////////////////////

} // namespace PatternOptimization
