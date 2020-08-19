////////////////////////////////////////////////////////////////////////////////
// PatternOptimizationJob.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Reads in JSON files specifying an optimization job. These job
//      specifications include a target material and the initial pattern
//      parameters from which to start optimization.
//      Format:
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
#ifndef PATTERNOPTIMIZATIONJOB_HH
#define PATTERNOPTIMIZATIONJOB_HH

#include <MeshFEM/Materials.hh>
#include <inflators/Inflator.hh>
#include <nonstd/optional.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace PatternOptimization  {

class JobBase {

public:
    virtual ~JobBase() { }

    size_t numParams() const;

    // Verifies that the correct number of parameters were specified in the job
    // (must match inflator). For non-parametric inflators (like the
    // BoundaryPerturbationInflator) an incorrect number of parameters is
    // tolerated with a warning, in which case the initial parameters are taken
    // to be all zero.
    // Returns the (possibly modified) initial parameters
    std::vector<Real> validatedInitialParams(const InflatorBase &inflator) const;

    // Export job file as a json dictionary
    virtual nlohmann::json getJson() const = 0;

    void writeJobFile(const std::string &jobFile) const;

    virtual void writeJobFile(std::ostream &os) const = 0;

public:
    std::vector<Real> initialParams, radiusBounds, translationBounds;
    std::vector<bool> paramsMask;
    std::vector<std::string> metaParams;
    std::vector<Real> blendingBounds = { 10.0, 100.0 };
    std::vector<Real> metaBounds = { 0.01, 0.99 };
    std::vector<Real> custom1Bounds = { 0.01, 0.99 };
    std::vector<Real> custom2Bounds = { 0.01, 0.99 };
    std::vector<Real> custom3Bounds = { 0.01, 0.99 };
    std::vector<Real> custom4Bounds = { 0.01, 0.99 };
    std::vector<Real> custom5Bounds = { 0.01, 0.99 };
    std::vector<Real> custom6Bounds = { 0.01, 0.99 };
    std::vector<Real> custom7Bounds = { 0.01, 0.99 };
    std::vector<Real> custom8Bounds = { 0.01, 0.99 };
    // The ground-truth parameters can be stored here--they are written to the
    // job file for reference.
    std::vector<Real> trueParams;
    std::vector<std::string> parameterConstraints;
    std::map<size_t, Real> varLowerBounds, varUpperBounds;
    nonstd::optional<Real> targetVolume;
    std::vector<Real> targetParams;
    size_t numberCustomTypes = 0;
};

// -----------------------------------------------------------------------------

template<size_t _N>
class Job : public JobBase {
public:
    virtual ~Job() = default;

    // Export job file as a json dictionary
    virtual nlohmann::json getJson() const override;

    virtual void writeJobFile(std::ostream &os) const override;

    Materials::Constant<_N> targetMaterial;
};

// -----------------------------------------------------------------------------

// Create new job based on a previously exported one
std::unique_ptr<JobBase> jobFromJson(const nlohmann::json &job);
std::unique_ptr<JobBase> parseJobFile(const std::string &jobFile);

// Create new job for a given wire mesh, and additional arguments
template<typename WireMesh>
std::unique_ptr<JobBase> jobFromWireMesh(const WireMesh &wm, const nlohmann::json &args);

} // namespace PatternOptimization

#endif /* end of include guard: PATTERNOPTIMIZATIONJOB_HH */
