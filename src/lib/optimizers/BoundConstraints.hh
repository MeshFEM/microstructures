#ifndef BOUNDCONSTRAINTS_HH
#define BOUNDCONSTRAINTS_HH

#include <inflators/Inflator.hh>
#include <vector>
#include <map>

namespace PatternOptimization {

struct BoundConstraints {
    // lb, ub: manually configured lower and upper bounds. All other bounds are
    // set using the defaults in radiusBounds[01] and translationBounds[01],
    // selected based on parameter type as determined by inflator.
    BoundConstraints(const InflatorBase &inflator,
                     const std::vector<Real> &radiusBounds,
                     const std::vector<Real> &translationBounds,
                     const std::vector<Real> &blendingBounds,
                     const std::vector<Real> &metaBounds,
                     const std::vector<Real> &custom1Bounds, // custom types, allowing us to define bounds for types different than radius, offset and blending
                     const std::vector<Real> &custom2Bounds,
                     const std::vector<Real> &custom3Bounds,
                     const std::vector<Real> &custom4Bounds,
                     const std::vector<Real> &custom5Bounds,
                     const std::vector<Real> &custom6Bounds,
                     const std::vector<Real> &custom7Bounds,
                     const std::vector<Real> &custom8Bounds,
                     const std::map<size_t, Real> &lb,
                     const std::map<size_t, Real> &ub) {
        size_t nparams = inflator.numParameters();

        // All params get bounds by default.
        hasLowerBound.assign(nparams, true);
        hasUpperBound.assign(nparams, true);
        lowerBound.assign(nparams, 0.0);
        upperBound.assign(nparams, 0.0);

        // Explicitly specified bounds override default type-based bounds.
        // Set lower bounds
        for (size_t p = 0; p < nparams; ++p) {
            if (lb.count(p)) { lowerBound.at(p) = lb.at(p); continue; }
            switch (inflator.parameterType(p)) {
                case ParameterType::Thickness: lowerBound.at(p) =      radiusBounds.at(0); break;
                case ParameterType::Offset:    lowerBound.at(p) = translationBounds.at(0); break;
                case ParameterType::Blending:  lowerBound.at(p) =    blendingBounds.at(0); break;
                case ParameterType::Meta:      lowerBound.at(p) =        metaBounds.at(0); break;
                case ParameterType::Custom1:      lowerBound.at(p) =        custom1Bounds.at(0); break;
                case ParameterType::Custom2:      lowerBound.at(p) =        custom2Bounds.at(0); break;
                case ParameterType::Custom3:      lowerBound.at(p) =        custom3Bounds.at(0); break;
                case ParameterType::Custom4:      lowerBound.at(p) =        custom4Bounds.at(0); break;
                case ParameterType::Custom5:      lowerBound.at(p) =        custom5Bounds.at(0); break;
                case ParameterType::Custom6:      lowerBound.at(p) =        custom6Bounds.at(0); break;
                case ParameterType::Custom7:      lowerBound.at(p) =        custom7Bounds.at(0); break;
                case ParameterType::Custom8:      lowerBound.at(p) =        custom8Bounds.at(0); break;
                default: assert(false);
            }
        }

        for (size_t p = 0; p < nparams; ++p) {
            if (ub.count(p)) { upperBound.at(p) = ub.at(p); continue; }
            switch (inflator.parameterType(p)) {
                case ParameterType::Thickness: upperBound.at(p) =      radiusBounds.at(1); break;
                case ParameterType::Offset:    upperBound.at(p) = translationBounds.at(1); break;
                case ParameterType::Blending:  upperBound.at(p) =    blendingBounds.at(1); break;
                case ParameterType::Meta:      upperBound.at(p) =        metaBounds.at(1); break;
                case ParameterType::Custom1:      upperBound.at(p) =        custom1Bounds.at(1); break;
                case ParameterType::Custom2:      upperBound.at(p) =        custom2Bounds.at(1); break;
                case ParameterType::Custom3:      upperBound.at(p) =        custom3Bounds.at(1); break;
                case ParameterType::Custom4:      upperBound.at(p) =        custom4Bounds.at(1); break;
                case ParameterType::Custom5:      upperBound.at(p) =        custom5Bounds.at(1); break;
                case ParameterType::Custom6:      upperBound.at(p) =        custom6Bounds.at(1); break;
                case ParameterType::Custom7:      upperBound.at(p) =        custom7Bounds.at(1); break;
                case ParameterType::Custom8:      upperBound.at(p) =        custom8Bounds.at(1); break;
                default: assert(false);
            }
        }
    }

    std::vector<int> generateFilterMap(const std::vector<bool> &parametersMask) {
        size_t originalIdx = 0;
        std::vector<int> result;

        for (; originalIdx < parametersMask.size(); originalIdx++) {
            if (!parametersMask[originalIdx]) {
                result.push_back(originalIdx);
            }
        }

        return result;
    }

    BoundConstraints(const InflatorBase &inflator,
                     const std::vector<bool> &paramsMask,
                     const std::vector<Real> &radiusBounds,
                     const std::vector<Real> &translationBounds,
                     const std::vector<Real> &blendingBounds,
                     const std::vector<Real> &metaBounds,
                     const std::vector<Real> &custom1Bounds,
                     const std::vector<Real> &custom2Bounds,
                     const std::vector<Real> &custom3Bounds,
                     const std::vector<Real> &custom4Bounds,
                     const std::vector<Real> &custom5Bounds,
                     const std::vector<Real> &custom6Bounds,
                     const std::vector<Real> &custom7Bounds,
                     const std::vector<Real> &custom8Bounds,
                     const std::map<size_t, Real> &lb,
                     const std::map<size_t, Real> &ub) {
        size_t nparams = inflator.numParameters();

        std::vector<int> paramMap = generateFilterMap(paramsMask);

        assert(nparams == paramMap.size());

        // All params get bounds by default.
        hasLowerBound.assign(nparams, true);
        hasUpperBound.assign(nparams, true);
        lowerBound.assign(nparams, 0.0);
        upperBound.assign(nparams, 0.0);

        // Explicitly specified bounds override default type-based bounds.
        // Set lower bounds
        for (size_t p = 0; p < nparams; ++p) {
            if (lb.count(paramMap[p])) { lowerBound.at(p) = lb.at(paramMap[p]); continue; }
            switch (inflator.parameterType(p)) {
                case ParameterType::Thickness: lowerBound.at(p) =      radiusBounds.at(0); break;
                case ParameterType::Offset:    lowerBound.at(p) = translationBounds.at(0); break;
                case ParameterType::Blending:  lowerBound.at(p) =    blendingBounds.at(0); break;
                case ParameterType::Meta:      lowerBound.at(p) =        metaBounds.at(0); break;
                case ParameterType::Custom1:      lowerBound.at(p) =        custom1Bounds.at(0); break;
                case ParameterType::Custom2:      lowerBound.at(p) =        custom2Bounds.at(0); break;
                case ParameterType::Custom3:      lowerBound.at(p) =        custom3Bounds.at(0); break;
                case ParameterType::Custom4:      lowerBound.at(p) =        custom4Bounds.at(0); break;
                case ParameterType::Custom5:      lowerBound.at(p) =        custom5Bounds.at(0); break;
                case ParameterType::Custom6:      lowerBound.at(p) =        custom6Bounds.at(0); break;
                case ParameterType::Custom7:      lowerBound.at(p) =        custom7Bounds.at(0); break;
                case ParameterType::Custom8:      lowerBound.at(p) =        custom8Bounds.at(0); break;
                default: assert(false);
            }
        }

        for (size_t p = 0; p < nparams; ++p) {
            if (ub.count(paramMap[p])) { upperBound.at(p) = ub.at(paramMap[p]); continue; }
            switch (inflator.parameterType(p)) {
                case ParameterType::Thickness: upperBound.at(p) =      radiusBounds.at(1); break;
                case ParameterType::Offset:    upperBound.at(p) = translationBounds.at(1); break;
                case ParameterType::Blending:  upperBound.at(p) =    blendingBounds.at(1); break;
                case ParameterType::Meta:      upperBound.at(p) =        metaBounds.at(1); break;
                case ParameterType::Custom1:      upperBound.at(p) =        custom1Bounds.at(1); break;
                case ParameterType::Custom2:      upperBound.at(p) =        custom2Bounds.at(1); break;
                case ParameterType::Custom3:      upperBound.at(p) =        custom3Bounds.at(1); break;
                case ParameterType::Custom4:      upperBound.at(p) =        custom4Bounds.at(1); break;
                case ParameterType::Custom5:      upperBound.at(p) =        custom5Bounds.at(1); break;
                case ParameterType::Custom6:      upperBound.at(p) =        custom6Bounds.at(1); break;
                case ParameterType::Custom7:      upperBound.at(p) =        custom7Bounds.at(1); break;
                case ParameterType::Custom8:      upperBound.at(p) =        custom8Bounds.at(1); break;
                default: assert(false);
            }
        }
    }

    std::vector<Real> lowerBound, upperBound;
    std::vector<bool> hasLowerBound, hasUpperBound;
};

}

#endif /* end of include guard: BOUNDCONSTRAINTS_HH */
