#ifndef OPTIMIZERCONFIG_HH
#define OPTIMIZERCONFIG_HH

#include <boost/optional.hpp>
#include <limits>

namespace PatternOptimization {

struct OptimizerConfig {
    size_t niters = std::numeric_limits<int>::max(); // not size_t's max b/c some optimizers use signed ints
    size_t lbfgs_memory = 0; // use full BFGS by default
    double gd_step = 0.0001;

    // Tolerance for the tensor fitting objective
    boost::optional<double> tensor_fit_tolerance;
    boost::optional<double> constraint_tolerance;
};

}


#endif /* end of include guard: OPTIMIZERCONFIG_HH */
