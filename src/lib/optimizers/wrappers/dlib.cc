#include "dlib.hh"

#if HAS_DLIB
#include <dlib/optimization.h>
#include <vector>

using dlib_vector = dlib::matrix<double, 0, 1>;
using namespace PatternOptimization;

// Evaluates the objective by inflating the wire mesh and homogenizing.
// The iterate stored internally also knows how to evaluate the gradient
// efficiently, so our GradientEvaluator below just accesses it.
struct DLibObjectiveEvaluator {
    DLibObjectiveEvaluator(IterateManagerBase &imanager) : m_imanager(imanager) { };

    double operator()(const dlib_vector &x) const {
        std::vector<Real> x_vec(m_imanager.numParameters());
        for (size_t p = 0; p < x_vec.size(); ++p)
            x_vec[p] = x(p);

        auto &it = m_imanager.get(x_vec.size(), &x_vec[0]);
        return it.evaluate();
    }
private:
    // Iterate manager is mutable so that operator() can be const as dlib requires
    IterateManagerBase &m_imanager;
};

// Extracts gradient from the iterate constructed by DLibObjectiveEvaluator.
struct DLibGradientEvaluator {
    DLibGradientEvaluator(const IterateManagerBase &imanager) : m_imanager(imanager) { }

    dlib_vector operator()(const dlib_vector &x) const {
        std::vector<Real> x_vec(m_imanager.numParameters());
        for (size_t p = 0; p < x_vec.size(); ++p)
            x_vec[p] = x(p);
        if (m_imanager.get().paramsDiffer(x_vec.size(), &x_vec[0]))
            throw std::runtime_error("Objective must be evaluated first");
        ScalarField<Real> gradp(m_imanager.get().gradp());

        dlib_vector result(x_vec.size());
        for (size_t p = 0; p < x_vec.size(); ++p)
            result(p) = gradp[p];
        return result;
    }
private:
    const IterateManagerBase &m_imanager;
};

// Hack to get notified at the end of each iteration: subclass the stop
// strategy.
class ReportingStopStrategy : public dlib::objective_delta_stop_strategy {
    typedef dlib::objective_delta_stop_strategy Base;
public:
    ReportingStopStrategy(double min_delta, unsigned long max_iter,
                          IterateManagerBase &imanager, const std::string &outPath)
        : Base(min_delta, max_iter), m_iter(0), m_im(imanager) {
        m_im.setOutPath(outPath);
    }

    template <typename T>
    bool should_continue_search(const T& x, const double funct_value,
        const T& funct_derivative) {

        std::vector<double> params(x.nr());
        for (size_t i=0; i<params.size(); i++) {
            params[i] = x(i);
        }

        m_im.updateAndReport(params);

        ++m_iter;
        return Base::should_continue_search(x, funct_value, funct_derivative);
    }

private:
    size_t m_iter;
    IterateManagerBase &m_im;
};

void optimize_dlib_bfgs(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    if (!im.isParametric()) throw std::runtime_error("BFGS only works with parametric inflators");
    DLibObjectiveEvaluator obj(im);
    DLibGradientEvaluator grad(im);

    size_t nParams = im.numParameters();
    // convert initial parameter vector
    dlib_vector optParams(nParams);
    for (size_t p = 0; p < nParams; ++p)
        optParams(p) = params[p];

    dlib_vector lowerBounds(nParams), upperBounds(nParams);
    for (size_t p = 0; p < nParams; ++p) {
        if (!(bds.hasLowerBound.at(p) && bds.hasUpperBound.at(p)))
            throw std::runtime_error("Currently all parameters must be bounded for bfgs.");
        lowerBounds(p) = bds.lowerBound.at(p);
        upperBounds(p) = bds.upperBound.at(p);
    }

    size_t max_size = oconfig.lbfgs_memory;
    if (max_size == 0) {
        dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
                ReportingStopStrategy(1e-16, oconfig.niters, im, outPath),
                obj, grad, optParams, lowerBounds, upperBounds);
    }
    else {
        dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(max_size),
                ReportingStopStrategy(1e-16, oconfig.niters, im, outPath),
                obj, grad, optParams, lowerBounds, upperBounds);
    }

    // convert solution
    for (size_t p = 0; p < nParams; ++p)
        params[p] = optParams(p);
}

dlib_vector scalarFieldToDlibVector(ScalarField<Real> input) {
    dlib_vector result(input.size());
    for (size_t p = 0; p < input.size(); ++p)
        result(p) = input[p];
    return result;
}

ScalarField<Real> dlibVectorToScalarField(dlib_vector input) {
    ScalarField<Real> result(input.nr());
    for (size_t p = 0; p < input.nr(); ++p)
        result[p] = input(p);
    return result;
}

ScalarField<Real> computeBfgsDirection(ScalarField<Real> params, ScalarField<Real> grad, dlib::bfgs_search_strategy &strategy) {
    ScalarField<Real> result;

    dlib_vector xDlib = scalarFieldToDlibVector(params);
    dlib_vector gradDlib = scalarFieldToDlibVector(grad);
    dlib_vector directionDlib = strategy.get_next_direction(xDlib, 0.0, gradDlib);
    ScalarField<Real> direction = dlibVectorToScalarField(directionDlib);

    // Normalize according to grad norm
    //result = direction * (grad.norm() / direction.norm());

    return direction;
}

void optimize_dlib_custom_bfgs(ScalarField<Real> &params,
                        const PatternOptimization::BoundConstraints &bds,
                        PatternOptimization::IterateManagerBase &im,
                        const PatternOptimization::OptimizerConfig &oconfig,
                        const std::string &outPath)
{
    double step_size = oconfig.gd_step;
    int success_iterations = 0;
    int failure_iterations = 0;
    std::vector<Real> oldParamsCopy;
    std::vector<Real> directionCopy;
    double val = std::numeric_limits<Real>::max();
    double bestVal = std::numeric_limits<Real>::max();
    double previousVal = std::numeric_limits<Real>::max();
    double minimumStep = 1e-15;
    size_t maxTries = 100;
    size_t tries = 0;

    im.setOutPath(outPath);

    dlib::bfgs_search_strategy strategy = dlib::bfgs_search_strategy();

    for (size_t i = 0; i < oconfig.niters; ++i) {
        //std::cout << "[Dlib] New bfgs iteration" << std::endl;

        if (step_size < minimumStep) {
            if (tries < maxTries) {
                step_size = oconfig.gd_step;
                tries++;
            }
            else {
                return;
            }
        }

        try {
            std::vector<Real> paramsCopy;
            params.getFlattened(paramsCopy);

            auto &it = im.get(params.size(), params.data());

            ScalarField<Real> direction = computeBfgsDirection(params, -it.steepestDescent(), strategy);
            //std::cout << "[Dlib] Direction norm: " << direction.norm() << std::endl;

            val = it.evaluate();
            if (val <= (previousVal + 1e-10)) {
                //std::cout << "[Dlib] Previous good val: " << previousVal << std::endl;
                previousVal = val;

                // Save direction and params
                direction.getFlattened(directionCopy);
                params.getFlattened(oldParamsCopy);

                //std::cout << "[Dlib] New good val: " << previousVal << std::endl;

                if (val < bestVal) {
                    bestVal = val;
                    //std::cout << "[Dlib] Best iteration: " << i << std::endl;
                    //std::cout << "[Dlib] Best value: " << val << std::endl;
                }
                success_iterations++;
                failure_iterations=0;
            }
            else {
                success_iterations=0;
                failure_iterations++;
            }

            // Report
            im.updateAndReport(paramsCopy);

            ScalarField<Real> perturbation = direction * step_size;
            params += perturbation;

            // Apply bound constraints
            for (size_t p = 0; p < im.numParameters(); ++p) {
                if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
                if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
            }
        }
        catch (...) {
            //std::cout << "[Dlib] Exploded!" << std::endl;
            step_size = step_size / 2;
            success_iterations = 0;
            failure_iterations = 0;

            // Recover!
            ScalarField<Real> oldParams(oldParamsCopy);
            ScalarField<Real> direction(directionCopy);
            ScalarField<Real> newPerturbation = direction * step_size;
            params = oldParams + newPerturbation;

            // Apply bound constraints
            for (size_t p = 0; p < im.numParameters(); ++p) {
                if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
                if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
            }
        }

        if (success_iterations >= 1) {
            step_size *= 2.0;
            success_iterations = 0;
            failure_iterations = 0;
        }
        else if (failure_iterations >= 1) {
            //std::cout << "[Dlib] Failure" << std::endl;
            step_size = step_size / 2;
            success_iterations = 0;
            failure_iterations = 0;

            // Recover!
            ScalarField<Real> oldParams(oldParamsCopy);
            ScalarField<Real> direction(directionCopy);
            ScalarField<Real> newPerturbation = direction * step_size;
            params = oldParams + newPerturbation;

            // Apply bound constraints
            for (size_t p = 0; p < im.numParameters(); ++p) {
                if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
                if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
            }
        }

        //std::cout << "[Dlib] Step size is now: " << step_size << std::endl;
        //std::cout << "[Dlib] End of iteration val: " << previousVal << std::endl;

    }
}

#else // !HAS_DLIB

void optimize_dlib_bfgs(ScalarField<Real> &/* params */,
        const PatternOptimization::BoundConstraints &/* bds */,
        PatternOptimization::IterateManagerBase &/* im */,
        const PatternOptimization::OptimizerConfig &/* oconfig */,
        const std::string &/* outPath */)
{
    throw std::runtime_error("Built without DLib");
}

void optimize_dlib_custom_bfgs(ScalarField<Real> &/* params */,
        const PatternOptimization::BoundConstraints &/* bds */,
        PatternOptimization::IterateManagerBase &/* im */,
        const PatternOptimization::OptimizerConfig &/* oconfig */,
        const std::string &/* outPath */)
{
    throw std::runtime_error("Built without DLib");
}

#endif // HAS_DLIB
