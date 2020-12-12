#include "nlopt.hh"

#if HAS_NLOPT

#include <nlopt.hpp>
#include <limits>

#include <vector>
#include <cassert>

#include "../Constraint.hh"

#define EQ_TOL 1e-2
#define INEQ_TOL 1e-14
#define MAX_ITERATIONS_WITHOUT_IMPROVEMENT 100

using namespace PatternOptimization;

std::string nloptPrettyPrint(nlopt::result res) {
    switch (res) {
        case NLOPT_FAILURE: return "NLOPT_FAILURE";
        case NLOPT_INVALID_ARGS: return "NLOPT_INVALID_ARGS";
        case NLOPT_OUT_OF_MEMORY: return "NLOPT_OUT_OF_MEMORY";
        case NLOPT_ROUNDOFF_LIMITED: return "NLOPT_ROUNDOFF_LIMITED";
        case NLOPT_FORCED_STOP: return "NLOPT_FORCED_STOP";
        case NLOPT_SUCCESS: return "NLOPT_SUCCESS";
        case NLOPT_STOPVAL_REACHED: return "NLOPT_STOPVAL_REACHED";
        case NLOPT_FTOL_REACHED: return "NLOPT_FTOL_REACHED";
        case NLOPT_XTOL_REACHED: return "NLOPT_XTOL_REACHED";
        case NLOPT_MAXEVAL_REACHED: return "NLOPT_MAXEVAL_REACHED";
        case NLOPT_MAXTIME_REACHED: return "NLOPT_MAXTIME_REACHED";
        default: return "UNDEF";
    }
}

struct NLOptState {
    NLOptState(IterateManagerBase &im, nlopt::opt &opt) : m_im(im), opt(opt) { }

    void manualTerminationCheck(const PatternOptimization::IterateBase &it, Real /* costVal */) const {
        if (tensor_fit_tolerance) {
            // When printability constraints are present, don't allow early
            // termination unless they are satisfied (maximum violation is small)
            if (it.hasConstraint("Printability")) {
                const auto &pconstraint = it.evaluatedConstraint("Printability");
                if (pconstraint.values.max() > 1e-16) return;
            }

            double relFrobDistSq = std::numeric_limits<double>::max();
            try {
                relFrobDistSq = it.evaluateNormalized("JS");
            }
            catch(...) { std::cerr << "ERROR: no tensor_fit_tolerance requires a JS objective term" << std::endl; }
            if (relFrobDistSq < *tensor_fit_tolerance)
                opt.force_stop();
        }

        if (m_im.numberIterations() > m_im.bestIteration() + MAX_ITERATIONS_WITHOUT_IMPROVEMENT) {
            std::cout << "Stopping nlopt cause solution did not improve after " << MAX_ITERATIONS_WITHOUT_IMPROVEMENT << std::endl;
            opt.force_stop();
        }
    }

    boost::optional<double> tensor_fit_tolerance;
    IterateManagerBase &m_im;
    nlopt::opt &opt;
};

struct NLOptConstraintEvaluator {
    NLOptConstraintEvaluator(NLOptState &s, size_t idx)
        : state(s), constraintIndex(idx) { }

    NLOptState &state;
    size_t constraintIndex;
};

double costFunc(const std::vector<double> &x, std::vector<double> &grad, void *optStateVoid) {
    auto optState = reinterpret_cast<NLOptState *>(optStateVoid);
    assert(optState);

    try {
        auto &it = optState->m_im.get(x.size(), &x[0]);

        Real val = it.evaluate();

        if (!grad.empty()) {
            assert(grad.size() == x.size());
            ScalarField<Real> gp = it.gradp();
            assert(gp.domainSize() == grad.size());
            for (size_t p = 0; p < gp.domainSize(); ++p)
                grad[p] = gp[p];
        }

        optState->m_im.updateAndReport(x);

        if (it.shouldReport()) // only run termination check if this is a valid iterate (not an estimate)
            optState->manualTerminationCheck(it, val);

        double result = it.evaluate();

        return result;
    }
    catch (...) {
        std::cerr << "Iteration failed, setting value as maximum!" << std::endl;
        Real errorVal = std::numeric_limits<Real>::max();
        grad.resize(x.size());
        return errorVal;
    }
}

void constraintFunc(unsigned m, double *result, unsigned n, const double *x,
                    double *gradient, void *cevalVoid) {
    auto ceval = reinterpret_cast<NLOptConstraintEvaluator *>(cevalVoid);
    assert(ceval);

    auto &it = ceval->state.m_im.get(n, x);
    const EvaluatedConstraint &c = it.evaluatedConstraint(ceval->constraintIndex);
    assert(c.values.domainSize() == m);

    for (size_t i = 0; i < m; ++i)
        result[i] = c.values[i];

    // Unless the gradient is requested, we are done.
    if (gradient == NULL) return;

    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            gradient[i * n + j] = c.jacobian(i, j);
}

void optimize_nlopt_slsqp(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    nlopt::opt opt(nlopt::LD_SLSQP, params.domainSize());

    const size_t nParams = params.domainSize();
    for (size_t p = 0; p < nParams; ++p) {
        if (!(bds.hasLowerBound.at(p) && bds.hasUpperBound.at(p)))
            throw std::runtime_error("Currently all parameters must be bounded.");
    }

    opt.set_lower_bounds(bds.lowerBound);
    opt.set_upper_bounds(bds.upperBound);

    NLOptState state(im, opt);
    state.m_im.setOutPath(outPath);
    state.tensor_fit_tolerance = oconfig.tensor_fit_tolerance;

    opt.set_min_objective(costFunc, (void *) &state);

    // Must create iterate to query the constraints.
    // This iterate for the initial parameters will be reused by the first
    // optimization iteration.
    std::vector<std::unique_ptr<NLOptConstraintEvaluator>> cevals;
    auto &it = im.get(nParams, &params[0]);
    for (size_t i = 0; i < it.numConstraints(); ++i) {
        const auto &c = it.evaluatedConstraint(i);
        const size_t m = c.dimension();

        // Use lower tolerance on equality constraints
        Real tolVal = (c.type == ConstraintType::EQUALITY) ? EQ_TOL : INEQ_TOL;

        std::vector<Real> tol(m, tolVal);
        cevals.push_back(Future::make_unique<NLOptConstraintEvaluator>(state, i));
        auto ceval = cevals.back().get();
        switch (c.type) {
            case ConstraintType::  EQUALITY: opt.  add_equality_mconstraint(constraintFunc, (void *) ceval, tol); break;
            case ConstraintType::INEQUALITY: opt.add_inequality_mconstraint(constraintFunc, (void *) ceval, tol); break;
            default: throw std::runtime_error("Unexpected constraint type");
        }
    }

    opt.set_maxeval(oconfig.niters);
    opt.set_xtol_rel(1e-16);

    // nlopt's C++ interface works on std::vectors
    std::vector<double> x(nParams);
    for (size_t p = 0; p < nParams; ++p)
        x[p] = params[p];

    double minf = 0;
    try {
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "NLOpt terminated with result " << nloptPrettyPrint(result) << std::endl;
    }
    catch (nlopt::forced_stop) {
        std::cout << "NLOpt forced stop" << std::endl;
    }
    catch (std::exception &e) {
        std::cout << "NLOpt threw exception: " << e.what() << std::endl;
    }

    // Convert the solution back.
    for (size_t p = 0; p < nParams; ++p)
        params[p] = x[p];
}

void optimize_nlopt_lbfgs(ScalarField<Real> &params,
                          const PatternOptimization::BoundConstraints &bds,
                          PatternOptimization::IterateManagerBase &im,
                          const PatternOptimization::OptimizerConfig &oconfig,
                          const std::string &outPath)
{
    nlopt::opt opt(nlopt::LD_LBFGS, params.domainSize());

    // Bounding parameters
    const size_t nParams = params.domainSize();
    for (size_t p = 0; p < nParams; ++p) {
        if (!(bds.hasLowerBound.at(p) && bds.hasUpperBound.at(p))) {
            std::cerr << "Currently all parameters must be bounded." << std::endl;
            throw std::runtime_error("Currently all parameters must be bounded.");
        }
    }

    opt.set_lower_bounds(bds.lowerBound);
    opt.set_upper_bounds(bds.upperBound);

    // Setting custom state (where we can force terminate the optimization)
    NLOptState state(im, opt);
    state.m_im.setOutPath(outPath);
    state.tensor_fit_tolerance = oconfig.tensor_fit_tolerance;

    // Setting the energy
    opt.set_min_objective(costFunc, (void *) &state);

    // Must create iterate to query the constraints.
    // In the case of lbfgs, we expect no constraints
    std::vector<std::unique_ptr<NLOptConstraintEvaluator>> cevals;
    auto &it = im.get(nParams, &params[0]);
    if (it.numConstraints() > 0) {
        std::cerr << "Constraints should not be used with lbfgs." << std::endl;
        throw std::runtime_error("Constraints should not be used with lbfgs.");
    }

    opt.set_maxeval(oconfig.niters);
    opt.set_xtol_rel(1e-16);

    // nlopt's C++ interface works on std::vectors
    std::vector<double> x(nParams);
    for (size_t p = 0; p < nParams; ++p)
        x[p] = params[p];

    double minf = 0;
    try {
        nlopt::result result = opt.optimize(x, minf);
        std::cout << "NLOpt terminated with result " << nloptPrettyPrint(result) << std::endl;
    }
    catch (nlopt::forced_stop) {
        std::cout << "NLOpt forced stop" << std::endl;
    }
    catch (std::exception &e) {
        std::cout << "NLOpt threw exception: " << e.what() << std::endl;
    }

    // Convert the solution back.
    for (size_t p = 0; p < nParams; ++p)
        params[p] = x[p];
}

#else // !HAS_NLOPT

void optimize_nlopt_slsqp(ScalarField<Real> &/* params */,
        const PatternOptimization::BoundConstraints &/* bds */,
        PatternOptimization::IterateManagerBase &/* im */,
        const PatternOptimization::OptimizerConfig &/* oconfig */,
        const std::string &/* outPath */) {
    throw std::runtime_error("Built without NLopt");
}

void optimize_nlopt_lbfgs(ScalarField<Real> &/* params */,
        const PatternOptimization::BoundConstraints &/* bds */,
        PatternOptimization::IterateManagerBase &/* im */,
        const PatternOptimization::OptimizerConfig &/* oconfig */,
        const std::string &/* outPath */) {
    throw std::runtime_error("Built without NLopt");
}


#endif // HAS_NLOPT
