#include "ceres.hh"

using namespace PatternOptimization;

#if HAS_CERES
#include <ceres/ceres.h>
#include <glog/logging.h>

// Wraps (nonlinear) least squares-supporting objective terms
// (those with residuals and jacobians)
struct CeresCostWrapper : public ceres::CostFunction {
    typedef ceres::CostFunction Base;
    CeresCostWrapper(const std::string &quantity, IterateManagerBase &itManager)
        : m_quantity(quantity), m_itManager(itManager)
    {
        const auto currIt = m_itManager.getPtr();
        if (!currIt) throw std::runtime_error("Initial iterate must be inflated before building cost wrapper");

        // Validate that this is a NLLS objective term.
        const auto &term = dynamic_cast<const EvaluatedObjectiveTermNLLS &>(currIt->evaluatedObjectiveTerm(m_quantity));

        // We assume the number of residuals/parameters stays constant during
        // optimization.
        const size_t nParams = m_itManager.numParameters();
        Base::set_num_residuals(term.numResiduals());
        // We put all the pattern parameters in a single parameter block.
        Base::mutable_parameter_block_sizes()->assign(1, nParams);
    }

    virtual bool Evaluate(double const * const *parameters,
            double *residuals, double **jacobians) const {
        const size_t nParams = parameter_block_sizes()[0];
        const auto &it = m_itManager.get(nParams, parameters[0]);

        const auto &term = dynamic_cast<const EvaluatedObjectiveTermNLLS &>(it.evaluatedObjectiveTerm(m_quantity));
        const size_t nResiduals = term.numResiduals();

        Real sqrt_weight = sqrt(term.normalizedWeight());
        for (size_t r = 0; r < nResiduals; ++r)
            residuals[r] = sqrt_weight * term.residual(r);

        if (jacobians == NULL) return true;

        for (size_t r = 0; r < nResiduals; ++r) {
            for (size_t p = 0; p < nParams; ++p)
                jacobians[0][r * nParams + p] = sqrt_weight * term.jacobian(r, p);
        }

        return true;
    }

    virtual ~CeresCostWrapper() { }

private:
    const std::string m_quantity;
    IterateManagerBase &m_itManager;
};

class IterationCallback : public ceres::IterationCallback {
public:
    IterationCallback(IterateManagerBase &itManager, ScalarField<Real> &params,
                      const std::string &outPath)
        : m_itManager(itManager), m_params(params), m_outPath(outPath), m_iter(0) {}
    ceres::CallbackReturnType operator()(const ceres::IterationSummary &/*sum*/) {
        const auto &curr = m_itManager.get(m_params.size(), &m_params[0]);
        curr.writeDescription(std::cout);
        std::cout << std::endl;

        if (m_outPath != "")
            curr.writeMeshAndFields(m_outPath + "_" + std::to_string(m_iter));

        ++m_iter;
        return ceres::SOLVER_CONTINUE;
    }

    virtual ~IterationCallback() { }
private:
    IterateManagerBase &m_itManager;
    ScalarField<Real> &m_params;
    std::string m_outPath;
    size_t m_iter;
};

void optimize_ceres_lm(ScalarField<Real> &params,
                       const BoundConstraints &bds,
                       IterateManagerBase &im,
                       const OptimizerConfig &oconfig,
                       const std::string &outPath)
{
    if (!im.isParametric()) throw std::runtime_error("Ceres optimizers only work with parametric inflators");

    const size_t nParams = params.domainSize();
    const IterateBase &it = im.get(nParams, &params[0]);
    ceres::Problem problem;
    for (std::string &name : it.objectiveTermNames()) {
        auto *costFun = new CeresCostWrapper(name, im);
        problem.AddResidualBlock(costFun, NULL, params.data());
    }

    for (size_t p = 0; p < nParams; ++p) {
        if (bds.hasLowerBound.at(p)) problem.SetParameterLowerBound(params.data(), p, bds.lowerBound.at(p));
        if (bds.hasUpperBound.at(p)) problem.SetParameterUpperBound(params.data(), p, bds.upperBound.at(p));
    }

    ceres::Solver::Options options;
    options.update_state_every_iteration = true;
    IterationCallback cb(im, params, outPath);
    options.callbacks.push_back(&cb);
    options.max_num_iterations = oconfig.niters;
    // options.minimizer_type = ceres::LINE_SEARCH;
    // options.line_search_direction_type = ceres::BFGS;
    // options.trust_region_strategy_type = ceres::DOGLEG;
    // options.dogleg_type = ceres::SUBSPACE_DOGLEG;
    // options.use_nonmonotonic_steps = true;
    // options.minimizer_progress_to_stdout = true;
    // options.initial_trust_region_radius = 0.01;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
}

void optimize_ceres_dogleg(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    if (!im.isParametric()) throw std::runtime_error("Ceres optimizers only work with parametric inflators");

    const size_t nParams = params.domainSize();
    const IterateBase &it = im.get(nParams, &params[0]);
    ceres::Problem problem;
    for (std::string &name : it.objectiveTermNames()) {
        auto *costFun = new CeresCostWrapper(name, im);
        problem.AddResidualBlock(costFun, NULL, params.data());
    }

    for (size_t p = 0; p < nParams; ++p) {
        if (bds.hasLowerBound.at(p)) problem.SetParameterLowerBound(params.data(), p, bds.lowerBound.at(p));
        if (bds.hasUpperBound.at(p)) problem.SetParameterUpperBound(params.data(), p, bds.upperBound.at(p));
    }

    ceres::Solver::Options options;
    options.max_num_iterations = oconfig.niters;
    options.update_state_every_iteration = true;
    IterationCallback cb(im, params, outPath);
    options.callbacks.push_back(&cb);

    options.trust_region_strategy_type = ceres::DOGLEG;  // MHS: this seems to be needed in order
    options.dogleg_type = ceres::SUBSPACE_DOGLEG;        // to get similar (rotated) pattern
    options.use_nonmonotonic_steps = true;               // for F [a 0; 0 b] and [b 0; 0 a]

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
}

#else // !HAS_CERES

void optimize_ceres_lm(ScalarField<Real> &params,
                       const BoundConstraints &bds,
                       IterateManagerBase &im,
                       const OptimizerConfig &oconfig,
                       const std::string &outPath)
{
    throw std::runtime_error("Built without Google ceres-solver");
}

void optimize_ceres_dogleg(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    throw std::runtime_error("Built without Google ceres-solver");
}

#endif // HAS_CERES
