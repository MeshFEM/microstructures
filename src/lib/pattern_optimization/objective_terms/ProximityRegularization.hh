#ifndef OBJECTIVETERMPROXIMITYREGULARIZATION_HH
#define OBJECTIVETERMPROXIMITYREGULARIZATION_HH

#include "../ObjectiveTerm.hh"

namespace PatternOptimization {
namespace ObjectiveTerms {

// alpha * (p - p0).(p - p0) / 2
template<size_t N>
struct ProximityRegularization : public NLLSObjectiveTerm<N> {
    using SField = ScalarField<Real>;
    using VField = VectorField<Real, N>;

    ProximityRegularization(const std::vector<Real> &p,
                            const std::vector<Real> &targetParams)
        : m_currParams(p), m_targetParams(targetParams) { }

    virtual SField gradp(const std::vector<VField> &/*bdrySVels*/) const {
        SField result(m_currParams);
        result -= m_targetParams;
        result *= this->m_weight;
        return result;
    }

    virtual Real evaluate() const { return (this->m_weight / 2.0) * (m_targetParams.values() - m_currParams.values()).squaredNorm(); }

    virtual SField residual() const {
        size_t nParams = m_targetParams.domainSize();
        SField result(nParams);
        Real sqrtWeight = std::sqrt(this->m_weight);
        for (size_t i = 0; i < nParams; ++i)
            result[i] = sqrtWeight * (m_currParams[i] - m_targetParams[i]);
        return result;
    }

    virtual Eigen::MatrixXd jacobian(const std::vector<VField> &/*bdrySVels*/) const {
        size_t nParams = m_targetParams.domainSize();
        Eigen::MatrixXd result(nParams, nParams);
        result.setIdentity();
        result *= std::sqrt(this->m_weight);
        return result;
    }

    virtual ~ProximityRegularization() { }

private:
    SField m_currParams, m_targetParams;
};

// Configuration to be applied by iterate factory
struct IFConfigProximityRegularization : public IFConfig {
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        if (!normalizations.isSet("ProximityRegularization"))
            normalizations.set("ProximityRegularization", 1.0); // TODO? Base on parameter range?

        auto pr = Future::make_unique<ProximityRegularization<_Iterate::_N>>(it->params(), targetParams);
        pr->setWeight(weight);

        pr->setNormalization(normalizations["ProximityRegularization"]);
        it->addObjectiveTerm("ProximityRegularization", std::move(pr));
    }
    Real weight = 0.0;
    std::vector<Real> targetParams;
};


}}

#endif /* end of include guard: OBJECTIVETERMPROXIMITYREGULARIZATION_HH */
