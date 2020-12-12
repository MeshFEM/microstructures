#ifndef CONSTRAINT_TENSORFIT_HH
#define CONSTRAINT_TENSORFIT_HH

#include <optimizers/Constraint.hh>

namespace PatternOptimization {
namespace Constraints {

template<class _Sim>
struct TensorFit : Constraint<_Sim::N> {
    using Base = Constraint<_Sim::N>;
    using Base::m_componentDifferentials;

    static constexpr size_t N = _Sim::N;
    using  OForm = ScalarOneForm<N>;
    using SField = ScalarField<Real>;
    using ETensor = ElasticityTensor<Real, N>;
    using VField = VectorField<Real, N>;

    template<class _Iterate>
    TensorFit(const ETensor &targetS, const _Iterate &it, bool orthotropicSymmetry)
        : Base(ConstraintType::EQUALITY), m_baseCellOps(it.baseCellOps())
    {
        const auto S = it.complianceTensor();
        m_diffS = S - targetS;

        // For reporting only
        auto targetC = targetS.inverse();
        m_CDist = (it.elasticityTensor() - targetC).frobeniusNormSq();
        m_CNormSq = targetC.frobeniusNormSq();

        BENCHMARK_START_TIMER_SECTION("HETDD");
        auto dChVol  = m_baseCellOps.homogenizedElasticityTensorDiscreteDifferential();
        BENCHMARK_STOP_TIMER_SECTION("HETDD");

        BENCHMARK_START_TIMER_SECTION("Convert");
        auto dChBdry = m_baseCellOps.diff_bdry_from_diff_vol(dChVol);
        auto dShBdry = compose([&](const ETensor &e) { ETensor result = S.doubleDoubleContract(e); result *= -1.0; return result; }, dChBdry);

        // If our pattern has orthotropic symmetry, there's no need to add
        // constraints for the zero entries. In fact, this will cause nlopt to
        // fail as the gradient of these constraints is exactly zero.
        // We create a constraint/residual for only the entries of the tensor we
        // care about.
        // TODO: add this to NLLS tensor fit objective
        for (size_t i = 0; i < LinearIndexer<ETensor>::size(); ++i) {
            size_t ij, kl;
            LinearIndexer<ETensor>::linearIndexTo2D(i, ij, kl);
            assert(kl >= ij);

            if (m_ignoreShear && ((ij >= N) || (kl >= N))) continue;
            if (orthotropicSymmetry)
                if (((ij >= N)  || (kl >= N)) && (ij != kl)) continue;
            m_entryForResidual.emplace_back(ij, kl);
        }

        for (size_t i = 0; i < numResiduals(); ++i) {
            m_componentDifferentials.push_back(
                    compose([&](const ETensor &e) { return weightedETensorEntry(e, i); }, dShBdry));
        }

        BENCHMARK_STOP_TIMER_SECTION("Convert");
    }

    // Get the weighted tensor entry corresponding to residual i.
    // Weight each term of each term of (linearly indexed) tensor e so that the
    // squared norm of the resulting vector corresponds to the Frobenius norm of
    // e. (Also, the shear components of e can optionally be zeroed out).
    Real weightedETensorEntry(const ETensor &e, size_t i) const {
        size_t ij, kl;
        std::tie(ij, kl) = m_entryForResidual.at(i);

        Real weight = 1.0;
        if (kl != ij) weight *= sqrt(2); // Account for lower triangle
        if (ignoringShear()) {
            if (ij >= N) weight = 0.0; // Zero out shear components
            if (kl >= N) weight = 0.0; // Zero out shear components
        }
        else {
            if (ij >= N) weight *= sqrt(2); // Left shear doubler
            if (kl >= N) weight *= sqrt(2); // Right shear doubler
        }
        return weight * e.D(ij, kl);
    }

    // The (ij, kl)th constraint residual (kl >= ij). The terms are weighted so
    // that the squared norm of the residual vector corresponds to the Frobenius
    // norm of the rank 4 tensor difference S - S^*.
    virtual SField evaluate() const override {
        SField result(numResiduals());
        for (size_t r = 0; r < numResiduals(); ++r)
            result[r] = weightedETensorEntry(m_diffS, r);

        return result;
    }

    bool ignoringShear() const { return m_ignoreShear; }
    void setIgnoreShear(bool ignore) { m_ignoreShear = ignore; }

    size_t numResiduals() const { return m_entryForResidual.size(); }

    virtual void writeFields(MSHFieldWriter &/* writer */) const override {
        try {
            // TODO: individual component descent output? JS descent?
        }
        catch (const std::exception &e) {
            std::cerr << "Couldn't write TensorFit fields due to exception: " << e.what() << std::endl;
        }
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const override {
        // os << "JS:\t" << 0.5 * evaluate().squaredNorm() << std::endl; // Also dump the analogous NLLS objective value
        // os << "Rel elasticity tensor dist: " << sqrt(m_CDist / m_CNormSq) << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~TensorFit() { }

private:
    // Differentials (one-forms) of each component of the compliance tensor
    const BaseCellOperations<_Sim> &m_baseCellOps;
    bool m_ignoreShear = false;

    ETensor m_diffS;
    Real m_CDist, m_CNormSq;
    // The (ij, kl) flattened elasticity tensor index pairs corresponding to
    // each constraint/residual.
    std::vector<std::pair<size_t, size_t>> m_entryForResidual;
};

// Configuration to be applyed by iterate factory
template<class _Sim>
struct IFConfigTensorFit : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &/* normalizations */) const {
        BENCHMARK_START_TIMER_SECTION("Tensor fit constraint");
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        auto tf = Future::make_unique<TensorFit<_Sim>>(targetS, *it, orthotropicSymmetry);

        tf->setIgnoreShear(ignoreShear);

        it->addConstraint("TensorFit", std::move(tf));
        BENCHMARK_STOP_TIMER_SECTION("Tensor fit constraint");
    }
    ElasticityTensor<Real, N> targetS;
    bool ignoreShear = false;
    bool orthotropicSymmetry = true;
};

}} // end namespace PatternOptimization::Constraints

#endif /* end of include guard: CONSTRAINT_TENSORFIT_HH */
