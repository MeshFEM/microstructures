#ifndef OBJECTIVETERMJS_HH
#define OBJECTIVETERMJS_HH

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"
#include "../SDConversions.hh"

#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <MeshFEM/ElasticityTensor.hh>
#include <MeshFEM/LinearIndexer.hh>

#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/MSHFieldWriter.hh>

#include <stdexcept>
#include <iostream>

namespace PatternOptimization {
namespace ObjectiveTerms {

// NLLS: 1/2 ||S^H - S^*||_F^2
// The individual residual components are entries in the upper triangle of
// flattened tensor (S^H - S^*), weighted so that the squared residual vector
// norm computes the full rank 4 tensor Frobenius norm.
//
// For more accurate estimation, we estimate the residual and evaluate the
// objective in terms of it.
template<class _Sim>
struct TensorFit : NLLSObjectiveTerm<_Sim::N> {
    using Base = NLLSObjectiveTerm<_Sim::N>;

    static constexpr size_t N = _Sim::N;
    using  OForm = ScalarOneForm<N>;
    using SField = ScalarField<Real>;
    using ETensor = ElasticityTensor<Real, N>;
    using VField = VectorField<Real, N>;
    using LI = LinearIndexer<ETensor>;

    template<class _Iterate>
    TensorFit(const ETensor &targetS, const _Iterate &it) : m_baseCellOps(it.baseCellOps()) {
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
        this->m_differential = compose([&](const ETensor &e) { return m_diffS.quadrupleContract(e); }, dShBdry);

        double diff_norm = sqrt(this->m_differential.asVectorField().frobeniusNormSq());
        std::cout << "||Tensor Fit Shape Derivative||: " << diff_norm  << std::endl;

#if 0
        {
            auto dJVol = compose([&](const ETensor &e) { return -m_diffS.quadrupleContract(S.doubleDoubleContract(e)); }, dChVol);
            MSHFieldWriter volWriter("volfields.msh", m_baseCellOps.mesh());
            volWriter.addField("dJVol", dJVol.asVectorField(), DomainType::PER_NODE);

            MSHBoundaryFieldWriter bdryWriter("bdryfields.msh", m_baseCellOps.mesh());
            const auto &sim = it.simulator();
            SField isInternal(sim.mesh().numBoundaryElements());
            isInternal.clear();
            for (auto be : sim.mesh().boundaryElements())
                isInternal(be.index()) = be->isInternal ? 1.0 : 0.0;
            bdryWriter.addField("isInternal", isInternal, DomainType::PER_ELEMENT);
        }
#endif

        for (size_t i = 0; i < LI::size(); ++i) {
            m_component_differentials.push_back(
                    compose([&](const ETensor &e) { return weightedETensorEntry(e, i); }, dShBdry));
        }
        BENCHMARK_STOP_TIMER_SECTION("Convert");
    }

    virtual Real evaluate() const override { return 0.5 * m_diffS.frobeniusNormSq(); }

    bool ignoringShear() const { return m_ignoreShear; }
    void setIgnoreShear(bool ignore) { m_ignoreShear = ignore; }

    static constexpr size_t numResiduals() { return LI::size(); }

	// The (ij, kl)th residual (kl >= ij) for the nonlinear least squares (a
    // single term of the Frobenius distance). The terms are weighted so
    // that the squared norm of the residual vector corresponds to the
    // Frobenius norm of the rank 4 tensor difference S - S^*.
    virtual SField residual() const override {
        SField result(numResiduals());
        for (size_t r = 0; r < numResiduals(); ++r)
            result(r) = weightedETensorEntry(m_diffS, r);

        return result;
    }

    // Derivative of residual(ij, kl) wrt parameter p:
    // d/dp (S_ijkl - target_ijkl) = d/dp S_ijkl = <gradS_ijkl, vn_p>
    // The terms are weighted in accordance with the residual weighting above.
    virtual Eigen::MatrixXd jacobian(const std::vector<VField> &bdrySVels) const override {
        const size_t np = bdrySVels.size();
        Eigen::MatrixXd result(numResiduals(), np);
        for (size_t p = 0; p < np; ++p) {
            for (size_t r = 0; r < numResiduals(); ++r)
                result(r, p) = m_component_differentials.at(r)[bdrySVels[p]];
        }
        return result;
    }

    // Weight each term of each term of (linearly indexed) tensor e so that the
    // squared norm of the resulting vector corresponds to the Frobenius norm of
    // e. (Also, the shear components of e can optionally be zeroed out).
    Real weightedETensorEntry(const ETensor &e, size_t i) const {
        size_t ij, kl;
        LI::linearIndexTo2D(i, ij, kl);
        assert(kl >= ij);

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

    virtual void writeFields(MSHFieldWriter &writer) const override {
        try {
            // VField differential(m_baseCellOps.mesh().numVertices());
            // differential.clear();
            // for (auto bv : m_baseCellOps.mesh().boundaryVertices()) {
            //     differential(bv.volumeVertex().index()) =
            //         this->m_differential.asVectorField()(bv.index());
            // }
            // writer.addField("JS Differential", differential, DomainType::PER_NODE);

            auto bdryVel = m_baseCellOps.descent_from_diff_bdry(this->m_differential);
            VField xferBdryVel(m_baseCellOps.mesh().numVertices());
            xferBdryVel.clear();
            for (auto v : m_baseCellOps.mesh().vertices()) {
                auto bv = v.boundaryVertex();
                if (!bv) continue;
                xferBdryVel(v.index()) = bdryVel(bv.index());
            }
            writer.addField("JS Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
        }
        catch (const std::exception &e) {
            std::cerr << "Couldn't write TensorFit fields due to exception: " << e.what() << std::endl;
        }
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const override {
        os << "Rel elasticity tensor dist: " << sqrt(m_CDist / m_CNormSq) << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~TensorFit() { }

private:
    // Differentials (one-forms) of each component of the compliance tensor
    const BaseCellOperations<_Sim> &m_baseCellOps;
    std::vector<OForm> m_component_differentials;
    bool m_ignoreShear = false;
    ETensor m_diffS;
    Real m_CDist, m_CNormSq;
};

// Configuration to be applied by iterate factory
template<class _Sim>
struct IFConfigTensorFit : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        BENCHMARK_START_TIMER_SECTION("Tensor fit");
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        auto tf = Future::make_unique<TensorFit<_Sim>>(targetS, *it);
        tf->setWeight(weight);

        // Default JS normalization is (twice) target tensor's squared Frobenius norm
        // This means the normalized term is the relative squared frobenius distance.
        if (!normalizations.isSet("JS"))
            normalizations.set("JS", 2.0 / targetS.frobeniusNormSq());

        tf->setNormalization(normalizations["JS"]);
        tf->setIgnoreShear(ignoreShear);
        it->addObjectiveTerm("JS", std::move(tf));
        BENCHMARK_STOP_TIMER_SECTION("Tensor fit");
    }
    Real weight = 1.0;
    ElasticityTensor<Real, N> targetS;
    bool ignoreShear = false;
};

}}

#endif /* end of include guard: OBJECTIVETERMJS_HH */
