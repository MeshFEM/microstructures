////////////////////////////////////////////////////////////////////////////////
// IsotropicFit.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Minimizes the homogenized tensor's distance to isotropy:
//          min_shape min_C^iso ||C^H - C^iso||^2_F
//  Due to the envelope theorem, the change in C^iso can be neglected when
//  computing derivatives. Thus we can formulate the problem as a NLLS
//  fit to the projection of C^H onto the isotropic space and use essentially
//  the same code as for TensorFit.
//
//  Note: this version can cause the tensor to shrink if convergence isn't fast
//  enough, especially with slsqp.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/04/2017 00:34:00
////////////////////////////////////////////////////////////////////////////////
#ifndef ISOTROPICFIT_HH
#define ISOTROPICFIT_HH

#include <MeshFEM/TensorProjection.hh>

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"
#include "../SDConversions.hh"

#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <MeshFEM/ElasticityTensor.hh>
#include <MeshFEM/LinearIndexer.hh>

#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/MSHFieldWriter.hh>

namespace PatternOptimization {
namespace ObjectiveTerms {

template<class _Sim>
struct IsotropicFit : NLLSObjectiveTerm<_Sim::N> {
    using Base = NLLSObjectiveTerm<_Sim::N>;

    static constexpr size_t N = _Sim::N;
    using  OForm = ScalarOneForm<N>;
    using SField = ScalarField<Real>;
    using ETensor = ElasticityTensor<Real, N>;
    using VField = VectorField<Real, N>;
    using LI = LinearIndexer<ETensor>;

    template<class _Iterate>
    IsotropicFit(const _Iterate &it, bool useFixedTarget = false, const ETensor &fixedTargetC = ETensor())
			: m_baseCellOps(it.baseCellOps()) {
        const auto Ch = it.elasticityTensor();
        ETensor targetC = closestIsotropicTensor(Ch);
        if (useFixedTarget) targetC = fixedTargetC;

        m_diffC = Ch - targetC;

        m_CDistSq = m_diffC.frobeniusNormSq();
        m_CNormSq = Ch.frobeniusNormSq();

        BENCHMARK_START_TIMER_SECTION("HETDD");
        auto dChVol  = m_baseCellOps.homogenizedElasticityTensorDiscreteDifferential();
        BENCHMARK_STOP_TIMER_SECTION("HETDD");

        BENCHMARK_START_TIMER_SECTION("Convert");
        auto dChBdry = m_baseCellOps.diff_bdry_from_diff_vol(dChVol);
        this->m_differential = compose([&](const ETensor &e) { return m_diffC.quadrupleContract(e); }, dChBdry);

        for (size_t i = 0; i < LI::size(); ++i) {
            m_component_differentials.push_back(
                    compose([&](const ETensor &e) { return weightedETensorEntry(e, i); }, dChBdry));
        }
        BENCHMARK_STOP_TIMER_SECTION("Convert");
    }

    virtual Real evaluate() const override { return 0.5 * m_diffC.frobeniusNormSq(); }

    static constexpr size_t numResiduals() { return LI::size(); }

	// The (ij, kl)th residual (kl >= ij) for the nonlinear least squares (a
    // single term of the Frobenius distance). The terms are weighted so
    // that the squared norm of the residual vector corresponds to the
    // Frobenius norm of the rank 4 tensor difference S - S^*.
    virtual SField residual() const override {
        SField result(numResiduals());
        for (size_t r = 0; r < numResiduals(); ++r)
            result(r) = weightedETensorEntry(m_diffC, r);

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
        if (ij >=  N) weight *= sqrt(2); // Left shear doubler
        if (kl >=  N) weight *= sqrt(2); // Right shear doubler
        return weight * e.D(ij, kl);
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const override {
        os << "Rel elasticity tensor dist: " << sqrt(m_CDistSq / m_CNormSq) << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~IsotropicFit() { }

private:
    // Differentials (one-forms) of each component of the compliance tensor
    const BaseCellOperations<_Sim> &m_baseCellOps;
    std::vector<OForm> m_component_differentials;
    bool m_ignoreShear = false;
    ETensor m_diffC;
    Real m_CDistSq, m_CNormSq;
};

// Configuration to be applied by iterate factory
template<class _Sim>
struct IFConfigIsotropyFit : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) {
        BENCHMARK_START_TIMER_SECTION("Isotropic Fit");
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");

        if (useFixedTarget && !fixedTargetSet) {
            fixedTargetC = closestIsotropicTensor(it->elasticityTensor());
            std::cout << "JIso fitting to fixed tensor:" << std::endl;
            std::cout << fixedTargetC << std::endl;
            std::cout << std::endl;
            fixedTargetSet = true;
        }

        std::unique_ptr<IsotropicFit<_Sim>> ifit;
        if (useFixedTarget) ifit = Future::make_unique<IsotropicFit<_Sim>>(*it, useFixedTarget, fixedTargetC);
        else                ifit = Future::make_unique<IsotropicFit<_Sim>>(*it);

        ifit->setWeight(weight);

        // Default JIso normalization is (twice) the initial tensor's squared Frobenius norm
        if (!normalizations.isSet("JIso")) {
            Real scale = 2 / it->elasticityTensor().frobeniusNormSq();
            normalizations.set("JIso", scale);
        }


        ifit->setNormalization(normalizations["JIso"]);
        it->addObjectiveTerm("JIso", std::move(ifit));
        BENCHMARK_STOP_TIMER_SECTION("Isotropic Fit");
    }
    Real weight = 1.0;
    // Fit to the closest tensor to the original.
    bool useFixedTarget = false;
    bool fixedTargetSet = false;
    ElasticityTensor<Real, N> fixedTargetC;
};

}}

#endif /* end of include guard: ISOTROPICFIT_HH */
