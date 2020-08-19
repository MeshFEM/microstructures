////////////////////////////////////////////////////////////////////////////////
// IsotropicFitRel.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Minimizes the homogenized tensor's distance to isotropy:
//          min_shape min_C^iso ||C^H - C^iso||^2_F / ||C^H||^2_F
//  Due to the envelope theorem, the change in C^iso can be neglected when
//  computing derivatives. Thus we can formulate the problem as a relative
//  fit to the projection of C^H onto the isotropic space and use essentially
//  the same code as for TensorFit.
//
//  This relative version prevents shrinkage.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/04/2017 00:34:00
////////////////////////////////////////////////////////////////////////////////
#ifndef ISOTROPICFITREL_HH
#define ISOTROPICFITREL_HH


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
struct IsotropicFitRel : ObjectiveTerm<_Sim::N> {
    using Base = ObjectiveTerm<_Sim::N>;

    static constexpr size_t N = _Sim::N;
    using  OForm = ScalarOneForm<N>;
    using SField = ScalarField<Real>;
    using ETensor = ElasticityTensor<Real, N>;
    using VField = VectorField<Real, N>;
    using LI = LinearIndexer<ETensor>;

    template<class _Iterate>
    IsotropicFitRel(const _Iterate &it) : m_baseCellOps(it.baseCellOps()) {
        const auto Ch = it.elasticityTensor();
        ETensor targetC = closestIsotropicTensor(Ch);

        m_diffC = Ch - targetC;

        m_CDistSq = m_diffC.frobeniusNormSq();
        m_CNormSq = Ch.frobeniusNormSq();

        BENCHMARK_START_TIMER_SECTION("HETDD");
        auto dChVol  = m_baseCellOps.homogenizedElasticityTensorDiscreteDifferential();
        BENCHMARK_STOP_TIMER_SECTION("HETDD");

        BENCHMARK_START_TIMER_SECTION("Convert");
        auto dChBdry = m_baseCellOps.diff_bdry_from_diff_vol(dChVol);
        this->m_differential = compose([&](const ETensor &e) {
                return 2 * (m_diffC.quadrupleContract(e) / m_CNormSq -
                                (m_CDistSq / (m_CNormSq * m_CNormSq)) * Ch.quadrupleContract(e)
                           ); },
            dChBdry);

        BENCHMARK_STOP_TIMER_SECTION("Convert");
    }

    virtual Real evaluate() const override { return m_CDistSq / m_CNormSq; }

    virtual void writeDescription(std::ostream &os, const std::string &name) const override {
        os << "Rel elasticity tensor dist: " << sqrt(m_CDistSq / m_CNormSq) << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~IsotropicFitRel() { }

private:
    // Differentials (one-forms) of each component of the compliance tensor
    const BaseCellOperations<_Sim> &m_baseCellOps;
    bool m_ignoreShear = false;
    ETensor m_diffC;
    Real m_CDistSq, m_CNormSq;
};

// Configuration to be applied by iterate factory
template<class _Sim>
struct IFConfigIsotropyFitRel : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        BENCHMARK_START_TIMER_SECTION("Relative Isotropic Fit");
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        auto ifit = Future::make_unique<IsotropicFitRel<_Sim>>(*it);
        ifit->setWeight(weight);

        if (!normalizations.isSet("JIsoRel")) {
            normalizations.set("JIsoRel", 1.0);
        }

        ifit->setNormalization(normalizations["JIsoRel"]);
        it->addObjectiveTerm("JIsoRel", std::move(ifit));
        BENCHMARK_STOP_TIMER_SECTION("Relative Isotropic Fit");
    }
    Real weight = 1.0;
};

}}

#endif /* end of include guard: ISOTROPICFITREL_HH */
