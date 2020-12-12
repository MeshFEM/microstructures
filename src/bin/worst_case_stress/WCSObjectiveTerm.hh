#ifndef WCSOBJECTIVETERM_HH
#define WCSOBJECTIVETERM_HH

#include "WorstCaseStress.hh"
#include <pattern_optimization/ObjectiveTerm.hh>
#include <pattern_optimization/ObjectiveTermNormalizations.hh>
#include <pattern_optimization/IterateFactory.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <stdexcept>

namespace PatternOptimization {
namespace ObjectiveTerms {

template<class _Sim, class _WCSObjectiveType>
struct WorstCaseStress : ObjectiveTerm<_Sim::N> {
    using    Base = ObjectiveTerm<_Sim::N>;
    using   OForm = ScalarOneForm<_Sim::N>;
    using  VField = typename _Sim::VField;
    using SMatrix = typename _Sim::SMatrix;
    template<class _Iterate>
    WorstCaseStress(const _Iterate &it,
                    Real globalObjectivePNorm, Real globalObjectiveRoot, Real stressTarget,
                    const std::string &measure,
                    const SMatrix *macroLoad)
        : m_baseCellOps(it.baseCellOps())
    {
        // Configure objective
        m_wcs_objective.integrand.p        = globalObjectivePNorm;
        m_wcs_objective.integrand.target   = stressTarget;
        m_wcs_objective.p                  = globalObjectiveRoot;
        m_wcs_objective.numReflectedCopies = m_baseCellOps.numReflectedCells();

        const auto &mesh = it.simulator().mesh();

        // Worst case stress currently assumes that the base material is
        // constant, so we can read it off a single element.
        auto CBase = mesh.element(0)->E();
        if (measure == "frobenius") {
            if (!macroLoad) {
                m_wcs_objective.setPointwiseWCS(mesh,
                    worstCaseFrobeniusStress(CBase, it.complianceTensor(),
                        m_baseCellOps.macroStrainToMicroStrainTensors()));
            }
            else {
                m_wcs_objective.setPointwiseWCS(mesh,
                    fixedLoadFrobeniusStress(CBase, it.complianceTensor(),
                        m_baseCellOps.macroStrainToMicroStrainTensors(), *macroLoad));
            }
        }
        else if (measure == "vonmises") {
            if (!macroLoad) {
                m_wcs_objective.setPointwiseWCS(mesh,
                    worstCaseVonMisesStress(CBase, it.complianceTensor(),
                        m_baseCellOps.macroStrainToMicroStrainTensors()));
            }
            else {
                m_wcs_objective.setPointwiseWCS(mesh,
                    fixedLoadVonMisesStress(CBase, it.complianceTensor(),
                        m_baseCellOps.macroStrainToMicroStrainTensors(), *macroLoad));
            }
        }
        else if (measure == "maxnorm") {
            if (!macroLoad) {
                m_wcs_objective.setPointwiseWCS(mesh,
                    worstCaseMaxStress(CBase, it.complianceTensor(),
                        m_baseCellOps.macroStrainToMicroStrainTensors()));
            }
            else {
                m_wcs_objective.setPointwiseWCS(mesh,
                    fixedLoadPrincipalStress(CBase, it.complianceTensor(),
                        m_baseCellOps.macroStrainToMicroStrainTensors(), *macroLoad));
            }
        }
        else throw std::runtime_error("Unknown worst-case measure: " + measure);

        // Compute and store WCS's boundary differential one-form
        m_diff_vol = m_wcs_objective.adjointDeltaJ(m_baseCellOps);
        this->m_differential = m_baseCellOps.diff_bdry_from_diff_vol(m_diff_vol);
    }

    virtual Real evaluate() const override { return m_wcs_objective.evaluate(); }

    virtual void writeFields(MSHFieldWriter &writer) const override {
        ScalarField<Real> j = m_wcs_objective.integrandValues();

        // plot max stress field
        ScalarField<Real> sqrtStressField = m_wcs_objective.wcStress.sqrtStressMeasure();
        writer.addField("Pointwise WCS", sqrtStressField);

        // plot difference between max stress and target at each element
        ScalarField<Real> diff(sqrtStressField.domainSize());
        for (size_t i = 0; i < diff.domainSize(); ++i) {
            if (sqrtStressField[i] > sqrt(m_wcs_objective.integrand.target))
                diff[i] = sqrtStressField[i] - sqrt(m_wcs_objective.integrand.target);
            else
                diff[i] = 0.0;
        }
        writer.addField("WCS above target", diff);

        // writer.addField("j", j);

        ScalarField<Real> eigPrincipal(m_wcs_objective.wcStress.eigPrincipal),
                          eigSecondary(m_wcs_objective.wcStress.eigSecondary),
                          eigMult(m_wcs_objective.wcStress.size()),
                          dist(m_wcs_objective.wcStress.size());
        for (size_t i = 0; i < eigMult.domainSize(); ++i) {
            eigMult[i] = Real(m_wcs_objective.wcStress.eigAlgebraicMult.at(i));
            dist[i] = (eigPrincipal[i] - eigSecondary[i]) / eigPrincipal[i];
        }
        // writer.addField("Principal eigenvalue", eigPrincipal, DomainType::PER_ELEMENT);
        // writer.addField("Secondary eigenvalue", eigSecondary, DomainType::PER_ELEMENT);
        // writer.addField("Eigenvalue multiplicity", eigMult, DomainType::PER_ELEMENT);
        writer.addField("Eigenvalue relative distance", dist, DomainType::PER_ELEMENT);

#if 0
        try {
            writer.addField("WCS Steepest Descent VVel",
                            m_baseCellOps.descent_from_diff_vol(m_diff_vol),
                            DomainType::PER_NODE);
        }
        catch (const std::exception &e) {
            std::cerr << "Couldn't write volume descent velocity: " << e.what() << std::endl;
        }
#endif

        auto bdryVel = m_baseCellOps.descent_from_diff_bdry(this->m_differential);
        VField xferBdryVel(m_baseCellOps.mesh().numVertices());
        xferBdryVel.clear();
        for (auto v : m_baseCellOps.mesh().vertices()) {
            auto bv = v.boundaryVertex();
            if (!bv) continue;
            xferBdryVel(v.index()) = bdryVel(bv.index());
        }
        writer.addField("WCS Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const override {
        os << "Max Ptwise WCS:\t" << sqrt(m_wcs_objective.wcStress.stressMeasure().maxMag()) << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~WorstCaseStress() { }
private:
    OForm m_diff_vol; // per-volume-vertex differential
    const BaseCellOperations<_Sim> &m_baseCellOps;
    _WCSObjectiveType m_wcs_objective;
};

// Configuration to be applied by iterate factory
template<class _Sim, class _WCSObjectiveType = PthRootObjective<IntegratedWorstCaseObjective<_Sim::N, WCStressIntegrandLp>>>
struct IFConfigWorstCaseStress : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        BENCHMARK_START_TIMER_SECTION("WCS Term");
        auto wcs = Future::make_unique<WorstCaseStress<_Sim, _WCSObjectiveType>>(
                *it, globalObjectivePNorm, globalObjectiveRoot, target, measure, macroLoad.get());
        wcs->setWeight(weight);

        // WCS normalization is the initial worst case stress
        if (!normalizations.isSet("WCS")) {
            if (target > 0.0)
                // Protection in cases where a target is used, since wcs could be then equal to 0.0
                normalizations.set("WCS", 1.0);
            else
                normalizations.set("WCS", 1.0 / wcs->evaluate());
        }

        wcs->setNormalization(normalizations["WCS"]);
        it->addObjectiveTerm("WCS", std::move(wcs));
        BENCHMARK_STOP_TIMER_SECTION("WCS Term");
    }

    Real weight = 1.0;
    Real globalObjectivePNorm = 1.0;
    Real globalObjectiveRoot = 1.0;
    Real target = 0.0;
    std::string measure = std::string("frobenius");
    std::unique_ptr<typename _Sim::SMatrix> macroLoad;
};

}}

#endif /* end of include guard: WCSOBJECTIVETERM_HH */
