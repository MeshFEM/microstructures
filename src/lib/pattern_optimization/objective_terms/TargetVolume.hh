#ifndef OBJECTIVETERMJVOL_HH
#define OBJECTIVETERMJVOL_HH

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"
#include <MeshFEM/OneForm.hh>
#include <pattern_optimization/SDConversions.hh>

namespace PatternOptimization {
    namespace ObjectiveTerms {

        // weight * 1/2 * (V_cur - V_target)^2
        template<class _Sim>
        struct TargetVolumeTerm : public NLLSObjectiveTerm<_Sim::N> {
            using SField = ScalarField<Real>;
            using VField = VectorField<Real, _Sim::N>;
            using  OForm = ScalarOneForm<_Sim::N>;

            static Real compute_cost(Real targetVol, Real currentVol) {
                return 0.5 * (currentVol - targetVol) * (currentVol - targetVol);
            }

            static OForm compute_differential(Real targetVol, Real currentVol, _Sim &simulator, OForm &diff_vol) {
                OForm result = (currentVol - targetVol) * SDConversions::diff_bdry_from_diff_vol(diff_vol, simulator); // TODO: consider orthotropic cells?
                return result;
            }

            template<class _Iterate>
            TargetVolumeTerm(Real tvol, _Iterate &it) : m_targetVol(tvol) {
                const auto &sim = it.simulator();
                m_vol = sim.mesh().volume();

                OForm m_diff_vol = sim.deltaVolumeForm();
                this->m_differential = compute_differential(m_targetVol, m_vol, it.simulator(), m_diff_vol);
            }

            // Evaluate objective function (without weight). In this case, it means: (V_cur - V_target)^2
            virtual Real evaluate() const override {
                return compute_cost(m_targetVol, m_vol);
            }

            // Compute single residual for volume
            virtual SField residual() const override {
                SField result(1);
                result(0) = (m_vol - m_targetVol);

                return result;
            }

            // Derivative of residual wrt parameter p.
            // Compute the derivative of inner term. In our case, we only have one function, which is the volume
            virtual Eigen::MatrixXd jacobian(const std::vector<VField> &bdrySVels) const override {
                const size_t np = bdrySVels.size();
                Eigen::MatrixXd result(1, np); // only one residual term: V_current - V_target

                for (size_t p = 0; p < np; ++p) {
                    result(0, p) = this->m_differential[bdrySVels[p]];
                }

                return result;
            }

            void setTarget(Real v) { m_targetVol = v; }
            virtual ~TargetVolumeTerm() { }
        private:
            Real m_targetVol;
            Real m_vol;
        };

        // Configuration to be applyed by iterate factory
        template<class _Sim>
        struct IFConfigTargetVolume : public IFConfig {
            template<class _Iterate>
            void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
                auto tv = Future::make_unique<TargetVolumeTerm<_Sim>>(targetVolume, *it);
                tv->setWeight(weight);

                if (!normalizations.isSet("JVol"))
                    normalizations.set("JVol", 1.0); // TODO? What about target vol = 0?

                tv->setNormalization(normalizations["JVol"]);
                it->addObjectiveTerm("JVol", std::move(tv));
            }
            Real weight = 0.0;
            Real targetVolume = 0.0;
        };

    }}

#endif /* end of include guard: OBJECTIVETERMJVOL_HH */
