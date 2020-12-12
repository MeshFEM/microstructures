#ifndef PERIODICSMOOTHINGREGULARIZATION_H
#define PERIODICSMOOTHINGREGULARIZATION_H

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"
#include <MeshFEM/OneForm.hh>
#include <pattern_optimization/SDConversions.hh>
#include <Eigen/Sparse>

namespace PatternOptimization {
    namespace ObjectiveTerms {

        template<class _Sim>
        struct PeriodicSmoothingRegularizationTerm : public ObjectiveTerm<_Sim::N> {
            using SField = ScalarField<Real>;
            using VField = VectorField<Real, _Sim::N>;
            using  OForm = ScalarOneForm<_Sim::N>;

            static constexpr size_t N = _Sim::N;

            static Eigen::SparseMatrix<Real> computeLaplacianMatrix(_Sim &simulator) {
                std::vector<std::set<size_t>> adj(simulator.mesh().numVertices());

                std::unique_ptr<PeriodicCondition<N>> pc = Future::make_unique<PeriodicCondition<N>>(simulator.mesh());
                std::vector<std::vector<size_t>> correspondingVertices(simulator.mesh().numVertices());

                for (auto v : simulator.mesh().boundaryVertices()) {
                    size_t vNodeIndex = v.node().volumeNode().index();
                    std::vector<size_t> vNodeIndices = pc->identifiedNodes(vNodeIndex);
                    std::vector<size_t> vIndices;

                    for (size_t i = 0; i < vNodeIndices.size(); i++) {
                        size_t vIndex = simulator.mesh().node(vNodeIndices[i]).vertex().index();
                        vIndices.push_back(vIndex);
                    }

                    // Sort corresponding vertices. This way, always first in the list has lowest index
                    std::sort(vIndices.begin(), vIndices.end());
                    correspondingVertices[v.node().volumeNode().vertex().index()] = vIndices;
                }

                for (auto be : simulator.mesh().boundaryElements()) {
                    if (be->isInternal) {
                        continue;
                    }

                    for (auto v1 : be.vertices()) {
                        size_t v1Index = v1.node().volumeNode().vertex().index();

                        for (auto v2 : be.vertices()) {
                            size_t v2Index = v2.node().volumeNode().vertex().index();

                            if (v2 == v1)
                                continue;
                            else if (!adj[v1Index].count(v2Index)) {
                                adj[v1Index].insert(v2Index);
                            }
                        }
                    }
                }

                int numVertices = simulator.mesh().numVertices();
                Eigen::SparseMatrix<Real> L(numVertices, numVertices);
                std::vector<Eigen::Triplet<Real>> triplets;
                for (auto ve : simulator.mesh().boundaryVertices()) {
                    size_t veIndex = ve.node().volumeNode().vertex().index();

                    if (adj[veIndex].size() == 0) {
                        continue;
                    }

                    size_t number_neighbors = 0;
                    for (auto idx : correspondingVertices[veIndex]) {
                        for (auto neighbor : adj[idx]) {
                            number_neighbors++;
                        }
                    }

                    // For each corresponding vertex, add a row.
                    for (auto idx : correspondingVertices[veIndex]) {
                        triplets.push_back(Eigen::Triplet<Real>(veIndex, idx, -1.0 * adj[idx].size() / number_neighbors));

                        for (auto neighbor : adj[idx]) {
                            triplets.push_back(Eigen::Triplet<Real>(veIndex, neighbor, 1.0 / number_neighbors));
                        }
                    }
                }

                L.setFromTriplets(triplets.begin(), triplets.end());

                return L;
            }

            static Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> computePositionMatrix(_Sim &simulator, int cols = 2) {
                int numVertices = simulator.mesh().numVertices();
                Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> p(numVertices, cols);
                p.fill(0.0);
                for (auto v : simulator.mesh().vertices()) {
                    p.row(v.index()) = v.node()->p;
                }

                return p;
            }

            static Real computeCost(Eigen::SparseMatrix<Real> L, Eigen::MatrixXd p) {
                Eigen::MatrixXd resultingMatrix = L * p;

                Real result = 0.0;
                for (unsigned i=0; i<resultingMatrix.rows(); i++) {
                    result += resultingMatrix.row(i).squaredNorm();
                }

                return result;
            }

            static SField computeSmoothingField(Eigen::SparseMatrix<Real> L, Eigen::MatrixXd p) {
                Eigen::MatrixXd resultingMatrix = L * p;
                SField result(p.rows());

                for (unsigned i=0; i<resultingMatrix.rows(); i++) {
                    result[i] = resultingMatrix.row(i).squaredNorm();
                }

                return result;
            }

            static OForm computeVolumeDifferential(Eigen::SparseMatrix<Real> L, Eigen::MatrixXd p) {
                Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> diff_L = 2 * L.transpose() * L * p;
                Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> transpose = diff_L.transpose();

                assert(diff_L.cols() == 2 || diff_L.cols() == 3);
                assert(diff_L.rows() == L.cols());

                transpose.resize(transpose.cols() * transpose.rows(), 1);
                VField vectorField(transpose);
                OForm oneform(vectorField);

                return oneform;
            }

            static OForm computeDifferential(_Sim &simulator, Eigen::SparseMatrix<Real> L, Eigen::MatrixXd p) {
                OForm volumeDifferential = computeVolumeDifferential(L, p);

                return SDConversions::diff_bdry_from_diff_vol(volumeDifferential, simulator);
            }

            template<class _Iterate>
            PeriodicSmoothingRegularizationTerm(_Iterate &it) : m_sim(it.simulator()) {
                m_L = computeLaplacianMatrix(it.simulator());
                m_p = computePositionMatrix(it.simulator(), _Sim::N);
                m_smoothingField = computeSmoothingField(m_L, m_p);
                this->m_differential = computeDifferential(it.simulator(), m_L, m_p);
            }

            // Evaluate objective function (without weight)
            virtual Real evaluate() const override {
                Real result = computeCost(m_L, m_p);
                return result;
            }

            virtual void writeFields(MSHFieldWriter &writer) const override {
                writer.addField("Pointwise smoothing", m_smoothingField);

                auto bdryVel = SDConversions::descent_from_diff_bdry(this->m_differential, m_sim);
                VField xferBdryVel(m_L.cols());
                xferBdryVel.clear();
                for (auto v : m_sim.mesh().vertices()) {
                    auto bv = v.boundaryVertex();
                    if (!bv) continue;
                    xferBdryVel(v.index()) = bdryVel(bv.index());
                }
                writer.addField("Smoothing Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
            }

            virtual ~PeriodicSmoothingRegularizationTerm() { }
        private:
            _Sim &m_sim;
            Eigen::SparseMatrix<Real> m_L;
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> m_p;
            SField m_smoothingField;
        };

        // Configuration to be applied by iterate factory
        template<class _Sim>
        struct IFConfigPeriodicSmoothingRegularization: public IFConfig {
            template<class _Iterate>
            void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
                auto sr = Future::make_unique<PeriodicSmoothingRegularizationTerm<_Sim>>(*it);
                sr->setWeight(weight);

                if (!normalizations.isSet("PeriodicSmoothingRegularization"))
                    normalizations.set("PeriodicSmoothingRegularization", 1.0);

                sr->setNormalization(normalizations["PeriodicSmoothingRegularization"]);
                it->addObjectiveTerm("PeriodicSmoothingRegularization", std::move(sr));
            }
            Real weight = 0.0;
        };

    }
}

#endif //PERIODICSMOOTHINGREGULARIZATION_H