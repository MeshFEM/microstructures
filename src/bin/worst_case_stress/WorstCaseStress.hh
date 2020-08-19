////////////////////////////////////////////////////////////////////////////////
// WorstCaseStress.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computation and differentiation of microstructure worst-case stress
//      measures (both local and global).
//
//      For now, ***only piecewise constant per-element quantities are used***
//      (strains and stresses must first be averaged over each element).
//
//      WARNING: Cbase is assumed constant throughout the structure (i.e. single
//      base material). Variable materials can be added, but it will require
//      more data storage (or complicate the code).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/03/2015 21:28:18
////////////////////////////////////////////////////////////////////////////////
#ifndef WORSTCASESTRESS_HH
#define WORSTCASESTRESS_HH
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/VonMises.hh>
#include <MeshFEM/Parallelism.hh>

#include <pattern_optimization/SDConversions.hh>
#include <pattern_optimization/BaseCellOperations.hh>

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <tuple>

// Local alias of PeriodicHomogenization namespace.
namespace {
    namespace PH = PeriodicHomogenization;
}

template<size_t N>
using MinorSymmetricRank4Tensor = ElasticityTensor<Real, N, false>;

// Piecewise constant rank 4 tensor fields.
template<size_t N>
using MinorSymmetricRank4TensorField = std::vector<MinorSymmetricRank4Tensor<N>>;

// (Pointwise) Worst-case stress measure
// Constructed by worstCase*Stress() below.
// This supports three types of worst-case stress:
//      1) Frobenius norm
//      2) von Mises
//      3) Max stress
//      4) Stress trace
// (1) and (2) are identical except (2) measures deviatoric micro stress.
// This means that for (2) F maps macro stress to appropriately (scaled)
// deviatoric micro stress, and Cbase is really V : Cbase, where V is the
// linear stress -> "von Mises stress tensor" map. (For 2D, this differs from
// the Deviatoric stress, as described in $CSGFEM/VonMises.hh.
//
// The max stress case (3) is the same as (1), except rank 4 tensor
// eigenvectors are and their eigenvalue are found (instead of eigenstrains).
// The eigenvectors n are converted into the corresponding worst-case macro
// stresses n n^T, and everything else (objective evaluation, shape
// derivatives) is identical.
template<size_t N>
struct WorstCaseStress {
    using SMF = SymmetricMatrixField<Real, N>;

    // d(wc stress)/d(kl^th cell problem strain) on element e
    SymmetricMatrixValue<Real, N> sensitivityToCellStrain(size_t kl, size_t e) const {
        assert(kl < flatLen(N));
        auto result = Cbase.doubleContract(F.at(e).doubleContract(wcMacroStress(e)));
        result *= 2 * (Sh.doubleContract(wcMacroStress(e)))[kl];
        return result;
    }

    // d(wc stress)/d(kl^th cell problem strain) rank 2 tensor field
    void sensitivityToCellStrain(size_t kl, SMF &result) const {
        assert(kl < flatLen(N));
        result.resizeDomain(size());
        for (size_t e = 0; e < size(); ++e) {
            result(e)  = Cbase.doubleContract(F.at(e).doubleContract(wcMacroStress(e)));
            result(e) *= 2 * (Sh.doubleContract(wcMacroStress(e)))[kl];
        }
    }

    // Get microscopic stress measure on element i.
    // Note: this is actually the squared Frobenius/von Mises/max stress
    Real operator()(size_t i) const {
        auto micro = F.at(i).doubleContract(wcMacroStress(i));
        return micro.doubleContract(micro);
    };

    // Get microscopic stress measure field.
    // Note: this is actually the squared Frobenius/von Mises/max stress
    ScalarField<Real> stressMeasure() const {
        ScalarField<Real> result(wcMacroStress.domainSize());
        for (size_t i = 0; i < result.domainSize(); ++i)
            result[i] = (*this)(i);
        return result;
    }

    // Get Frobenius/von Mises/max stress
    ScalarField<Real> sqrtStressMeasure() const {
        ScalarField<Real> result(wcMacroStress.domainSize());
        for (size_t i = 0; i < result.domainSize(); ++i)
            result[i] = sqrt((*this)(i));
        return result;
    }

#if 0 // Disabled since it doesn't compute the actual micro stress in the von Mises case...
    SMF wcMicroStress() const {
        SMF result(wcMacroStress.domainSize());
        for (size_t i = 0; i < F.size(); ++i)
            result(i) = F[i].doubleContract(wcMacroStress(i));
        return result;
    }
#endif

    size_t size() const {
        assert(F.size() == wcMacroStress.domainSize());
        return F.size();
    }

    ////////////////////////////////////////////////////////////////////////////
    // Discrete shape derivatives (Lagrangian). Useful for forward-mode diff.
    ////////////////////////////////////////////////////////////////////////////
    template<class Sim>
    ScalarField<Real> deltaStressMeasure(const Sim &sim,
                const std::vector<typename Sim::VField> &w,
                const typename Sim::VField &delta_p) const {
        size_t numElements = size();
        assert(numElements == sim.mesh().numElements());

        auto deltaSh = PH::deltaHomogenizedComplianceTensor(sim, w, delta_p);
        auto delta_w = PH::deltaFluctuationDisplacements(   sim, w, delta_p);
        auto deltaG  = PH::deltaMacroStrainToMicroStrainTensors(sim, w, delta_w, delta_p);

        MinorSymmetricRank4TensorField<N> deltaF; deltaF.reserve(deltaG.size());
        for (size_t e = 0; e < numElements; ++e) {
            deltaF.push_back(Cbase.doubleContract(deltaG[e].doubleContract(Sh)));
            deltaF.back() += Cbase.doubleContract(G[e].doubleContract(deltaSh));
        }

        ScalarField<Real> delta_s(numElements);
        for (size_t e = 0; e < numElements; ++e) {
            delta_s[e] = 2 * F[e].doubleContract(wcMacroStress(e)).doubleContract(
                    deltaF[e].doubleContract(wcMacroStress(e)));
        }

        return delta_s;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Public data members
    ////////////////////////////////////////////////////////////////////////////
    ElasticityTensor<Real, N> Sh;
    // Minor symmetry:
    //  For the worst-case von Mises analysis, Cbase maps strain to "von Mises
    //  stress." I.e. "Cbase = VonMises : C" where C is the actual printing
    //  base material. In this case, the tensor is *not* major symmetric.
    // TODO: make per-element field:
    //  We store a per-element tensor field rather than a single, homogenous
    //  tensor despite only considering solid/void structures. This is because
    //  the fixed macro load principal stress objective case can be implemented
    //  with a hack of modifying Cbase on each element using the principal
    //  stress direction
    // MinorSymmetricRank4TensorField<N> Cbase;
    MinorSymmetricRank4Tensor<N> Cbase;
    // Per-element macro->micro {stress, deviatoric stress} map.
    // See Worst Case Microstructure writeup for details.
    MinorSymmetricRank4TensorField<N> F, G;
    SMF wcMacroStress;
    std::vector<int>  eigAlgebraicMult;
    std::vector<Real> eigPrincipal, eigSecondary;
};

template<size_t N, bool _majorSymmCBase>
WorstCaseStress<N> worstCaseFrobeniusStress(
        const ElasticityTensor<Real, N, _majorSymmCBase> &Cbase, // Allow non major-symmetric to handle the von Mises case
        const ElasticityTensor<Real, N> &Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain)
{
    WorstCaseStress<N> result;
    result.Cbase = MinorSymmetricRank4Tensor<N>(Cbase);
    result.Sh = Sh;
    result.G = m2mStrain;

    const size_t numElems = m2mStrain.size();

    // macro -> micro stress tensor map
    result.F.reserve(numElems);
    for (size_t i = 0; i < numElems; ++i)
        result.F.emplace_back(Cbase.doubleContract(m2mStrain[i].doubleContract(Sh)));

    // worst-case macro stress
    result.wcMacroStress.resizeDomain(numElems);
    result.eigAlgebraicMult.resize(numElems);
    result.eigPrincipal.resize(numElems);
    result.eigSecondary.resize(numElems);
    for (size_t i = 0; i < numElems; ++i) {
        ElasticityTensor<Real, N> T = result.F[i].transpose().doubleContract(result.F[i]);
        SymmetricMatrixValue<Real, N> sigma;
        std::tie(sigma, result.eigPrincipal[i],
                 result.eigAlgebraicMult[i],
                 result.eigSecondary[i]) = T.maxEigenstrainMultiplicity();
        result.wcMacroStress(i) = sigma;
    }
    return result;
}

template<size_t N>
WorstCaseStress<N> worstCaseVonMisesStress(
        ElasticityTensor<Real, N> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain)
{
    auto V = vonMisesExtractor<N>();
    return worstCaseFrobeniusStress(V.doubleContract(Cbase), Sh, m2mStrain);
}

template<size_t N, bool _majorSymmCBase>
WorstCaseStress<N> fixedLoadFrobeniusStress(
        ElasticityTensor<Real, N, _majorSymmCBase> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain,
        const SymmetricMatrixValue<Real, N> &macroStress)
{
    // Initialize with worst-case macro load
    WorstCaseStress<N> result =
        worstCaseFrobeniusStress(Cbase, Sh, m2mStrain);
    const size_t numElems = result.size();

    // Replace macro stress on each element with the fixed macro load
    for (size_t i = 0; i < numElems; ++i)
        result.wcMacroStress(i) = macroStress;
    return result;
}

template<size_t N, bool _majorSymmCBase>
WorstCaseStress<N> fixedLoadVonMisesStress(
        ElasticityTensor<Real, N, _majorSymmCBase> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain,
        const SymmetricMatrixValue<Real, N> &macroStress)
{
    auto V = vonMisesExtractor<N>();
    return fixedLoadFrobeniusStress(V.doubleContract(Cbase), Sh, m2mStrain,
                                    macroStress);
}

template<size_t N>
WorstCaseStress<N> fixedLoadPrincipalStress(
        ElasticityTensor<Real, N> /* Cbase */,
        ElasticityTensor<Real, N> /* Sh */,
        const MinorSymmetricRank4TensorField<N> &/* m2mStrain */,
        const SymmetricMatrixValue<Real, N> &/* macroStress */)
{
    // Make WorstCaseStress::Cbase a per-element tensor field, replacing it with
    // (n x n x n x n) : Cbase, where n is the element's principal stress
    // direction.
    throw std::runtime_error("unimplemented");
}


template<size_t N>
WorstCaseStress<N> worstCaseMaxStress(
        ElasticityTensor<Real, N> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain)
{
    WorstCaseStress<N> result;
    result.Cbase = Cbase;
    result.Sh = Sh;
    result.G = m2mStrain;

    size_t numElems = m2mStrain.size();

    // TODO: transpose and then correct the shape derivative terms.
    // ALTERNATIVE: transpose all the things (swap Sh, Cbase, transpose G)
    result.F.reserve(numElems);
    for (size_t i = 0; i < numElems; ++i)
        result.F.emplace_back(Cbase.doubleContract(m2mStrain[i].doubleContract(Sh)));

    // worst-case macro stress
    result.wcMacroStress.resizeDomain(numElems);
    MinorSymmetricRank4TensorField<N> T(numElems);
    for (size_t e = 0; e < numElems; ++e) {
        // Actually major symmetric, but conversion is not yet supported
        T[e] = result.F[e].doubleContract(result.F[e].transpose());
    }

    std::ofstream tensorOut("tensors.txt");
    if (!tensorOut.is_open()) throw std::runtime_error("Couldn't open 'tensors.txt' for writing");
    for (size_t e = 0; e < numElems; ++e) {
        // Actually major symmetric, but conversion is not yet supported
        tensorOut << T[e] << std::endl;
    }

    throw std::runtime_error("Rank 4 Tensor Eigenvectors Unimplemented; tensors dumped to 'tensors.txt'");

    return result;
}

template<size_t N>
WorstCaseStress<N> worstCaseStressTrace(
        ElasticityTensor<Real, N> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain)
{
    WorstCaseStress<N> result;
    result.Cbase = Cbase;
    result.Sh = Sh;
    result.G = m2mStrain;

    size_t numElems = m2mStrain.size();

    // Analog to the macro -> micro stress tensor map for the trace objective:
    // we must actually transpose it!
    // TODO: correct the shape derivative terms.
    // ALTERNATIVE: transpose all the things (swap Sh, Cbase, transpose G)
    std::cerr << "WARNING: shape derivative for worstCaseStressTrace must be corrected (apply transposes)" << std::endl;
    result.F.reserve(numElems);
    for (size_t i = 0; i < numElems; ++i) {
        result.F.emplace_back(Cbase.doubleContract(m2mStrain[i].doubleContract(Sh)));
        result.F.back() = result.F.back().transpose();
    }

    result.wcMacroStress.resizeDomain(numElems);
    SymmetricMatrixValue<Real, N> identity;
    identity.clear();
    for (size_t i = 0; i < N; ++i)
        identity(i, i) = 1;

    for (size_t e = 0; e < numElems; ++e) {
        // worst-case "macro stress" for the trace objective is the identity,
        // though it doesn't actually correspond to a macro stress scenario due
        // to the transposition of F (compare to Frobenius/von Mises case).
        result.wcMacroStress(e) = identity;
    }

    return result;
}

// Global worst-case objective of the form
// int j(s, x) dV
// Currently, j (and s) are considered piecewise constant per-element.
template<size_t N, class Integrand>
struct IntegratedWorstCaseObjective {
    using SMF = typename WorstCaseStress<N>::SMF;

    IntegratedWorstCaseObjective() { }
    template<class Mesh> IntegratedWorstCaseObjective(const Mesh &m, const WorstCaseStress<N>  &wcs) { setPointwiseWCS(m, wcs); }
    template<class Mesh> IntegratedWorstCaseObjective(const Mesh &m,       WorstCaseStress<N> &&wcs) { setPointwiseWCS(m, std::move(wcs)); }

    template<class Mesh> void setPointwiseWCS(const Mesh &m, const WorstCaseStress<N>  &wcs) { wcStress = wcs;            integrand.init(m, wcStress); m_updateEvalCache(m); }
    template<class Mesh> void setPointwiseWCS(const Mesh &m,       WorstCaseStress<N> &&wcs) { wcStress = std::move(wcs); integrand.init(m, wcStress); m_updateEvalCache(m); }

    // Evaluate objective by integrating over m
    Real evaluate() const { return m_evalCache; }

    // NOTE: these are *NOT* scaled to account for possible reflected base cell
    // copies (computes just the pointwise value).
    ScalarField<Real> integrandValues() const {
        const size_t numElems = wcStress.size();
        ScalarField<Real> result(numElems);
        for (size_t i = 0; i < numElems; ++i)
            result(i) = integrand.j(wcStress(i), i);
        return result;
    }

    // Sensitivity of global objective integrand to cell problem fluctuation
    // strains (rank 2 tensor field tau^kl in the writeup).
    // NOTE: this is *NOT* scaled to account for possible reflected base cell
    // copies (computes just the pointwise value).
    void tau_kl(size_t kl, SMF &result) const {
        wcStress.sensitivityToCellStrain(kl, result);
        for (size_t e = 0; e < wcStress.size(); ++e)
            result(e) *= integrand.j_prime(wcStress(e), e);
    }

    // Compute objective's partial derivative with respect to the homogenized
    // *elasticity* tensor C^H. This is rank-4 tensor "gamma" from the writeup:
    // dJ/dC^H =
    // 2 * int_omega (j') (G^T : C^base : F : sigma^*) otimes sigma^* dV :: dS^H[v]
    // := -D :: dS^H[v]                 (D = -2 * int_omega (j')...)
    // =   D :: (S^H : dC^H[v] : S^H)
    // = (S^H : D : S^H) :: dC^H[v]     (S^H is major symmetric)
    // NOTE: this is *NOT* scaled to account for possible reflected base cell
    // copies (i.e. the integral is over a single base cell).
    template<class Sim>
    MinorSymmetricRank4Tensor<N> dJ_dCH(const Sim &sim) const {
        MinorSymmetricRank4Tensor<N> D;
        for (size_t col = 0; col < flatLen(N); ++col) {
            auto c = D.DColAsSymMatrix(col);
            for (auto e : sim.mesh().elements()) {
                size_t ei = e.index();
                SymmetricMatrixValue<Real, N> contrib =
                    wcStress.G[ei].transpose().doubleContract(
                        wcStress.Cbase.doubleContract(
                            wcStress.F[ei].doubleContract(
                                wcStress.wcMacroStress(ei))));
                contrib *= e->volume() * wcStress.wcMacroStress(ei)[col] *
                           integrand.j_prime(wcStress(ei), ei);
                c += contrib;
            }
        }
        D *= -2.0;
        return wcStress.Sh.doubleDoubleContract(D);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Discrete shape derivative of the integrated worst-case stress objective
    // under mesh vertex perturbation delta_p (forward mode).
    // J is an integral of a piecewise-constant per-element field, so it's easy
    // to differentiate:
    //      J = int_omega j(s) dV = sum_e j(s)_e vol(e)
    //     dJ = sum_e [j(s)_e dvol(e) + j'(s)_e ds_e vol(e)]
    ////////////////////////////////////////////////////////////////////////////
    template<class Sim>
    Real deltaJ(const Sim &sim,
                const std::vector<typename Sim::VField> &w,
                const typename Sim::VField &delta_p) const {
        // Homogenized elasticity tensor change term: gamma :: dCh[v]
        auto deltaCh = PH::deltaHomogenizedElasticityTensor(sim, w, delta_p);
        Real dJ = dJ_dCH(sim).quadrupleContract(deltaCh);
        const auto &mesh = sim.mesh();

        // int tau_kl : delta strain(w^kl) dV
        auto delta_w = PH::deltaFluctuationDisplacements(sim, w, delta_p);
        {
            SMF tau;
            for (size_t kl = 0; kl < delta_w.size(); ++kl) {
                // Sum is only over the "upper triangle" of integrands
                Real shearDoubler = (kl >= N) ? 2.0 : 1.0;
                tau_kl(kl, tau);
                auto delta_we = sim.deltaAverageStrainField(w[kl], delta_w[kl], delta_p);
                for (auto e : mesh.elements())
                    dJ += tau(e.index()).doubleContract(delta_we(e.index())) * e->volume() * shearDoubler;
            }
        }

        // Volume dilation term: int j div v dV
        std::vector<VectorND<N>> cornerPerturbations;
        for (auto e : mesh.elements()) {
            size_t ei = e.index();
            sim.extractElementCornerValues(e, delta_p, cornerPerturbations);
            dJ += integrand.j(wcStress(ei), ei) * e->relativeDeltaVolume(cornerPerturbations) * e->volume();
        }

        return dJ * numReflectedCopies;
    }

    // Compute the linear functional (one-form) dJ[v], represented as a
    // per-vertex vector field to be dotted with a **periodic** mesh vertex
    // velocity.
    //      dJ[v] = sum_i <delta_j[i], v[i]>
    //      (delta_j[i][c] = partial_derivative(J, v[i][c]))
    // where the sum is over all mesh vertices. (Vertices on the periodic
    // boundary get "half" contributions)
    //
    // -delta_j can be interpreted as a steepest descent direction in the
    // vertex position space R^(N*|v|), but it is not a good descent direction
    // for non-uniform meshes. For a better, mesh-independent direction, one
    // should use a different metric (E.g. ||v||_M = int <v(x), v(x)> dx = v[i]
    // M_ij v[j], where M is the deg 1 mass matrix. Or for Newton's method, the
    // Hessian.).
    template<class Sim>
    ScalarOneForm<Sim::N>
    adjointDeltaJ(const BaseCellOperations<Sim> &baseCellOps) const {
        // Dilation and delta strain terms
        const auto mesh = baseCellOps.mesh();
        size_t nv = mesh.numVertices();

        BENCHMARK_START_TIMER_SECTION("Adjoint Cell Problem");
        BENCHMARK_START_TIMER("Compute Tau_kl");
        // Cache tau_pq and adjoint solutions
        std::vector<SMF> tau(flatLen(N));
        for (size_t pq = 0; pq < flatLen(N); ++pq)
            tau_kl(pq, tau[pq]);
        BENCHMARK_STOP_TIMER("Compute Tau_kl");

        std::vector<VectorField<Real, N>> lambda;
        lambda.reserve(flatLen(N));
        for (size_t pq = 0; pq < flatLen(N); ++pq) {
            lambda.emplace_back(baseCellOps.solveAdjointCellProblem(pq,
                    baseCellOps.sim().perElementStressFieldLoad(tau[pq])));
        }

        const std::vector<VectorField<Real, N>> &w = baseCellOps.fluctuationDisplacements();

        BENCHMARK_STOP_TIMER_SECTION("Adjoint Cell Problem");

        BENCHMARK_START_TIMER_SECTION("Dilation and Delta strain");

        using OF = ScalarOneForm<N>;
        OF delta_j(nv);
        delta_j.clear();

#if USE_TBB
        tbb::combinable<OF> sum(delta_j);
#endif

        auto accumElementContrib = [&](size_t ei) {
            auto e = mesh.element(ei);
#if USE_TBB
            OF &result = sum.local();
#else
            OF &result = delta_j;
#endif
            // j contribution to dilation integrand
            // Part 1: int_w [j - strain(p^kl):stress^kl] grad(lambda_m) dx
            // Part 1a: int_w j grad(lambda_m) dx
            Real dilationIntegrand = integrand.j(wcStress(e.index()), e.index());
            const auto C = e->E();
            for (size_t pq = 0; pq < flatLen(N); ++pq) {
                // Sum is only over the "upper triangle" of integrands
                Real shearDoubler = (pq >= Sim::N) ? 2.0 : 1.0;
                using Strain = typename Sim::Strain;
                Strain strain_lambda, strain_u;
                e->strain(e, lambda[pq], strain_lambda);
                e->strain(e,      w[pq], strain_u);
                strain_u += Sim::SMatrix::CanonicalBasis(pq);

                using Stress = typename Sim::Stress;
                Stress stress_lambda, stress_u;
                for (size_t i = 0; i < Stress::numNodalValues; ++i) {
                    stress_lambda[i] = C.doubleContract(strain_lambda[i]);
                    stress_u     [i] = C.doubleContract(strain_u     [i]);
                }

                // pq^th contribution to dilation integrand.
                // Part 1b: int_w strain(p^kl):stress^kl * grad(lambda_m) dx
                dilationIntegrand -= Quadrature<N, 2 * Strain::Deg>::integrate(
                        [&] (const EvalPt<N> &p) {
                    return strain_lambda(p).doubleContract(stress_u(p));
                }) * shearDoubler;

                // Delta strain term
                // Part 2: int_w [grad(lambda_m) . (stress * p^kl + (strain(p^kl):C - tau) w^kl)] grad(phi_n) dx
                for (auto n : e.nodes()) {
                    auto gradPhi_n = e->gradPhi(n.localIndex());

                    Interpolant<VectorND<N>, N, Stress::Deg> glam_functional;
                    glam_functional = tau[pq](e.index()).contract(w[pq](n.index()));
                    for (size_t i = 0; i < glam_functional.size(); ++i) {
                        glam_functional[i] -= stress_lambda[i].contract(w[pq](n.index()))
                                            + stress_u[i].contract(lambda[pq](n.index()));
                    }

                    for (auto v_m : e.vertices()) {
                        auto gradLam_m = e->gradBarycentric().col(v_m.localIndex());
                        Interpolant<Real, N, Strain::Deg> mu_grad_lam_m;
                        for (size_t i = 0; i < mu_grad_lam_m.size(); ++i)
                            mu_grad_lam_m[i] = gradLam_m.dot(glam_functional[i]);
                        result(v_m.index()) -= Quadrature<N, 2 * Strain::Deg>::
                            integrate([&](const EvalPt<N> &pt) { return
                                    (mu_grad_lam_m(pt) * gradPhi_n(pt)).eval();
                                }, e->volume() * shearDoubler);
                    }
                }
            }
            for (auto v : e.vertices()) {
                result(v.index()) += dilationIntegrand * e->volume()
                    * e->gradBarycentric().col(v.localIndex());
            }
        };

#if USE_TBB
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, mesh.numElements()),
            [&](const tbb::blocked_range<size_t> &r) {
                for (size_t ei = r.begin(); ei < r.end(); ++ei) accumElementContrib(ei);
            });
        delta_j = sum.combine([](const OF &a, const OF &b) { return a + b; } );
#else
        for (auto e : mesh.elements()) { accumElementContrib(e.index()); }
#endif

        BENCHMARK_STOP_TIMER_SECTION("Dilation and Delta strain");

        // Gamma term
        BENCHMARK_START_TIMER_SECTION("Gamma Term");
        {
            MinorSymmetricRank4Tensor<N> gamma = dJ_dCH(baseCellOps.sim());
            auto dCh = baseCellOps.homogenizedElasticityTensorDiscreteDifferential();
            OF gtermNew = compose([&](const ElasticityTensor<Real, N> &e) { return e.quadrupleContract(gamma); }, dCh);
            delta_j += gtermNew;

            // using ETensorSD = PH::BEHTensorGradInterpolant<Sim>;
            // using GTermSD   = Interpolant<Real, ETensorSD::K, ETensorSD::Deg>;
            // std::vector<ETensorSD> sdCh = PH::homogenizedElasticityTensorGradient(w, sim);
            // std::vector<GTermSD> gammaTermFunctional(sdCh.size());
            // for (size_t i = 0; i < sdCh.size(); ++i) {
            //     for (size_t j = 0; j < ETensorSD::size(); ++j)
            //         gammaTermFunctional[i][j] = gamma.quadrupleContract(sdCh[i][j]);
            // }
            // auto gterm_bdry = SDConversions::diff_bdry_from_nsv_functional(gammaTermFunctional, mesh);
            // ScalarOneForm<N> gtermOld(mesh.numVertices());
            // gtermOld.clear();
            // for (auto bv : mesh.boundaryVertices())
            //     gtermOld(bv.volumeVertex().index()) += gterm_bdry(bv.index());
            // MSHFieldWriter writer("debug_gamma.msh", mesh);
            // writer.addField("gamma term new", gtermNew.asVectorField());
            // writer.addField("gamma term old", gtermOld.asVectorField());
        }
        BENCHMARK_STOP_TIMER_SECTION("Gamma Term");

        if (numReflectedCopies != 1) delta_j *= Real(numReflectedCopies);
        return delta_j;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Public data members
    ////////////////////////////////////////////////////////////////////////////
    WorstCaseStress<N> wcStress;
    Integrand integrand;
    // number of reflected copies (e.g. for orthocell)
    // NOTE: assumes the integrand is invariant to the reflection operation!!!
    // (Should be true for WCS)
    size_t numReflectedCopies = 1;

private:
    Real m_evalCache = 0.0;
    template<class Mesh> void m_updateEvalCache(const Mesh &m) {
        assert(m.numElements() == wcStress.size());
        m_evalCache = 0;
        for (auto e : m.elements())
            m_evalCache += integrand.j(wcStress(e.index()), e.index()) * e->volume();

        m_evalCache *= Real(numReflectedCopies);
    }
};

// int_omega worst_case_stress dV
// i.e. j(s, x) = s -> j' = 1.
struct WCStressIntegrandTotal {
    template<size_t N, class Mesh>
    void init(const Mesh &/* m */, const WorstCaseStress<N> &/* wcs */) { }

    // Derivative of global objective integrand wrt worst case stress.
    static Real j(Real wcStress, size_t /* x_i */) { return wcStress; }
    static Real j_prime(Real /* wcStress */, size_t /* x_i */) { return 1; }
};

// int_omega worst_case_stress^p dV
// i.e. j(s, x) = s^p -> j' = p s^(p - 1).
// if there is a target stress:
//   j(s) = (s-target)^p, if s >= target
//   j(s) = 0 , otherwise
//   If p >= 2, guaranteed to be smooth
struct WCStressIntegrandLp {
    template<size_t N, class Mesh>
    void init(const Mesh &/* m */, const WorstCaseStress<N> &/* wcs */) { }

    WCStressIntegrandLp() { }

    Real j(Real wcStress, size_t /* x_i */) const {
        if (wcStress >= target)
            return pow(wcStress - target, p);
        else
            return 0;
    }
    // Derivative of global objective integrand wrt worst case stress.
    Real j_prime(Real wcStress, size_t /* x_i */) const {
        if (wcStress >= target)
            return p * pow(wcStress - target, p - 1);
        else
            return 0;
    }

    Real p = 1.0;
    Real target = 0.0;
};

// Global max worst-case objective:
// max_{x in omega} s(x)
// We do approximate gradient computation by treating this as an integrated
// objective int_omega j(s, x) where
// j(s, x) = s * w(x)
// j'(s, x) = w(x)
// w(x) = 1.0 / V if worst case stress at s is at the max.
// V = total volume of all regions at max stress.
struct WCStressIntegrandLinf {
    template<size_t N, class Mesh>
    void init(const Mesh &m, const WorstCaseStress<N> &wcs) {
        assert(wcs.size() == m.numElements());
        assert(wcs.size() > 0);

        std::vector<size_t> maxStressElements(1, 0);
        Real maxStress = wcs(0);
        for (auto e : m.elements()) {
            Real val = wcs(e.index());
            if (maxStress < val) {
                maxStressElements.assign(1, e.index());
                maxStress = val;
            }
            else if (maxStress == val) {
                maxStressElements.push_back(e.index());
            }
        }

        Real volumeAtMaxStress = 0;
        for (size_t e : maxStressElements)
            volumeAtMaxStress += m.element(e)->volume();

        weightFunction.assign(m.numElements(), 0.0);
        for (size_t e : maxStressElements)
            weightFunction.at(e) = m.element(e)->volume() / volumeAtMaxStress;
    }

    Real j(Real s, size_t x_i) const { return s * weightFunction.at(x_i); }
    Real j_prime(Real /* s */, size_t x_i) const { return weightFunction.at(x_i); }

    std::vector<Real> weightFunction;
};

// Take the pth root of a WCS objective function.
// Based on chain rule, all derivatives are weighted by
// 1/p f^((1 / p) - 1), where f = SubObjective::eval();
template<class SubObjective>
struct PthRootObjective : public SubObjective {
using Base = SubObjective;
    // Forward all constructor args to SubObjective
    template<typename... Args>
    PthRootObjective(Args&&... args)
        : SubObjective(std::forward<Args>(args)...) { }

    Real evaluate() const {
        if (p == 1) return SubObjective::evaluate();
        return pow(SubObjective::evaluate(), 1.0 / p);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Discrete shape derivative under mesh vertex perturbation delta_p
    // (forward mode).
    ////////////////////////////////////////////////////////////////////////////
    template<class Sim>
    Real deltaJ(const Sim &sim,
                const std::vector<typename Sim::VField> &w,
                const typename Sim::VField &delta_p) const {
        Real pder = SubObjective::deltaJ(sim, w, delta_p);
        if (p == 1) return pder;
        return m_gradientScale() * pder;
    }

    template<class Sim>
    ScalarOneForm<Sim::N>
    adjointDeltaJ(const BaseCellOperations<Sim> &baseCellOps) const {
        ScalarOneForm<Sim::N> pder = SubObjective::adjointDeltaJ(baseCellOps);
        if (p == 1) return pder;
        pder *= m_gradientScale();
        return pder;
    }

    Real p = 1.0;
private:
    Real m_gradientScale() const {
        Real value = Base::evaluate();

        // Protect for the case where value is 0.0, which would make computation undefined.
        // When value is 0.0, derivative will also be 0.
        if (value > 0.0)
            return (1.0 / p) * pow(value, 1.0 / p - 1.0);
        else
            return 0.0;
    }
};


#endif /* end of include guard: WORSTCASESTRESS_HH */
