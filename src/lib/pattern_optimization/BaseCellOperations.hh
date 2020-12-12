////////////////////////////////////////////////////////////////////////////////
// BaseCellOperations.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Organizes the computations needed to run homogenization and shape
//      optimization on various types of microstructure base cells (triply
//      periodic, orthotropic, ...)
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  08/26/2016 03:08:14
////////////////////////////////////////////////////////////////////////////////
#ifndef BASECELLOPERATIONS_HH
#define BASECELLOPERATIONS_HH

#include "SDConversions.hh"
#include "SDConversionsOrthoCell.hh"
#include <MeshFEM/Future.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/BaseCellType.hh>
#include <MeshFEM/OrthotropicHomogenization.hh>
#include <vector>
#include <memory>
#include <stdexcept>

template<class _Sim>
class BaseCellOperations {
public:
    static constexpr size_t N = _Sim::N;
    using VField  = typename _Sim::VField;
    using ETensor = typename _Sim::ETensor;

    static std::unique_ptr<BaseCellOperations<_Sim>> construct(BaseCellType type, _Sim &sim) {
        // It seems the actual construction must happen in a helper function
        // (_constructBaseCellOps)--couldn't declare this static member function
        // as a friend in derived classes possibly due to _Sim being a dependent
        // type...
        auto bco = _constructBaseCellOps(type, sim);
        bco->m_cellType = type;
        bco->m_solveCellProblems(sim);
        return bco;
    }

    // Solve adjoint problem for the ij^th cell problem
    VField solveAdjointCellProblem(size_t ij, const VField &adjointRHS) const {
        return m_solveProbeSystem(ij, adjointRHS);
    }

    virtual ETensor             homogenizedElasticityTensor()                     const = 0;
    virtual size_t              numReflectedCells()                               const = 0;

    // Cache the homogenizedElasticityTensorDiscreteDifferential, since it may be needed
    // more than once (e.g. for TensorFit and WCS shape derivatives)
    OneForm<ETensor, N> homogenizedElasticityTensorDiscreteDifferential() const {
        if (!m_cachedHETDD)
            m_cachedHETDD = Future::make_unique<OneForm<ETensor, N>>(m_homogenizedElasticityTensorDiscreteDifferential());
        assert(m_cachedHETDD);
        return *m_cachedHETDD;
    }

    // Note: tensors are only computed for elements in the base cell. E.g., for
    // the orthotropic cell case, the tensors for elements in the other
    // reflected cells must be obtained by transforming the tensor for the
    // corresponding orthotropic base cell element.
    // However, for the worst-case stress, the transformations end up cancelling
    // out and the same worst-case stress quantities appear in every reflected
    // element (thus we needn't explicitly transform anything).
    std::vector<ElasticityTensor<Real, _Sim::N, false>>
    macroStrainToMicroStrainTensors() const {
        // The ordinary PeriodicHomogenization::macroStrainToMicroStrainTensors
        // computes the correct tensors (for elements inside the base cell)
        // regardless of base cell type.
        return PeriodicHomogenization::macroStrainToMicroStrainTensors(m_w_ij, m_sim);
    }

    const                _Sim & sim() const { return m_sim; }
    const typename _Sim::Mesh &mesh() const { return m_sim.mesh(); }

    const VField &w_ij(size_t ij)                         const { return m_w_ij.at(ij); }
    const std::vector<VField> &fluctuationDisplacements() const { return m_w_ij; }

    ////////////////////////////////////////////////////////////////////////////
    // Shape Derivative Conversions
    // These are template functions and cannot be made virtual; we must handle
    // dispatch manually...
    ////////////////////////////////////////////////////////////////////////////
    template<size_t _SDDeg>
    ScalarOneForm<N> diff_bdry_from_nsv_functional(const std::vector<Interpolant<Real, N - 1, _SDDeg>> &sd) const {
        if (m_cellType == BaseCellType::TriplyPeriodic) { return SDConversions         ::diff_bdry_from_nsv_functional(sd, mesh()); }
        if (m_cellType == BaseCellType::Orthotropic)    { return SDConversionsOrthoCell::diff_bdry_from_nsv_functional(sd, mesh()); }
        throw std::runtime_error("Invalid cell type.");
    }

    template<typename T = Real>
    OneForm<T, N> diff_bdry_from_diff_vol(const OneForm<T, N> &diff_vol) const {
        if (m_cellType == BaseCellType::TriplyPeriodic) { return SDConversions         ::diff_bdry_from_diff_vol(diff_vol, sim()); }
        if (m_cellType == BaseCellType::Orthotropic)    { return SDConversionsOrthoCell::diff_bdry_from_diff_vol(diff_vol, sim()); }
        throw std::runtime_error("Invalid cell type.");
    }

    VField descent_from_diff_bdry(const ScalarOneForm<N> &diff_bdry) const {
        if (m_cellType == BaseCellType::TriplyPeriodic) { return SDConversions         ::descent_from_diff_bdry(diff_bdry, sim()); }
        if (m_cellType == BaseCellType::Orthotropic)    { return SDConversionsOrthoCell::descent_from_diff_bdry(diff_bdry, sim()); }
        throw std::runtime_error("Invalid cell type.");
    }

    VField descent_from_diff_vol(const ScalarOneForm<N> &diff_vol) const {
        if (m_cellType == BaseCellType::TriplyPeriodic) { return SDConversions         ::descent_from_diff_vol(diff_vol, sim()); }
        if (m_cellType == BaseCellType::Orthotropic)    { return SDConversionsOrthoCell::descent_from_diff_vol(diff_vol, sim()); }
        throw std::runtime_error("Invalid cell type.");
    }

    virtual ~BaseCellOperations() { }

protected:
    BaseCellOperations(const _Sim &sim) : m_sim(sim) { }

    // Needs non-const sim
    virtual void m_solveCellProblems(_Sim &sim) = 0;
    virtual VField m_solveProbeSystem(size_t ij, const VField &rhs) const = 0;
    virtual OneForm<ETensor, N> m_homogenizedElasticityTensorDiscreteDifferential() const = 0;

    // Solved for by derived classes
    std::vector<VField> m_w_ij;
    const _Sim &m_sim;

    mutable std::unique_ptr<OneForm<ETensor, N>> m_cachedHETDD;

    BaseCellType m_cellType;
};

template<class _Sim>
std::unique_ptr<BaseCellOperations<_Sim>> constructBaseCellOps(BaseCellType type, _Sim &sim) { return BaseCellOperations<_Sim>::construct(type, sim); }

// Forward declare _constructBaseCellOps helper function so it can be friended.
template<class _Sim>
std::unique_ptr<BaseCellOperations<_Sim>> _constructBaseCellOps(BaseCellType type, _Sim &sim);

////////////////////////////////////////////////////////////////////////////////
// Triply Periodic Base Cell
////////////////////////////////////////////////////////////////////////////////
// For the traditional triply periodic homogenization, only a single system
// needs to be built, which is just stored in the simulator for historic
// reasons.
template<class _Sim>
class TriplyPeriodicBaseCellOperations : public BaseCellOperations<_Sim> {
public:
    static constexpr size_t N = _Sim::N;
    using Base    = BaseCellOperations<_Sim>;
    using VField  = typename Base::VField;
    using ETensor = typename Base::ETensor;

    virtual ETensor homogenizedElasticityTensor() const override {
        return PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(this->m_w_ij, this->m_sim);
    }


    virtual size_t numReflectedCells() const override { return 1; }

    virtual ~TriplyPeriodicBaseCellOperations() { }

protected:
    TriplyPeriodicBaseCellOperations(const _Sim &sim) : Base(sim) { }

    // Needs non-const sim
    virtual void m_solveCellProblems(_Sim &sim) override {
        PeriodicHomogenization::solveCellProblems(this->m_w_ij, sim, 1e-11);
    }

    // All probes involve the same linear system in the triply periodic case.
    virtual VField m_solveProbeSystem(size_t /* ij */, const VField &rhs) const override {
        return this->m_sim.solve(rhs);
    }

    virtual OneForm<ETensor, N> m_homogenizedElasticityTensorDiscreteDifferential() const override {
        return PeriodicHomogenization::homogenizedElasticityTensorDiscreteDifferential(this->m_w_ij, this->m_sim);
    }

    template<class _S2>
    friend std::unique_ptr<BaseCellOperations<_S2>> _constructBaseCellOps(BaseCellType type, _S2 &sim);
};

////////////////////////////////////////////////////////////////////////////////
// Orthotropic Base Cell
////////////////////////////////////////////////////////////////////////////////
// The different probing strains require different boundary conditions on the
// orthotropic base cell, meaning we need to store several probe systems.
template<class _Sim>
class OrthotropicBaseCellOperations : public BaseCellOperations<_Sim> {
public:
    static constexpr size_t N = _Sim::N;
    using Base    = BaseCellOperations<_Sim>;
    using VField  = typename Base::VField;
    using ETensor = typename Base::ETensor;

    virtual ETensor homogenizedElasticityTensor() const override {
        return PeriodicHomogenization::Orthotropic::homogenizedElasticityTensorDisplacementForm(this->m_w_ij, this->m_sim);
    }

    virtual size_t numReflectedCells() const override { return PeriodicHomogenization::Orthotropic::numReflectedCells(N); }

    virtual ~OrthotropicBaseCellOperations() { }

protected:
    OrthotropicBaseCellOperations(const _Sim &sim) : Base(sim) { }

    // Needs non-const sim
    virtual void m_solveCellProblems(_Sim &sim) override {
        m_probeSystems = PeriodicHomogenization::Orthotropic::solveCellProblems(this->m_w_ij, sim, 1e-11);
    }

    virtual VField m_solveProbeSystem(size_t ij, const VField &rhs) const override {
        return m_systemForProbe(ij).solve(rhs);
    }

    SPSDSystem<Real> &m_systemForProbe(size_t ij) const {
        if (ij >= flatLen(N)) throw std::runtime_error("probe condition index out of bounds");
        if (ij < N) return *m_probeSystems.at(0);
        else        return *m_probeSystems.at((ij - N) + 1);
    }

    virtual OneForm<ETensor, N> m_homogenizedElasticityTensorDiscreteDifferential() const override {
        return PeriodicHomogenization::Orthotropic::homogenizedElasticityTensorDiscreteDifferential(this->m_w_ij, this->m_sim);
    }

    std::vector<std::unique_ptr<SPSDSystem<Real>>> m_probeSystems;

    template<class _S2>
    friend std::unique_ptr<BaseCellOperations<_S2>> _constructBaseCellOps(BaseCellType type, _S2 &sim);
};

////////////////////////////////////////////////////////////////////////////////
// Helper method for allocating (but not initializing!) a base cell.
////////////////////////////////////////////////////////////////////////////////
template<class _Sim>
std::unique_ptr<BaseCellOperations<_Sim>>
_constructBaseCellOps(BaseCellType type, _Sim &sim) {
    std::unique_ptr<BaseCellOperations<_Sim>> bco;

    // Can't use make_unique since it would need to be a friend...
    if       (type == BaseCellType::TriplyPeriodic) bco = std::unique_ptr<TriplyPeriodicBaseCellOperations<_Sim>>(new TriplyPeriodicBaseCellOperations<_Sim>(sim));
    else if  (type == BaseCellType::Orthotropic   ) bco = std::unique_ptr<   OrthotropicBaseCellOperations<_Sim>>(new    OrthotropicBaseCellOperations<_Sim>(sim));
    else throw std::runtime_error("Unknown base cell type.");

    return bco;
}

#endif /* end of include guard: BASECELLOPERATIONS_HH */
