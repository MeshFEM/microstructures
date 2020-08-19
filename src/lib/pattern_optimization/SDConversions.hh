////////////////////////////////////////////////////////////////////////////////
// SDConversions.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Conversions between various forms of shape derivatives.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  05/02/2016 18:18:27
////////////////////////////////////////////////////////////////////////////////
#ifndef SDCONVERSIONS_HH
#define SDCONVERSIONS_HH

#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <MeshFEM/Functions.hh>
#include <MeshFEM/GaussQuadrature.hh>
#include <MeshFEM/MassMatrix.hh>
#include <isosurface_inflator/ShapeVelocityInterpolator.hh>
#include <vector>

namespace SDConversions {

// Compute a boundary differential from a NSV functional. This is
// the per-boundary-vertex vector field whose inner product with a (periodic)
// vector field v_b representing a piecewise linear shape velocity computes the
// same functional as sd[dot(v_b, n)].
// (Nonperiodic, vertices on the periodic cell get "half" (partial) contributions.)
template<class _FEMMesh, size_t _SDDeg>
ScalarOneForm<_FEMMesh::K> diff_bdry_from_nsv_functional(
        const std::vector<Interpolant<Real, _FEMMesh::K - 1, _SDDeg>> &sd,
        const _FEMMesh &mesh)
{
    constexpr size_t BdryK = _FEMMesh::K - 1;
    constexpr size_t N     = _FEMMesh::K;
    ScalarOneForm<N> result(mesh.numBoundaryVertices());
    result.clear();

    assert(sd.size() == mesh.numBoundaryElements());
    for (auto be : mesh.boundaryElements()) {
        if (be->isInternal) continue;
        const auto &nsv = sd.at(be.index());
        // Integrate against each boundary vertex's linear shape function
        for (auto bv : be.vertices()) {
            Interpolant<Real, BdryK, 1> bary_bv;
            bary_bv = 0.0;
            bary_bv[bv.localIndex()] = 1.0;
            result(bv.index()) += be->normal() *
                Quadrature<BdryK, _SDDeg + 1>::integrate(
                        [&](const EvalPt<BdryK> &p) { return
                            bary_bv(p) * nsv(p);
                        }, be->volume());
        }
    }

    return result;
}

// Compute a boundary differential one-form corresponding from a volume
// differential one-form "diff_vol". Evaluating this boundary form on a boundary
// velocity field is equivalent to evaluating the volume form on a smoothly
// (harmonically) extended volume velocity field:
//      <diff_bdry, v_bdry> := <diff_vol, v_interp>
//                           = <diff_vol, I * v_bdry>
//                           = <I^T diff_vol, v_bdry>
//  for all v_bdry ==> diff_bdry = I^T diff_vol
//  Where I is the linear bdry -> vol velocity interpolation operator.
template<class _Sim, typename T = Real>
OneForm<T, _Sim::N> diff_bdry_from_diff_vol(
        const OneForm<T, _Sim::N> &diff_vol,
        const _Sim &sim)
{
    ShapeVelocityInterpolator interpolator(sim);
    return interpolator.adjoint(sim, diff_vol);
}

// Compute from a boundary differential one-form, "diff_bdry", a
// mesh-independent steepest-descent boundary velocity field ("steepest" with
// respect to the boundary geometry's L^2 norm):
//      min dJ[g]                   (dJ[g] = sum_i diff_bdry[i] . g[i])
//       g          ==> diff_bdry + 2 l M g = 0
//   g^T M g = h^2
//      ==> g \propto -M^-1 diff_bdry
// Here M is the boundary mass matrix (g^T M g = int_bdry g(x) . g(x) dA(x),
// where g(x) is the interpolated scalar field corresponding to vector g).
// This g can also be interpreted as the "Riesz representative" for dJ.
// Periodicity details:
//    We seek a periodic g, and the boundary integral defining M is over
//    the "true" boundary (i.e. the boundary remaining after stitching
//    together identified vertices). We simplify the problem by reducing to
//    "periodic dof" variables.
//    Introducing matrix S that selects the dof value for each vertex:
//          g = S g_dof
//    (Note that S^T sums values over all identified vertices and places
//     them on the corresponding periodic dof.)
//
//    diff_bdry is non-periodic; it contains partial contributions on
//    each vertex (so it correctly evaluates the one-form when dotted with
//    a *periodic* discrete boundary nodal vector field). Thus the summed
//    quantity S^T diff_bdry gives the one-form acting on periodic
//    boundary node displacements.
//
//    Letting M_dof = S^T M S, we arrive at the equation:
//         g_dof \propto -  M_dof^-1 S^T diff_bdry
//    ==>  g     \propto -S M_dof^-1 S^T diff_bdry
// We let the caller decide on a normalization by simply computing:
// @return      -S M_dof^-1 S^T diff_bdry
template<class _Sim>
typename _Sim::VField descent_from_diff_bdry(
        const ScalarOneForm<_Sim::N> &diff_bdry,
        const _Sim &sim)
{
    // Extract boundary vertex periodicity from m_sim.
    // Unfortunately, sim's periodic boundary conditions are on volume
    // nodes, not boundary vertices, so we have to reindex to get
    // contiguous variables
    const auto mesh = sim.mesh();
    size_t nbv = mesh.numBoundaryVertices();
    std::vector<size_t> dofForBdryVertex; dofForBdryVertex.reserve(nbv);

    constexpr size_t NO_VAR = std::numeric_limits<size_t>::max();
    std::vector<size_t> dofForSimNodalDoF(sim.numDoFs(), NO_VAR);
    size_t numDoFs = 0;
    for (auto bv : mesh.boundaryVertices()) {
        size_t simDoF = sim.DoF(bv.volumeVertex().node().index());
        size_t &dof = dofForSimNodalDoF.at(simDoF);
        if (dof == NO_VAR)
            dof = numDoFs++;
        dofForBdryVertex.push_back(dof);
    }

    // Determine which variables correspond to true boundary vertices
    // (if any non-periodic boundary element touches them), and extract an
    // array for testing periodic boundary element membership.
    std::vector<bool> isTrueBoundaryVertex(mesh.numBoundaryVertices(), false),
                      isPeriodicBE(mesh.numBoundaryElements(), false);
    for (auto be : mesh.boundaryElements()) {
        isPeriodicBE[be.index()] = be->isInternal;
        if (be->isInternal) continue;
        for (size_t i = 0; i < be.numVertices(); ++i)
            isTrueBoundaryVertex.at(be.vertex(i).index()) = true;
    }

    // Determine dofs corresponding to false boundary vertices (to be constrained)
    std::vector<bool> dofVisited(numDoFs, false); // Add each dof only once
    std::vector<size_t> falseBoundaryDoFs;
    for (auto bv : mesh.boundaryVertices()) {
        if (isTrueBoundaryVertex.at(bv.index())) continue;

        size_t dof = dofForBdryVertex.at(bv.index());
        if (dofVisited.at(dof)) continue;
        falseBoundaryDoFs.push_back(dof);
        dofVisited[dof] = true;
    }

    // Build M_dof, ignoring contributions from periodic boundary elements
    auto M = MassMatrix::construct<1>(mesh.boundary(), false, isPeriodicBE);
    // ==> S M S^T
    M.reindexVariables(numDoFs, dofForBdryVertex);
    SPSDSystem<Real> M_dof(M);

    // Enforce zero velocity on false boundary vertices
    // (Since we skipped periodic bdry elements' contributions, the system
    //  would be singular otherwise).
    M_dof.fixVariables(falseBoundaryDoFs,
                       std::vector<Real>(falseBoundaryDoFs.size(), 0.0));

    typename _Sim::VField g(nbv);
    for (size_t c = 0; c < _Sim::N; ++c) {
        // Apply -S^T
        std::vector<Real> neg_S_t_diff_bdry(numDoFs, 0.0);
        for (auto bv : mesh.boundaryVertices()) {
            neg_S_t_diff_bdry.at(dofForBdryVertex.at(bv.index()))
                -= diff_bdry(bv.index())[c];
        }

        auto g_dof = M_dof.solve(neg_S_t_diff_bdry);

        // Apply S
        for (auto bv : mesh.boundaryVertices())
            g(bv.index())[c] = g_dof.at(dofForBdryVertex.at(bv.index()));
    }

    return g;
}

// Compute from a volume differential one-form, "diff_vol", a mesh-independent
// steepest-descent volume velocity field ("steepest" with respect to the period
// cell geometry's L^2 norm).
// This amounts to computing a periodic velocity field, g, whose integrated
// inner product over the mesh against velocity v gives the change in
// objective (i.e. a "Riesz representative" for dJ):
//      int_omega g(x) . v(x) dx := dJ[v] = - sum_i diff_vol[i] . v[i]
// where v(x) := v[k] * phi^k(x),   (periodic)
//       g(x) := g[k] * phi^k(x),   (periodic)
// and phi_k is a scalar vertex shape function.
//
// Correctly handling the periodicity requires care:
// Notice the sum over i in the one-form evaluation is over **all** vertices, so
// it's effectively summing all partial contributions in diff_vol[i] for
// periodically identified vertices.
//
// Introducing "selection" matrix S that distributes the reduced periodic DoF
// values to the corresponding vertices (v = S v_dof):
//   dJ[v] = sum_i diff_vol[i] . (S v_dof)[i] = sum_i (S^T diff_vol)[i] . v_dof[i]
// (here, S^T sums the partial contributions each identified vertex into a
//  per-dof total).
// Likewise, we can write v(x) = (S v_dof)[k] phi^k(x),
//                        g(x) = (S g_dof)[k] phi^k(x) ==>
//      int_omega g(x) . v(x) dx := -dJ[v] = -sum_i diff_vol[i] . v[i]
//    = g_dof[k] (int_omega S_sk phi^s(x) phi^t(x) S_tl) . v_dof[l]
//   := g_dof[k] . v_dof[l] [M_dof]_kl
// with reduced mass matrix M_dof = S^T M S.
// Finally, for this inner product to equal -dJ[v] for all v, we see:
//      sum_i (M_dof g_dof)[i] . v_dof[i] = -sum_i (S^T diff_vol)[i] . v_dof[i]
//  ==> g_dof = -M_dof^{-1} S^T diff_vol
// and we can recover the velocity field coefficients g = S g_dof
template<class _Sim>
typename _Sim::VField descent_from_diff_vol(
        const ScalarOneForm<_Sim::N> &diff_vol,
        const _Sim &sim)
{
    // Determine S, the map from reduced periodic DoFs to vertices
    // Unfortunately, sim's periodic boundary conditions are on nodes, not
    // just vertices, so we have to reindex to get contiguous variables
    const auto &mesh = sim.mesh();
    size_t nv = mesh.numVertices();
    constexpr size_t NO_VAR = std::numeric_limits<size_t>::max();
    std::vector<size_t> dofForSimNodalDoF(sim.numDoFs(), NO_VAR);
    std::vector<size_t> dofForVertex;
    dofForVertex.reserve(nv);
    size_t numDoFs = 0;
    for (auto v : mesh.vertices()) {
        size_t &dof = dofForSimNodalDoF.at(sim.DoF(v.node().index()));
        if (dof == NO_VAR) dof = numDoFs++;
        dofForVertex.push_back(dof);
    }

    auto M = MassMatrix::construct<1>(mesh);
    // ==> S M S^T
    M.reindexVariables(numDoFs, dofForVertex);
    SPSDSystem<Real> M_dof(M);

    typename _Sim::VField g(nv);
    for (size_t c = 0; c < _Sim::N; ++c) {
        // Apply -S^T
        std::vector<Real> neg_S_t_diff_vol(numDoFs, 0.0);
        for (auto v : mesh.vertices())
            neg_S_t_diff_vol.at(dofForVertex[v.index()]) -= diff_vol(v.index())[c];

        auto g_dof = M_dof.solve(neg_S_t_diff_vol);

        // Apply S
        for (auto v : mesh.vertices())
            g(v.index())[c] = g_dof.at(dofForVertex.at(v.index()));
    }

    return g;
}

}

#endif /* end of include guard: SDCONVERSIONS_HH */
