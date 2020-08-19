////////////////////////////////////////////////////////////////////////////////
// SDConversionsOrthoCell.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Conversions between various forms of shape derivatives.
//      Just like SDConversions.hh, but for orthotropic base cells.
//
//      In the future we could make this a class specialized by base cell type.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  08/26/2016 16:07:27
////////////////////////////////////////////////////////////////////////////////
#ifndef SDCONVERSIONSORTHOCELL_HH
#define SDCONVERSIONSORTHOCELL_HH

#include "SDConversions.hh"
#include "ShapeVelocityInterpolatorOrthoCell.hh"
#include <MeshFEM/OrthotropicHomogenization.hh>
#include <MeshFEM/MassMatrix.hh>
#include <stdexcept>

namespace SDConversionsOrthoCell {

// Conversion of nsv functionals is identical to the periodic case (though the
// discussion there about half/partial contributions and periodic velicity
// fields is not relevant).
template<class _FEMMesh, size_t _SDDeg>
ScalarOneForm<_FEMMesh::K> diff_bdry_from_nsv_functional(
        const std::vector<Interpolant<Real, _FEMMesh::K - 1, _SDDeg>> &sd,
        const _FEMMesh &mesh)
{
    return SDConversions::diff_bdry_from_nsv_functional(sd, mesh);
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
    ShapeVelocityInterpolatorOrthoCell<_Sim::N> interpolator(sim);
    return interpolator.adjoint(sim, diff_vol);
}

// Compute the steepest descent boundary velocity field from a one-form on the
// boundary:
//      min      diff_bdry[vb]
//  ||vb||_M = 1
// The mass-matrix M is used as the metric.
//  ==> vb = normalize(- M^-1 diff_bdry)
// The mass matrix integrated over the orthotropic cell is exactly
// 1/numReflectedCells that which would be computed by reflecting and using
// periodic DoFs.
template<class _Sim>
typename _Sim::VField descent_from_diff_bdry(
        const ScalarOneForm<_Sim::N> &diff_bdry,
        const _Sim &sim)
{
    constexpr size_t N = _Sim::N;
    const auto &mesh = sim.mesh();
    size_t NR = PeriodicHomogenization::Orthotropic::numReflectedCells(N);
    if (diff_bdry.domainSize() != mesh.numBoundaryVertices())
        throw std::runtime_error("Invalid diff_bdry size");

    // Shouldn't actually need to skip elements since we can assume diff_bdry is
    // already zero on the internal boundaries...
    std::vector<bool> skipElem(mesh.numBoundaryElements());
    for (auto be : mesh.boundaryElements())
        skipElem.at(be.index()) = be->isInternal;

    auto M = MassMatrix::construct<1>(mesh.boundary(), /* lumped = */ false, skipElem);

    // Account for reflections... (though scale factor only affects
    // normalization and is probably irrelevant)
    M *= Real(NR);

    // Determine which variables correspond to true boundary vertices (as
    // opposed to internal periodic vertices).
    // Vertices are true boundary vertices if any non-internal boundary
    // element touches them. Only these will be free variables.
    std::vector<bool> isTrueBoundaryVertex(mesh.numBoundaryVertices(), false);
    for (auto be : mesh.boundaryElements()) {
        if (be->isInternal) continue;
        for (auto bv : be.vertices())
            isTrueBoundaryVertex.at(bv.index()) = true;
    }

    typename _Sim::VField result(mesh.numBoundaryVertices());
    result.clear();
    for (size_t c = 0; c < N; ++c) {
        SPSDSystem<Real> MSys(M);

        std::vector<size_t> fixedVars;
        // Non-true boundary vertices will have zero velocity. Also true
        // boundary vertices on the cell faces will still have zero velocity
        // perpendicular to the cell face.
        for (auto bv : mesh.boundaryVertices()) {
            if (!isTrueBoundaryVertex.at(bv.index())) {
                fixedVars.push_back(bv.index());
            }
            else {
                PeriodicBoundaryMatcher::FaceMembership<N> fm(
                        bv.volumeVertex().node()->p, mesh.boundingBox(), 1e-11);
                if (fm.onMinOrMaxFace(c)) fixedVars.push_back(bv.index());
            }
        }
        MSys.fixVariables(fixedVars, std::vector<Real>(fixedVars.size()));

        std::vector<Real> rhs(mesh.numBoundaryVertices());
        for (auto bv : mesh.boundaryVertices())
            rhs.at(bv.index()) = diff_bdry(bv.index())[c];
        auto x = MSys.solve(rhs);
        // d = -M^{-1} g
        for (auto bv : mesh.boundaryVertices())
            result(bv.index())[c] = -x.at(bv.index());
    }

    // Normalize:
    //  min_{||d||_M = 1}  g^T d ==> lambda d = -M^-1 g
    //  ||lambda d||_M = 1 ==> lambda^2 d^T M D = 1
    //  ==> lambda = sqrt(1/(d^T M d)) = 1/sqrt(-g^T d)
    result /= sqrt(-diff_bdry[result]);

    return result;
}

// Compute the steepest descent volume velocity field from a
// one-form in the volume.
template<class _Sim>
typename _Sim::VField descent_from_diff_vol(
        const ScalarOneForm<_Sim::N> &/* diff_vol */,
        const _Sim &/* sim */)
{
    throw std::runtime_error("descent_from_diff_vol unimplemented");
}

}

#endif /* end of include guard: SDCONVERSIONSORTHOCELL_HH */
