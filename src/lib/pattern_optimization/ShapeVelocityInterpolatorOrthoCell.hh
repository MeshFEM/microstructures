////////////////////////////////////////////////////////////////////////////////
// ShapeVelocityInterpolatorOrthoCell.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Just like ShapeVelocityInterpolator.hh but for orthotropic base cells.
//      Here, reflectional symmetry of the geometry (which must be respected by
//      all velocity fields) requires that the velocity component perpendicular
//      to each face be zero. However, unlike the triply periodic base cell
//      case, there are no periodicity conditions and no partial contributions.
//
//      There are no explicit periodicity conditions for
//      the interpolation, though.
//
//      Again, shape velocities of the "internal mesh boundary nodes" are not
//      read/applied as Dirichlet  constraints (they should be zero anyway, but
//      the inflator is not required to guarantee this).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  08/26/2016 10:59:57
////////////////////////////////////////////////////////////////////////////////
#ifndef SHAPEVELOCITYINTERPOLATORORTHOCELL_HH
#define SHAPEVELOCITYINTERPOLATORORTHOCELL_HH

#include <MeshFEM/Laplacian.hh>
#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <MeshFEM/LinearIndexer.hh>
#include <MeshFEM/PeriodicBoundaryMatcher.hh>
#include <limits>

template<size_t N>
class ShapeVelocityInterpolatorOrthoCell {
    using FM = PeriodicBoundaryMatcher::FaceMembership<N>;
public:
    // WARNING: assumes periodic homogenization has been run on sim so that the
    // be->internal field has been written.
    template<class Sim>
    ShapeVelocityInterpolatorOrthoCell(const Sim &sim) {
        static_assert(Sim::N == N, "Incorrect simulator dimension");
        const auto &mesh = sim.mesh();
        L = Laplacian::construct<1>(mesh);

        // Determine the boundary vertices' cell face memberships--the velocity
        // components perpendicular to the cell faces will need to be pinned.
        const auto &cell = mesh.boundingBox();
        for (auto bv : mesh.boundaryVertices())
            m_bdryVtxFM.emplace_back(bv.volumeVertex().node()->p, cell, 1e-11);

        // Determine which boundary vertices are on the true boundary.
        // Vertices are true boundary vertices if any non-internal boundary
        // element touches them.
        m_isTrueBoundaryVertex.assign(mesh.numBoundaryVertices(), false);
        for (auto be : mesh.boundaryElements()) {
            if (be->isInternal) continue;
            for (auto bv : be.vertices())
                m_isTrueBoundaryVertex.at(bv.index()) = true;
        }
    }

    // Smoothly interpolate per-boundary-vertex shape velocity vector field,
    // bdry_svel.
    template<class Sim>
    typename Sim::VField interpolate(const Sim &sim,
                                     const typename Sim::VField &bdry_svel) const {
        static_assert(Sim::N == N, "Incorrect simulator dimension");
        const auto &mesh = sim.mesh();
        const size_t nv = mesh.numVertices();

        typename Sim::VField result;
        result.resizeDomain(nv); // zero-inits

        // Interpolate each component of the boundary velocity.
        for (size_t c = 0; c < N; ++c) {
            size_t numInvalidComponents = 0;

            std::unique_ptr<SPSDSystem<Real>> Lsys;
            std::vector<size_t> fixedVars;
            std::tie(Lsys, fixedVars) = m_buildSystem(c, sim);

            // Extract boundary variable values from the per-bdry-vertex field
            std::vector<Real> fixedVarValues;
            for (const size_t vvi : fixedVars) {
                const size_t bvi = mesh.vertex(vvi).boundaryVertex().index();
                Real val = bdry_svel(bvi)[c];
                if (m_bdryVtxFM.at(bvi).onMinOrMaxFace(c)) {
                    if (std::abs(val) > 1e-9) ++numInvalidComponents;
                    val = 0.0;
                }
                fixedVarValues.push_back(val);
            }

            if (numInvalidComponents) {
                std::cerr << "WARNING: " << numInvalidComponents
                          << " cellface-perpendicular components are nonzero"
                          << std::endl;
            }

            Lsys->fixVariables(fixedVars, fixedVarValues);
            auto x = Lsys->solve(std::vector<Real>(nv, 0.0));

            for (size_t vi = 0; vi < nv; ++vi)
                result(vi)[c] = x.at(vi);
        }

        return result;
    }

    // Find the discrete boundary velocity one-form, dvb, that satisfies
    //    dvb[vb] := dv[interpolate(vb)]
    // given the volume velocity one-form dv, for all periodic boundary
    // velocity fields vb.
    // In other words, applies [-(L_bi L_ii^{-1}) I] (the transpose of the matrix
    // applied for interpolation) to discrete volume velocity one-form dv.
    // Here i ranges over the internal vertex coordinate variables (omitting those
    // that are fixed to zero)
    // Conceptual matrices:
    //    B: Extract true boundary vars from all vars,  B^T: distribute boundary vars to full var vector (zeros on internal vars)
    //    C: Extract      internal vars from all vars,  C^T: distribute internal vars to full var vector (zeros on boundary vars)
    template<class Sim, class T = Real>
    OneForm<T, N> adjoint(const Sim &sim, const OneForm<T, N> &dv) const {
        static_assert(Sim::N == N, "Incorrect simulator dimension");
        const auto &mesh = sim.mesh();
        const size_t nv = mesh.numVertices();
        assert(dv.domainSize() == nv);

        OneForm<T, N> dvb(mesh.numBoundaryVertices());
        dvb.clear();

        for (size_t c = 0; c < N; ++c) {
            std::unique_ptr<SPSDSystem<Real>> Lsys;
            std::vector<size_t> fixedVars;
            std::tie(Lsys, fixedVars) = m_buildSystem(c, sim);

            // All fixed vars in the forward interpolation system are fixed to
            // zero in the adjoint system (then Lsys.solve() will actually apply
            // C^T L_ii^{-1} C)
            Lsys->fixVariables(fixedVars, std::vector<Real>(fixedVars.size()));

            // Interpolate cc-th component of output tensor (just one component
            // for scalar-valued one-forms)
            using LI = LinearIndexer<T>;
            for (size_t cc = 0; cc < LI::size(); ++cc) {
                std::vector<Real> dv_c_cc(nv);
                for (size_t i = 0; i < nv; ++i) dv_c_cc[i] = LI::index(dv(i)[c], cc);

                // internalContrib := L C^T L_ii^{-1} C dv
                auto internalContrib = L.apply(Lsys->solve(dv_c_cc));

                // Compute B(-internalContrib + dv) = [-(L_bi L_ii^{-1})   I] dv
                for (size_t vi : fixedVars) {
                    size_t bvi = mesh.vertex(vi).boundaryVertex().index();

                    // fixedVars consists of true boundary vars + perpendicular cell
                    // components. We only want to compute the differential for true
                    // boundary vertices, leaving the other components zero.
                    if (m_bdryVtxFM.at(bvi).onMinOrMaxFace(c)) continue;

                    LI::index(dvb(bvi)[c], cc) = dv_c_cc[vi] - internalContrib[vi];
                }
            }
        }

        return dvb;
    }

private:
    // Laplacian without boundary conditions.
    TripletMatrix<> L;

    // Build a new SPSDSystem laplacian (not fixing anything yet), and determine the
    // vertices whose coresponding vertices must be fixed for the cth-coordinate
    // interpolation.
    template<class Sim>
    std::pair<std::unique_ptr<SPSDSystem<Real>>, std::vector<size_t>>
    m_buildSystem(size_t c, const Sim &sim) const {
        std::vector<size_t> fixedVertices;
        for (auto bv : sim.mesh().boundaryVertices()) {
            const size_t bvi = bv.index();
            const size_t vvi = bv.volumeVertex().index();
            if (m_isTrueBoundaryVertex.at(bvi) || m_bdryVtxFM.at(bvi).onMinOrMaxFace(c))
                fixedVertices.push_back(vvi);
        }
        return make_pair(Future::make_unique<SPSDSystem<Real>>(L),
                         std::move(fixedVertices));
    }

    std::vector<FM>   m_bdryVtxFM;
    std::vector<bool> m_isTrueBoundaryVertex;
};

#endif /* end of include guard: SHAPEVELOCITYINTERPOLATORORTHOCELL_HH */
