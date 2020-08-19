////////////////////////////////////////////////////////////////////////////////
// ShapeVelocityInterpolator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      While the change in shape is determined entirely by the boundary
//      (normal) velocity field, it is often more accurate for the purposes of
//      discrete shape optimization to consider smooth perturbations of the
//      whole mesh instead of just the boundary. To do this, we interpolate the
//      boundary velocity by solving a Laplace equation.
//
//      To maintain periodicity, this interpolated shape velocity must be
//      periodic (but we don't need its boundary to stay clamped to the square
//      periodic cell). We achieve this by enforcing periodic boundary
//      conditions in the Laplace solve on the period cell.
//
//      Shape velocities of the "internal periodic vertices" (mesh boundary
//      vertices that would become internal after the mesh is periodically
//      tiled) are not read/applied as Dirichlet constraints (they should be
//      zero anyway, but the inflator is not required to guarantee this).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  03/19/2016 21:13:18
////////////////////////////////////////////////////////////////////////////////
#ifndef SHAPEVELOCITYINTERPOLATOR_HH
#define SHAPEVELOCITYINTERPOLATOR_HH

#include <MeshFEM/Laplacian.hh>
#include <limits>

#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <MeshFEM/LinearIndexer.hh>

class ShapeVelocityInterpolator {
    static constexpr size_t NO_VAR = std::numeric_limits<size_t>::max(),
                              NONE = std::numeric_limits<size_t>::max();
public:
    template<class Sim>
    ShapeVelocityInterpolator(const Sim &sim) {
        const auto &mesh = sim.mesh();
        L = Laplacian::construct<1>(mesh);

        // Apply same periodic boundary conditions as sim.
        // Unfortunately, sim's periodic boundary conditions are on nodes, not
        // just vertices, so we have to reindex to get contiguous variables
        size_t nv = mesh.numVertices();
        m_varForVertex.reserve(nv);
        std::vector<size_t> varForDoF(sim.numDoFs(), size_t(NO_VAR));
        size_t numVars = 0;
        for (auto v : mesh.vertices()) {
            size_t d = sim.DoF(v.node().index());
            size_t &var = varForDoF.at(d);
            if (var == NO_VAR)
                var = numVars++;
            m_varForVertex.push_back(var);
        }

        // Apply periodic boundary conditions to L
        L.reindexVariables(numVars, m_varForVertex);

        // Determine which variables correspond to true boundary vertices (as
        // opposed to internal periodic vertices).
        // Vertices are true boundary vertices if any non-internal boundary
        // element touches them.
        std::vector<bool> isTrueBoundaryVertex(mesh.numBoundaryVertices(), false);
        for (auto be : mesh.boundaryElements()) {
            if (be->isInternal) continue;
            for (auto bv : be.vertices())
                isTrueBoundaryVertex.at(bv.index()) = true;
        }

        // Determine vars corresponding to true boundary vertices (to be constrained)
        std::vector<size_t> bdryVarIdxForVar(numVars, size_t(NONE)); // Has a variable has been marked as bdry?
        m_bdryVarIdxForBdryVtx.assign(mesh.numBoundaryVertices(), size_t(NONE));
        m_numBdryVtxsPerBdryVar.clear();

        for (auto bv : mesh.boundaryVertices()) {
            if (isTrueBoundaryVertex.at(bv.index())) {
                size_t vari = m_varForVertex.at(bv.volumeVertex().index());
                if (bdryVarIdxForVar.at(vari) == NO_VAR) {
                    // vari is newly discovered as boundary variable.
                    bdryVarIdxForVar.at(vari) = m_bdryVars.size();
                    m_bdryVars.push_back(vari);
                    m_numBdryVtxsPerBdryVar.push_back(0);
                }
                size_t bvarIdx = bdryVarIdxForVar.at(vari);
                ++m_numBdryVtxsPerBdryVar.at(bvarIdx);
                m_bdryVarIdxForBdryVtx.at(bv.index()) = bvarIdx;
            }
        }
    }

    // Periodically interpolate per-boundary-vertex shape velocity vector field,
    // bdry_svel.
    // WARNING: non-periodic input fields will be made periodic by choosing the
    // value from an arbitrary one of the periodically identified vertices.
    template<class Sim>
    typename Sim::VField interpolate(const Sim &sim,
                                     const typename Sim::VField &bdry_svel) const {
        constexpr size_t N = Sim::N;
        const auto &mesh = sim.mesh();

        typename Sim::VField result;
        result.resizeDomain(mesh.numVertices()); // zero-inits

        std::vector<Real> bdryVarValues, zero(L.m), x;
        // Interpolate each component of the boundary velocity.
        for (size_t c = 0; c < N; ++c) {
            bdryVarValues.assign(m_bdryVars.size(), 0.0);
            // Extract boundary variable values from the per-bdry-vertex field
            for (auto bv : mesh.boundaryVertices()) {
                size_t bvarIdx = m_bdryVarIdxForBdryVtx.at(bv.index());
                if (bvarIdx != NONE) bdryVarValues.at(bvarIdx) = bdry_svel(bv.index())[c];
            }

            // TODO: allow changing the fixed variable constraint RHS without
            // rebuilding the system
            SPSDSystem<Real> Lsys(L);
            Lsys.fixVariables(m_bdryVars, bdryVarValues);
            Lsys.solve(zero, x);
            for (size_t vi = 0; vi < mesh.numVertices(); ++vi)
                result(vi)[c] = x.at(m_varForVertex.at(vi));
        }

        return result;
    }

    // Find the discrete boundary velocity one-form, dvb, that satisfies
    //    dvb[vb] := dv[interpolate(vb)]
    // given the volume velocity one-form dv, for all periodic boundary
    // velocity fields vb.
    // In other words, applies [-L_bi L_ii^{-1} I] (the transpose of the matrix
    // applied for interpolation) to discrete volume velocity one-form dv.
    // Periodicity details:
    //    dv is non-periodic; it contains partial contributions on each
    //    identified periodic vertex (i.e. it gives the correct answer when
    //    dotted with a *periodic* nodal vector field).
    //    We introduce the following matrices mapping between vertex/dof field types:
    //       S^T: Sum identified vertices onto DoFs,                 S: copy DoFs to identified vertices
    //         B: Extract boundary DoFs from all DoFs,             B^T: distribute boundary DoFs to full DoF vector (zeros on internal DoFs)
    //         C: Extract internal DoFs from all DoFs,             C^T: distribute internal DoFs to full DoF vector (zeros on boundary DoFs)
    //         A: Average identified bdry vertices onto bdry DoFs, A^T: fractional distribution of bdry DoFs to identified bdry vertices
    //    dvb[vb] := dv[interpolate(vb)] = dv[S [-L_bi L_ii^-1 I]^T A vb]
    //             = lambda . vb
    //    ==> dvb = lambda = A^T [-L_bi L_ii^-1 I] S^T dv = A^T B  S^T dv - A^T L_bi L_ii^{-1} C S^T dv
    //                                                    = A^T B (S^T dv - L C^T    L_ii^{-1} C S^T dv)
    // The final step of applying A^T means that we create a periodic output
    // boundary field with 1/N of the contribution to each of the N identified
    // vertices.
    template<class Sim, class T = Real>
    OneForm<T, Sim::N> adjoint(const Sim &sim, const OneForm<T, Sim::N> &dv) const {
        assert(dv.domainSize() == sim.mesh().numVertices());
        using OF = OneForm<T, Sim::N>;

        // If there are no true boundary vertices, there is no adjoint velocity
        if (m_bdryVars.size() == 0) {
            OF result(sim.mesh().numBoundaryVertices());
            result.clear();
            return result;
        }
        SPSDSystem<Real> Lsys(L);

        // Fix boundary vars to zero
        // Now Lsys.solve() actually applies C^T L_ii^{-1} C
        Lsys.fixVariables(m_bdryVars, std::vector<Real>(m_bdryVars.size()));

        std::vector<Real> S_t_dv, C_t_Lii_inv_C_S_t_dv;
        OF dvb(sim.mesh().numBoundaryVertices());
        dvb.clear();
        for (size_t c = 0; c < Sim::N; ++c) {
            // Interpolate each component of the output tensor (e.g. just one
            // component for scalar-valued one-forms).
            using LI = LinearIndexer<T>;
            for (size_t cc = 0; cc < LI::size(); ++cc) {
                // Sum identified values onto the DoFs
                S_t_dv.assign(L.m, 0.0);
                for (size_t i = 0; i < m_varForVertex.size(); ++i)
                    S_t_dv.at(m_varForVertex[i]) += LI::index(dv(i)[c], cc);

                Lsys.solve(S_t_dv, C_t_Lii_inv_C_S_t_dv);
                auto dofValues = L.apply(C_t_Lii_inv_C_S_t_dv);
                // dofValues = S^T dv - L C^T    L_ii^{-1} C S^T dv
                for (size_t i = 0; i < dofValues.size(); ++i)
                    dofValues[i] = S_t_dv[i] - dofValues[i];

                // Apply A^T B: read (fractional) value for each boundary vertex.
                for (auto bv : sim.mesh().boundaryVertices()) {
                    size_t bvarIdx = m_bdryVarIdxForBdryVtx.at(bv.index());
                    if (bvarIdx == NONE) continue;
                    LI::index(dvb(bv.index())[c], cc) = dofValues.at(m_bdryVars.at(bvarIdx))
                                                      / Real(m_numBdryVtxsPerBdryVar.at(bvarIdx));
                }
            }
        }

        return dvb;
    }

private:
    // Periodic Laplacian
    TripletMatrix<> L;
    std::vector<size_t> m_varForVertex,
                        m_bdryVars,              // Indices of variables corresponding to true boundary vertices
                        m_bdryVarIdxForBdryVtx,  // Index into m_bdryVars corresponding to a boundary vertex
                        m_numBdryVtxsPerBdryVar; // Number of boundary vertices linked to the var m_bdryVars[i]
};

#endif /* end of include guard: SHAPEVELOCITYINTERPOLATOR_HH */
