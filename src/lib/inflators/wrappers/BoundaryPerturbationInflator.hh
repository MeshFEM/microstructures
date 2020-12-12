////////////////////////////////////////////////////////////////////////////////
// BoundaryPerturbationInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Inflator taking an initial mesh and using the boundary vertex position
//      offsets as parameters.
//
//      Vertices on the periodic boundary are allowed to move, but are
//      constrained to stay on their original periodic boundaries (i.e. cannot
//      leave or enter a new one--that would mean a topology change). Also,
//      they are constrained to move consistently with their periodic pair(s).
//      WARNING: vertices that start on a corner are constrained to stay there!
//
//      Only the positions of the "true" boundary vertices are specified as
//      parameters. All other position variables are solved for using a uniform
//      graph Laplacian. "True" boundary vertices are those with at least one
//      non-periodic boundary element incident. In other words, they are
//      the vertices that are still on the boundary after the period cell's
//      identified faces have been stitched together. Note that the true
//      boundary vertices lying on P periodic boundaries (0<=P<=N) will only
//      have N-P parameters due to the periodicity constraints, and these will
//      be shared by all 2^P identified vertices.
//
//      First, periodicity constraints are enforced by constructing a reduced
//      set of variables. Then then we construct a set of parameters (one
//      parameter for each "true" boundary variable) from the variable set. The
//      i^th coordinate of each vertex is expressed in terms of variable
//      indices by array m_varForCoordinate[i], and the parameter corresponding
//      to each variable is given by m_paramForVariable[i] (which is NONE for
//      dependent variables).
//
//      Given a set of parameters, vertex coordinates are computed by:
//          For each coordinate i in 0 to N-1
//              1) Solve uniform Laplacian with parameter values (and periodic
//                 face coordinates) as Dirichlet constraints.
//              2) Decode the variables into vertex coordinates using
//                 m_paramForVariable[i]
//
//      Given a set of vertex coordinates/offsets, parameter values/offsets are
//      computed by:
//      For each coordinate i in 0 to N-1
//          1) Loop over vertices and extract variable values using
//             m_varForCoordinate[i]. Validate consistent extraction (i.e.
//             input coordinates should respectd periodic constraints)
//          2) Loop over variables and extract parameters using
//             m_paramForVariable[i]
//
//      There are 2 * N special variables, indexed 0..2N - 1 that hold the
//      periodic face coordinates:
//      0: min x, 1: min y[, 2: min z], N: max x, N + 1: max y[, N + 2: max z]
//
//  Alternate versions:
//    - Leave internal vertices unperturbed.
//    - Solve for internal vertex perturbation using Laplacian w/ Dirichlet
//      boundary conditions.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/19/2015 02:14:38
////////////////////////////////////////////////////////////////////////////////
#ifndef BOUNDARYPERTURBATIONINFLATOR_HH
#define BOUNDARYPERTURBATIONINFLATOR_HH
#include <vector>
#include <array>
#include <limits>

#include <MeshFEM/GlobalBenchmark.hh>

#include "../Inflator.hh"
#include <MeshFEM/UniformLaplacian.hh>
#include <MeshFEM/FEMMesh.hh>
#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/BoundaryConditions.hh>

template<size_t N>
class BoundaryPerturbationInflator : public Inflator<N> {
public:
    using Base = Inflator<N>;
    using Base::m_vertices;
    using Base::m_elements;

    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    using Mesh = FEMMesh<N, 1, VectorND<N>>;

    BoundaryPerturbationInflator(const std::string &meshPath,
                                 Real epsilon = 1e-5);

    BoundaryPerturbationInflator(const std::vector<MeshIO::IOVertex>  &inVertices,
                                 const std::vector<MeshIO::IOElement> &inElements,
                                 bool periodic = true,
                                 Real epsilon = 1e-5);

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override;
public:

    // Extract custom bounds given a maxim offset and the tolerance.
    // This is important for periodic cases where you don't want vertices to move outise of
    // your original base cell.
    std::vector<std::vector<Real>> customBounds(Real maxOffset = 0.1, Real tol = 1e-3) const {
        std::vector<std::vector<Real>> result(numParameters());
        Real upperBound, lowerBound;

        // initialize
        for (unsigned param = 0; param < numParameters(); param++) {
            result[param] = {-maxOffset, maxOffset};
        }

        for (size_t d = 0; d < N; ++d) {
            for (size_t vari = 0; vari < m_numVars[d]; ++vari) {
                size_t p = m_paramForVariable[d][vari];
                if (p != NONE) {
                    if ((m_origParams[p] + maxOffset) > (m_bbox.maxCorner[d] - tol)) {
                        upperBound = m_bbox.maxCorner[d] - tol - m_origParams[p];
                    }
                    else {
                        upperBound = maxOffset;
                    }

                    if ((m_origParams[p] - maxOffset) < (m_bbox.minCorner[d] + tol)) {
                        lowerBound = m_bbox.minCorner[d] + tol - m_origParams[p];
                    }
                    else {
                        lowerBound = - maxOffset;
                    }

                    result[p] = {lowerBound, upperBound};
                }
            }
        }

        return result;
    }

    // Returns volume of current mesh
    Real volume() const {
        return m_mesh->volume();
    }

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity
    ////////////////////////////////////////////////////////////////////////////
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const override {
        std::vector<VectorField<Real, N>> result(numParameters());

        for (unsigned param = 0; param < numParameters(); param++) {
            VectorField<Real, N> velocityForP(m_mesh->numVertices());
            velocityForP.clear();
            for (auto vv : m_mesh->vertices()) {
                for (size_t d = 0; d < N; ++d) {
                    auto var = m_varForCoordinate[d].at(vv.index());
                    size_t p = m_paramForVariable[d].at(var);
                    if (p == param)
                        velocityForP(vv.index())[d] = 1.0;
                    else
                        velocityForP(vv.index())[d] = 0.0;
                }
            }

            result[param] = velocityForP;
        }

        return result;
    }

    // Read off the parameter values from a particular per-boundary-vertex
    // vector field, verifying its consistency with the periodic constraints
    // If guaranteeing consistent boundary vertex enumerations across multiple
    // FEMMesh instances becomes a problem, we could change this to take a
    // per-volume-vertex field.
    virtual ScalarField<Real> paramsFromBoundaryVField(const VectorField<Real, N> &values) const override;

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; } // vertices coordinates can be considered parameters
    virtual size_t numParameters() const override { return m_numParams; }
    virtual std::vector<Real> defaultParameters() const override { return std::vector<Real>(m_numParams); }
    virtual ParameterType parameterType(size_t /* p */) const override {
        return ParameterType::Offset;
    }
    virtual bool isPrintable(const std::vector<Real> &/* p */) override {
        // TODO
        return true;
    }
    virtual void setReflectiveInflator(bool /*use*/)  override { } // needed when inflator does not implement symmetry

    ////////////////////////////////////////////////////////////////////////////
    // Boundary perturbation-specific
    ////////////////////////////////////////////////////////////////////////////
    void setNoPerturb(bool noPerturb) { m_noPerturb = noPerturb; }

    // Boundary vertex normal vector field (0 on periodic boundary) Uses area
    // for averaging.
    VectorField<Real, N> boundaryVertexNormals() const {
        VectorField<Real, N> result(m_mesh->numBoundaryVertices());
        result.clear();
        for (auto be : m_mesh->boundaryElements()) {
            if (m_isPeriodicBE.at(be.index())) continue;
            for (size_t c = 0; c < be.numVertices(); ++c)
                result(be.vertex(c).index()) += be->volume() * be->normal();
        }
        for (size_t i = 0; i < result.domainSize(); ++i) {
            Real norm = result(i).norm();
            if (norm > 1e-6)
                result(i) /= norm;
        }
        return result;
    }

    // Get the boundary vector field corresponding to "params" (i.e. the
    // inverse of paramsFromBoundaryVField)
    virtual VectorField<Real, N> boundaryVFieldFromParams(const ScalarField<Real> &params) const override;

    const Mesh &mesh() const { return *m_mesh; }

    bool isInsideBoundaryCondition(size_t vi, std::vector<CondPtr<N> > &bconds) const {
        for (CondPtr<N> bc : bconds) {
            if (bc->containsPoint(m_mesh->vertex(vi).node()->p))
                return true;
        }

        return false;
    }

    // Obtain parameter indices related to vertices around a given input point.
    std::vector<int> pointToParametersIndices(const VectorND<N> inputPoint, double tolerance = 1e-8) {
        std::vector<int> result;

        for (size_t vi = 0; vi < m_mesh->numVertices(); ++vi) {
            VectorND<N> currentPoint = m_mesh->vertex(vi).node()->p;

            if ((currentPoint - inputPoint).squaredNorm() < tolerance) {

                for (size_t d = 0; d < N; d++) {
                    size_t var = m_varForCoordinate[d].at(vi);
                    size_t p = m_paramForVariable[d].at(var);

                    if (p != NONE)
                        result.push_back(p);
                }
            }
        }

        return result;
    }

    // Obtain point related to parameter index
    VectorND<N> parameterIndexToPoint(size_t param) {
        VectorND<N> empty;

        for (auto vv : m_mesh->vertices()) {
            for (size_t d = 0; d < N; ++d) {
                auto var = m_varForCoordinate[d].at(vv.index());
                size_t p = m_paramForVariable[d].at(var);
                if (p == param)
                    return vv.node()->p;
            }
        }

        return empty;
    }

    virtual ~BoundaryPerturbationInflator() { }

private:
    // Vector of indices into the coordinate variables
    std::array<std::vector<size_t>, N> m_varForCoordinate;
    std::array<std::vector<size_t>, N> m_paramForVariable;
    std::vector<bool> m_isPeriodicBE;
    std::array<std::vector<bool>, N> m_bcVertexVariable; // say if variable is from boundary condition node
    std::array<std::vector<Real>, N> m_bcVertexValue; // saves variable constant value when inside a boundary condition region

    size_t m_numParams;
    std::array<size_t, N> m_numVars;

    // The original boundary coordinates pre-perturbation (i.e., the
    // coordinates to which offset parameters are applied)
    ScalarField<Real> m_origParams;

    // Avoid perturbing interior vertices when parameter vector is zero.
    // (Useful for remeshing gradient descent.)
    bool m_noPerturb = false;
    bool m_isPeriodicMesh = true;

    std::unique_ptr<Mesh> m_mesh;
    BBox<VectorND<N>> m_bbox;
    std::vector<CondPtr<N> > m_bconds;

    void m_setMesh(const std::vector<MeshIO::IOVertex>  &inVertices,
                   const std::vector<MeshIO::IOElement> &inElements, Real epsilon);

    // Non periodic case
    void m_setNonPeriodicMesh(const std::vector<MeshIO::IOVertex>  &inVertices,
                   const std::vector<MeshIO::IOElement> &inElements);
};

#endif /* end of include guard: BOUNDARYPERTURBATIONINFLATOR_HH */
