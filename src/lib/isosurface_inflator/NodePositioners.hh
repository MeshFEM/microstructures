////////////////////////////////////////////////////////////////////////////////
// NodePositioners.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Classes wrapping nodes' positional degrees of freedom in the symmetry
//      base unit. DoFs are constrained to keep nodes inside the base unit.
//
//      Positioners for two base units types are implemented: box (for triply
//      periodic and orthotropic symmetry), and tet (for cubic symmetry).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/20/2015 11:44:59
////////////////////////////////////////////////////////////////////////////////
#ifndef NODEPOSITIONERS_HH
#define NODEPOSITIONERS_HH

#include "InflatorTypes.hh"
#include "FuzzySign.hh"
#include <MeshFEM/Geometry.hh>

// Offsets must keep a point on the base cube region (face, edge, node, cube
// interior) on which it started.
// Position variables are always in the range [0, 1]
template<typename Real, typename TOL>
struct BoxNodePositioner {
    BoxNodePositioner(const BBox<Point3<Real>> &cell, const Point3<Real> &p) { forPoint(cell, p); }

    void forPoint(const BBox<Point3<Real>> &cell, const Point3<Real> &p) {
        for (size_t c = 0; c < 3; ++c) {
            // Node can move in any coordinate it doesn't share with a
            // min/max cell face.
            if      (isZero<TOL>(p[c] - cell.minCorner[c])) m_ctype[c] = ComponentType::MinFace;
            else if (isZero<TOL>(p[c] - cell.maxCorner[c])) m_ctype[c] = ComponentType::MaxFace;
            else                                            m_ctype[c] = ComponentType::Free;
        }
        m_dimensions = cell.dimensions();
        m_minCorner = cell.minCorner;
    }

    size_t numDoFs() const { return (m_ctype[0] == ComponentType::Free) +
                                    (m_ctype[1] == ComponentType::Free) +
                                    (m_ctype[2] == ComponentType::Free); }

    // Get the 3D position corresponding to the input degrees of freedom.
    // Supports a different type since this method will be autodiff-ed
    template<typename Real2>
    Point3<Real2> getPosition(const Real2 *dofs) const {
        Vector3<Real2> pos;
        size_t d = 0;
        for (size_t i = 0; i < 3; ++i) {
            pos[i] = m_minCorner[i];
            if (m_ctype[i] == ComponentType::Free)    pos[i] += dofs[d++] * m_dimensions[i];
            if (m_ctype[i] == ComponentType::MaxFace) pos[i] += m_dimensions[i];
        }
        assert(d == numDoFs());
        return pos;
    }

    // Get the map from [params 1] to (x, y, z) position (the last column is an
    // constant translation).
    // "paramOffset" is the index of the first position parameter for this
    // node.
    Eigen::Matrix3Xd
    getPositionMap(const size_t paramsSize, const size_t paramOffset) const {
        assert(paramOffset + numDoFs() <= paramsSize);
        Eigen::Matrix3Xd posMap(3, paramsSize + 1);
        const size_t constTransCol = paramsSize;
        posMap.setZero();
        posMap.col(constTransCol) = m_minCorner;
        size_t d = paramOffset;
        for (size_t i = 0; i < 3; ++i) {
            if (m_ctype[i] == ComponentType::   Free) posMap(i,           d++) += m_dimensions[i];
            if (m_ctype[i] == ComponentType::MaxFace) posMap(i, constTransCol) += m_dimensions[i];
        }
        assert(d == numDoFs() + paramOffset);

        return posMap;
    }

    // Get the degrees of freedom corresponding to a 3D position. Assumes that
    // the position lies within the space spanned by this positioner.
    template<typename Real2>
    void getDoFsForPoint(const Point3<Real> &p, Real2 *dofs) const {
        size_t d = 0;
        for (size_t i = 0; i < 3; ++i) {
            if (m_ctype[i] == ComponentType::Free)
                dofs[d++] = (p[i] - m_minCorner[i]) / m_dimensions[i];
            else {
                // Verify that p lies in the expected subspace
                if (m_ctype[i] == ComponentType::MinFace)
                    assert(isZero<TOL>(p[i] - m_minCorner[i]));
                if (m_ctype[i] == ComponentType::MaxFace)
                    assert(isZero<TOL>(p[i] - m_minCorner[i] - m_dimensions[i]));
            }

        }
        assert(d == numDoFs());
    }

private:
    enum class ComponentType { MinFace, MaxFace, Free };
    ComponentType m_ctype[3];
    Point3<Real> m_minCorner;
    Vector3<Real> m_dimensions;
};

// Offsets must keep a point on the base tetrahedron simplex on which it
// started (face, node, edge, tet interior). The positional degrees of
// freedom are barycentric coordinates for this simplex.
// Position variables are always in the range [0, 1], but there is an
// additional linear inequality constraint restricting their sum to 1.
template<typename Real, typename TOL>
struct TetNodePositioner {
    typedef Eigen::Matrix<Real, 4, 1> BaryCoords;

    // (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
    static int baseTetCornerPosition(size_t corner, size_t component) { return component < corner; }
    static BaryCoords barycentricCoordinates(const Point3<Real> &p) {
        return BaryCoords(1 - p[0], p[0] - p[1], p[1] - p[2], p[2]);
    }

    TetNodePositioner(const BBox<Point3<Real>> & /* cell */, const Point3<Real> &p) { forPoint(p); }
    TetNodePositioner(const Point3<Real> &p) { forPoint(p); }

    // Assumes point is within the base unit!
    void forPoint(const Point3<Real>    &p) { forPoint(barycentricCoordinates(p)); }
    void forPoint(const BaryCoords &lambda) {
        m_numAffectedBarycoords = 0;
        for (size_t c = 0; c < 4; ++c) {
            if (!isZero<TOL>(lambda[c]))
                m_affectedBaryCoordIndices[m_numAffectedBarycoords++] = c;
        }
        assert(m_numAffectedBarycoords > 0);
    }

    size_t numDoFs() const { return m_numAffectedBarycoords - 1; }

    // Get the 3D position corresponding to the input degrees of freedom.
    // Supports a different type since this method will be autodiff-ed
    template<typename Real2>
    Point3<Real2> getPosition(const Real2 *dofs) const {
        Vector3<Real2> pos(Vector3<Real2>::Zero());
        Real2 lastCoordinate = 0;
        for (size_t i = 0; i < numDoFs(); ++i) {
            size_t corner = m_affectedBaryCoordIndices[i];
            pos[0] += dofs[i] * baseTetCornerPosition(corner, 0);
            pos[1] += dofs[i] * baseTetCornerPosition(corner, 1);
            pos[2] += dofs[i] * baseTetCornerPosition(corner, 2);
            lastCoordinate += dofs[i];
        }
        if (lastCoordinate > 1.0) {
            std::cerr << "WARNING: node left base tet (lastCoordinate = " << lastCoordinate << ")" << std::endl;
            std::cerr << "DOFs: " << std::endl;
            for (size_t i = 0; i < numDoFs(); ++i)
                std::cerr << dofs[i] << "\t";
            std::cerr << std::endl;
        }
        lastCoordinate = 1.0 - lastCoordinate;
        size_t corner = m_affectedBaryCoordIndices[m_numAffectedBarycoords - 1];
        pos[0] += lastCoordinate * baseTetCornerPosition(corner, 0);
        pos[1] += lastCoordinate * baseTetCornerPosition(corner, 1);
        pos[2] += lastCoordinate * baseTetCornerPosition(corner, 2);
        return pos;
    }

    // Get the map from [params 1] to (x, y, z) position (the last column is an
    // constant translation).
    // "paramOffset" is the index of the first position parameter for this
    // node.
    // Parameters are the nonzero barycentric coordinates (except for the last,
    // redundant one.)
    Eigen::Matrix3Xd
    getPositionMap(const size_t paramsSize, const size_t paramOffset) const {
        assert(paramOffset + numDoFs() <= paramsSize);

        Eigen::Matrix3Xd posMap(3, paramsSize + 1);
        const size_t constTransCol = paramsSize;
        posMap.setZero();

        size_t lastCorner = m_affectedBaryCoordIndices[m_numAffectedBarycoords - 1];
        // last nonzero barycentric coordinate is 1.0 - sum(previous)
        posMap(0, constTransCol) = baseTetCornerPosition(lastCorner, 0);
        posMap(1, constTransCol) = baseTetCornerPosition(lastCorner, 1);
        posMap(2, constTransCol) = baseTetCornerPosition(lastCorner, 2);

        for (size_t d = 0; d < numDoFs(); ++d) {
            size_t corner = m_affectedBaryCoordIndices[d];
            posMap(0, paramOffset + d) = baseTetCornerPosition(corner, 0) - baseTetCornerPosition(lastCorner, 0);
            posMap(1, paramOffset + d) = baseTetCornerPosition(corner, 1) - baseTetCornerPosition(lastCorner, 1);
            posMap(2, paramOffset + d) = baseTetCornerPosition(corner, 2) - baseTetCornerPosition(lastCorner, 2);
        }
        return posMap;
    }

    // Get the degrees of freedom corresponding to a 3D position. Assumes that
    // the position lies within the space spanned by this positioner.
    template<typename Real2>
    void getDoFsForPoint(const Point3<Real> &p, Real2 *dofs) const {
        auto b = barycentricCoordinates(p);
        // Verify that b is in the base tet.
        for (size_t i = 0; i < 4; ++i)
            assert(isPositive<TOL>(b[i])), assert(isNegative<TOL>(b[i] - 1));
        Real simplexWeight = 0.0;
        for (size_t i = 0; i < m_numAffectedBarycoords; ++i) {
            simplexWeight += b[m_affectedBaryCoordIndices[i]];
            if (i < numDoFs()) dofs[i] = b[m_affectedBaryCoordIndices[i]];
        }
        // Verify that the point is within the expected simplex
        assert(isZero<TOL>(simplexWeight - 1.0));
    }

private:
    size_t m_numAffectedBarycoords;
    size_t m_affectedBaryCoordIndices[4];
};

// Offsets must keep a point on the base prism simplex on which it
// started (face, node, edge, prism interior). The positional degrees of
// freedom are barycentric coordinates on the base triangle + z coordinate.
//    c
//   /|
//  / |
// a  +
//  \ |
//   \|
//    b
// - If a vertex is on a corner a,b,c then it has 0 dofs on the XY plane.
// - If a vertex is on an edge, then it is constrained to stay on the edge
// - Parameters are assigned to vertices on edge (b, c) to ensure reflectional symmetry is
//   preserved when parameters are linked (by `Symmetry::independentVertexPosition()`):
//   increasing the parameter shared by the vertex pair (v_upper, v_lower)
//   moves v_upper up and v_lower down.
template<typename Real, typename TOL, bool VerticalInterfaceSymmetry = false>
struct PrismNodePositioner {
    typedef Eigen::Matrix<Real, 3, 1> BaryCoords;

    // (0, 0), (1, -1), (1, 1)
    static int baseCornerPosition(size_t corner, size_t component) {
        static int coords [] = {
            0, 0,
            1, -1,
            1, 1
        };
        return coords[2 * corner + component];
    }

    static BaryCoords barycentricCoordinates(const Point3<Real> &p) {
        return BaryCoords(1 - p[0], (p[0] - p[1]) / 2.0, (p[0] + p[1]) / 2.0);
    }

    PrismNodePositioner(const BBox<Point3<Real>> &cell, const Point3<Real> &p) { forPoint(cell, p); }

    // Assumes point is within the base unit!
    void forPoint(const BBox<Point3<Real>> &cell, const Point3<Real> &p) {
        // XY plane
        BaryCoords lambda = barycentricCoordinates(p);
        m_numAffectedBarycoords = 0;
        for (size_t c = 0; c < 3; ++c) {
            if (!isZero<TOL>(lambda[c])) {
                m_affectedBaryCoordIndices[m_numAffectedBarycoords++] = c;
            }
        }
        assert(m_numAffectedBarycoords > 0);
        m_isInterfaceMidpoint = false;
        if (VerticalInterfaceSymmetry && isZero<TOL>(lambda[0]) && (m_numAffectedBarycoords == 2)) {
            // Point is on the vertical cell interface; we need to assign the
            // positional DoF in a way that ensures vertical symmetry across
            // y = 0 is maintained when symmetric pairs are linked to the same
            // variable.
            // We do this by defining the DoF as corner (1, -1)'s barycentric
            // coordinate if y < 0 and (1, 1)'s barycentric coordinate if y > 0.
            // In other words, we choose the larger of the two barycentric
            // coordinates as the variable.
            // If y = 0 (barycentric coordinates equal), we must remove the DoF.
            Real l = lambda[m_affectedBaryCoordIndices[0]];
            if (isZero<TOL>(l - 0.5)) {
                m_numAffectedBarycoords = 1; // remove the degree of freedom
                m_isInterfaceMidpoint = true; // flag this as the interface midpoint node so that we know how to position it.
            }
            if (l < 0.5) std::swap(m_affectedBaryCoordIndices[0], m_affectedBaryCoordIndices[1]);
        }

        // Z coordinate
        {
            // Node can move in any coordinate it doesn't share with a min/max cell face.
            if      (isZero<TOL>(p[2] - cell.minCorner[2])) m_ctype = ComponentType::MinFace;
            else if (isZero<TOL>(p[2] - cell.maxCorner[2])) m_ctype = ComponentType::MaxFace;
            else                                            m_ctype = ComponentType::Free;
        }
        m_minCorner = cell.minCorner[2];
        m_extent = cell.dimensions()[2];
    }

    size_t numDoFsXY() const { return m_numAffectedBarycoords - 1; }
    size_t numDoFsZ() const { return (m_ctype == ComponentType::Free); }
    size_t numDoFs() const { return numDoFsXY() + numDoFsZ(); }

    // Get the 3D position corresponding to the input degrees of freedom.
    // Supports a different type since this method will be autodiff-ed
    template<typename MyReal>
    Point3<MyReal> getPosition(const MyReal *dofs) const {
        Point3<MyReal> pos(Vector3<MyReal>::Zero());

        // XY plane
        MyReal lastCoordinate = 0;
        for (size_t i = 0; i < numDoFsXY(); ++i) {
            size_t corner = m_affectedBaryCoordIndices[i];
            pos[0] += dofs[i] * baseCornerPosition(corner, 0);
            pos[1] += dofs[i] * baseCornerPosition(corner, 1);
            lastCoordinate += dofs[i];
        }
        if (lastCoordinate > 1.0) {
            std::cerr << "WARNING: node left base triangle (lastCoordinate = " << lastCoordinate << ")" << std::endl;
            std::cerr << "DOFs: " << std::endl;
            for (size_t i = 0; i < numDoFsXY(); ++i) {
                std::cerr << dofs[i] << "\t";
            }
            std::cerr << std::endl;
        }
        lastCoordinate = 1.0 - lastCoordinate;
        size_t corner = m_affectedBaryCoordIndices[m_numAffectedBarycoords - 1];
        pos[0] += lastCoordinate * baseCornerPosition(corner, 0);
        pos[1] += lastCoordinate * baseCornerPosition(corner, 1);

        // The above code incorrectly positions the interface midpoint at
        // corner b; manually position it.
        if (VerticalInterfaceSymmetry && m_isInterfaceMidpoint) {
            assert(numDoFsXY() == 0);
            pos = Point3<MyReal>(1, 0, 0);
        }

        // Z coordinate
        size_t d = numDoFsXY();
        {
            pos[2] = m_minCorner;
            if (m_ctype == ComponentType::Free)    pos[2] += dofs[d++] * m_extent;
            if (m_ctype == ComponentType::MaxFace) pos[2] += m_extent;
        }

        assert(d == numDoFs());
        return pos;
    }

    // Get the map from [params 1] to (x, y, z) position (the last column is an
    // constant translation).
    // "paramOffset" is the index of the first position parameter for this
    // node.
    Eigen::Matrix3Xd
    getPositionMap(const size_t paramsSize, const size_t paramOffset) const {
        assert(paramOffset + numDoFs() <= paramsSize);

        Eigen::Matrix3Xd posMap(3, paramsSize + 1);
        const size_t constTransCol = paramsSize;
        posMap.setZero();

        // XY plane
        size_t lastCorner = m_affectedBaryCoordIndices[m_numAffectedBarycoords - 1];
        // last nonzero barycentric coordinate is 1.0 - sum(previous)
        posMap(0, constTransCol) = baseCornerPosition(lastCorner, 0);
        posMap(1, constTransCol) = baseCornerPosition(lastCorner, 1);

        for (size_t d = 0; d < numDoFsXY(); ++d) {
            size_t corner = m_affectedBaryCoordIndices[d];
            posMap(0, paramOffset + d) = baseCornerPosition(corner, 0) - baseCornerPosition(lastCorner, 0);
            posMap(1, paramOffset + d) = baseCornerPosition(corner, 1) - baseCornerPosition(lastCorner, 1);
        }

        // The above code incorrectly positions the interface midpoint at
        // corner b; manually position it.
        if (VerticalInterfaceSymmetry && m_isInterfaceMidpoint) {
            assert(numDoFsXY() == 0);
            posMap(0, constTransCol) = 1;
            posMap(1, constTransCol) = 0;
        }

        // Z coordinate
        posMap(2, constTransCol) = m_minCorner;
        size_t d = paramOffset + numDoFsXY();
        {
            if (m_ctype == ComponentType::   Free) posMap(2,           d++) += m_extent;
            if (m_ctype == ComponentType::MaxFace) posMap(2, constTransCol) += m_extent;
        }
        assert(d == numDoFs() + paramOffset);

        return posMap;
    }

    // Get the degrees of freedom corresponding to a 3D position. Assumes that
    // the position lies within the space spanned by this positioner.
    template<typename Real2>
    void getDoFsForPoint(const Point3<Real> &p, Real2 *dofs) const {
        // XY plane
        auto b = barycentricCoordinates(p);
        // Verify that b is in the base triangle.
        for (size_t i = 0; i < 3; ++i) {
            assert(isPositive<TOL>(b[i])), assert(isNegative<TOL>(b[i] - 1));
        }
        Real simplexWeight = 0.0;
        for (size_t i = 0; i < m_numAffectedBarycoords; ++i) {
            simplexWeight += b[m_affectedBaryCoordIndices[i]];
            if (i < numDoFsXY()) {
                dofs[i] = b[m_affectedBaryCoordIndices[i]];
            }
        }
        // Verify that the point is within the expected simplex
        assert(isZero<TOL>(simplexWeight - 1.0));

        // Z coordinate
        size_t d = numDoFsXY();
        {
            if (m_ctype == ComponentType::Free) {
                dofs[d++] = (p[2] - m_minCorner) / m_extent;
            } else {
                // Verify that p lies in the expected subspace
                if (m_ctype == ComponentType::MinFace) {
                    assert(isZero<TOL>(p[2] - m_minCorner));
                }
                if (m_ctype == ComponentType::MaxFace) {
                    assert(isZero<TOL>(p[2] - m_minCorner - m_extent));
                }
            }

        }
        assert(d == numDoFs());
    }

private:
    // Barycentric coordinate for (x, y) position
    size_t m_numAffectedBarycoords;
    size_t m_affectedBaryCoordIndices[4];

    // When enforcing VerticalInterfaceSymmetry, the midpoint has no variables
    // and must be positioned manually.
    bool   m_isInterfaceMidpoint;

    // For the z coordinate
    enum class ComponentType { MinFace, MaxFace, Free };
    ComponentType m_ctype;
    Real m_minCorner;
    Real m_extent;
};

#endif /* end of include guard: NODEPOSITIONERS_HH */
