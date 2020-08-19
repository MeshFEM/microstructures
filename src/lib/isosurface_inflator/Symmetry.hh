////////////////////////////////////////////////////////////////////////////////
// Symmetry.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Holds symmetry classes to be used as template parameters for the
//      WireMesh class. These provide the queries needed to determine the base
//      unit for a particular tiling symmetry and the degrees of freedom
//      parametrizing the symmetry-preserving subspace of patterns. They can
//      also generate (a finite subset of) the elements of the corresponding
//      symmetry group.
//
//      All symmetry queries are performed with a tolerance that is configurable
//      by template parameter.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/26/2015 17:55:53
////////////////////////////////////////////////////////////////////////////////
#ifndef SYMMETRY_HH
#define SYMMETRY_HH

#include "InflatorTypes.hh"
#include "Isometries.hh"
#include "NodePositioners.hh"
#include "FuzzySign.hh"
#include "AutomaticDifferentiation.hh"
#include <MeshFEM/Geometry.hh>
#include <ratio>
#include <vector>
#include <memory>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace Symmetry {

// Enum value is the node's number of positional degrees of freedom.
enum class NodeType : unsigned int { Vertex = 0, Edge = 1, Face = 2, Interior = 3 };

// Forward declarations of Symmetry types.
template<typename TOL = DEFAULT_TOL, size_t N = 3> struct NonPeriodic;
template<typename TOL = DEFAULT_TOL> struct TriplyPeriodic;
template<typename TOL = DEFAULT_TOL> struct DoublyPeriodic;
template<typename TOL = DEFAULT_TOL> struct Orthotropic;
template<typename TOL = DEFAULT_TOL> struct Diagonal;
template<typename TOL = DEFAULT_TOL> struct Cubic;
template<typename TOL = DEFAULT_TOL> struct Square;
template<typename TOL = DEFAULT_TOL> struct Null;

// We need a traits class for CRTP to look up the correct NodePositioner class.
// This traits class must be specialized for each symmetry type.
// TriplyPeriodic and Orthotropic symmetries have a box base cell,
// Cubic and Diagonal symmetries have a tet base cell.
template<class Sym> struct SymmetryTraits { };
template<typename TOL> struct SymmetryTraits<NonPeriodic<TOL, 2>> { template<typename Real> using NodePositioner = BoxNodePositioner<Real, TOL>; };
template<typename TOL> struct SymmetryTraits<NonPeriodic<TOL, 3>> { template<typename Real> using NodePositioner = BoxNodePositioner<Real, TOL>; };
template<typename TOL> struct SymmetryTraits<TriplyPeriodic<TOL>> { template<typename Real> using NodePositioner = BoxNodePositioner  <Real, TOL>; };
template<typename TOL> struct SymmetryTraits<DoublyPeriodic<TOL>> { template<typename Real> using NodePositioner = BoxNodePositioner  <Real, TOL>; };
template<typename TOL> struct SymmetryTraits<Orthotropic<TOL>>    { template<typename Real> using NodePositioner = BoxNodePositioner  <Real, TOL>; };
template<typename TOL> struct SymmetryTraits<Diagonal<TOL>>       { template<typename Real> using NodePositioner = PrismNodePositioner<Real, TOL, true>; };
template<typename TOL> struct SymmetryTraits<Cubic<TOL>>          { template<typename Real> using NodePositioner = TetNodePositioner  <Real, TOL>; };
template<typename TOL> struct SymmetryTraits<Square<TOL>>         { template<typename Real> using NodePositioner = TetNodePositioner  <Real, TOL>; };

// Implements some of the shared interface of the symmetry classes
template<class Sym>
struct SymmetryCRTP {
    template<typename Real>
    using NodePositioner = typename SymmetryTraits<Sym>::template NodePositioner<Real>;

    template<typename Real>
    static NodePositioner<Real> nodePositioner(const Point3<Real> &p) {
        assert(Sym::inBaseUnit(p));
        return NodePositioner<Real>(Sym::template representativeMeshCell<Real>(), p);
    }

    template<typename Real>
    static NodeType nodeType(const Point3<Real> &p) {
        return static_cast<NodeType>(nodePositioner(p).numDoFs());
    }
};

template<typename T, bool _isAutoDiffType = IsAutoDiffType<T>::value>
struct OptionalFMod2 {
    static void run(T &val) {
        T q = (int) (val / 2); // round toward zero
        // if (abs(q) > 1000) {
        //     std::cout << "queried far point " << p.transpose() << std::endl;
        // }
        val -= 2 * q; // p[c] reduced to (-2, 2)
    }
};

template<typename T>
struct OptionalFMod2<T, true> {
    static void run(T &/* val */) { /* nop */ }
};

////////////////////////////////////////////////////////////////////////////////
// Symmetry class definitions
////////////////////////////////////////////////////////////////////////////////
// Base unit is a full period cell: [-1, 1]^3
template<typename TOL, size_t N>
struct NonPeriodic : SymmetryCRTP<NonPeriodic<TOL, N>> {
    typedef TOL Tolerance;
    // Disambiguate CRTP instances
    typedef SymmetryCRTP<NonPeriodic<TOL, N>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    static constexpr double tolerance = double(TOL::num) / double(TOL::den);

    // TODO: where is this used exactly and how this changes execution? Should I scale every non periodic structure to the cube position?
    // In non periodic structures, it does not matter what is the representative mesh cell.
    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() {
        float max = 1.0; //changes performance, since resolution seems to be kept the same

        if (N == 2) {
            return BBox<Point3<Real>>(Point3<Real>(-max, -max, 0),
                                      Point3<Real>(max, max, 0));
        }
        else {
            return BBox<Point3<Real>>(Point3<Real>(-max, -max, -max),
                                      Point3<Real>(max, max, max));
        }
    }

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        return p;
    }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &/*p*/) {
        return true;
    }

    template<typename Real>
    static bool inMeshingCell(const Point3<Real> &p) {
        return inBaseUnit(p);
    }

    template<typename Real>
    static Point3<Real> independentVertexPosition(Point3<Real> p) {
        return p;
    }

    // Only identity is used for non periodic wires
    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group(1, Isometry()); // Add identity element
        return group;
    }
};

// Base unit is a full period cell: [-1, 1]^3
template<typename TOL>
struct TriplyPeriodic : SymmetryCRTP<TriplyPeriodic<TOL>> {
    typedef TOL Tolerance;
    // Disambiguate CRTP instances
    typedef SymmetryCRTP<TriplyPeriodic<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    static constexpr double tolerance = double(TOL::num) / double(TOL::den);
    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() {
        return BBox<Point3<Real>>(Point3<Real>(-1, -1, -1),
                                  Point3<Real>(1, 1, 1));
    }

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        for (size_t c = 0; c < 3; ++c) {
            // Quickly reduce to [-2, 2] for plain scalar types to accelerate
            // far point sampling.
            // This integer conversion-based reduction will not work with
            // autodiff types, but the subsequent brute-force reduction to
            // [-1, 1] should be fast since we should never be querying outside
            // the base cell when using autodiff to compute shape
            // velocities/normals).
            OptionalFMod2<Real>::run(p[c]);

            // Now reduce to [-1, 1] with tolerance...
            while (p[c] >  1.0 + tolerance) p[c] -= 2.0;
            while (p[c] < -1.0 - tolerance) p[c] += 2.0;
        }
        return p;
    }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        return (isNegative<TOL>(std::abs(p[0]) - 1.0)) &&
               (isNegative<TOL>(std::abs(p[1]) - 1.0)) &&
               (isNegative<TOL>(std::abs(p[2]) - 1.0));
    }

    template<typename Real>
    static bool inMeshingCell(const Point3<Real> &p) {
        return inBaseUnit(p);
    }

    // Find the location of the independent vertex linked to p. For vertices in
    // the base cell's interior, this is just the vertex position itself. For
    // vertices on the period cell face(s), this is the corresponding location
    // on the minimum face(s).
    template<typename Real>
    static Point3<Real> independentVertexPosition(Point3<Real> p) {
        assert(inBaseUnit(p));
        if (isZero<TOL>(std::abs(p[0] - 1.0))) p[0] = -1.0;
        if (isZero<TOL>(std::abs(p[1] - 1.0))) p[1] = -1.0;
        if (isZero<TOL>(std::abs(p[2] - 1.0))) p[2] = -1.0;
        return p;
    }

    // Note that the group of translational symmetries is infinite, but for our
    // purposes (determining incident edges from neighboring cells), the
    // isometries that take us to adjacent cells are sufficient
    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group(1, Isometry()); // Start with identity element
        // Add in translational symmetries getting us to the 26 adjacent cells
        for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
        for (int z = -1; z <= 1; ++z) {
            if ((x == 0) && (y == 0) && (z == 0)) continue;
            // Note that the base cell is size 2 ([-1, 1]^3)
            group.emplace_back(Isometry::translation(2.0 * x, 2.0 * y, 2.0 * z));
        }
        return group;
    }
};

// Hack: same as TriplyPeriodic (since the isometries involving z don't have
// any effect on the z = 0 signed distance slice), except we don't want to
// assign z translation parameters. By setting the representative cell to have zero
// height, points will be constrained to stay on the midplane.
template<typename TOL>
struct DoublyPeriodic : public TriplyPeriodic<TOL>, SymmetryCRTP<DoublyPeriodic<TOL>> {
    using CRTP = SymmetryCRTP<DoublyPeriodic<TOL>>;
    using CRTP::nodePositioner;
    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() {
        return BBox<Point3<Real>>(Point3<Real>(-1, -1, 0),
                                  Point3<Real>(1, 1, 0));
    }
};

// Base unit is the positive octant: [0, 1]^3
// Symmetry group D_2h x Translations
template<typename TOL>
struct Orthotropic : public TriplyPeriodic<TOL>, SymmetryCRTP<Orthotropic<TOL>> {
    typedef TOL Tolerance;
    // Disambiguate CRTP instances
    typedef SymmetryCRTP<Orthotropic<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    using TriplyPeriodic<TOL>::tolerance;

    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() {
        return BBox<Point3<Real>>(Point3<Real>(0, 0, 0),
                                  Point3<Real>(1, 1, 1));
    }

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        p = TriplyPeriodic<TOL>::mapToBaseUnit(p);
        for (size_t c = 0; c < 3; ++c)
            if (p[c] < 0) p[c] = -p[c]; // std::abs is problematic for autodiff
        return p;
    }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        return TriplyPeriodic<TOL>::inBaseUnit(p) &&
                isPositive<TOL>(p[0]) && isPositive<TOL>(p[1]) &&
                isPositive<TOL>(p[2]);
    }

    template<typename Real>
    static bool inMeshingCell(const Point3<Real> &p) {
        return inBaseUnit(p);
    }

    // All vertices in the orthotropic base unit are independent.
    template<typename Real>
    static Point3<Real> independentVertexPosition(Point3<Real> p) {
        assert(inBaseUnit(p));
        return p;
    }

    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group;
        std::vector<Isometry> parentGroup = TriplyPeriodic<TOL>::symmetryGroup();
        for (const Isometry &p : parentGroup) {
            group.push_back(p); // Identity reflection

            // Single axis
            group.push_back(p.compose(Isometry::reflection(Axis::X)));
            group.push_back(p.compose(Isometry::reflection(Axis::Y)));
            group.push_back(p.compose(Isometry::reflection(Axis::Z)));

            // Double axis
            group.push_back(p.compose(Isometry::reflection(Axis::X)).compose(Isometry::reflection(Axis::Y)));
            group.push_back(p.compose(Isometry::reflection(Axis::X)).compose(Isometry::reflection(Axis::Z)));
            group.push_back(p.compose(Isometry::reflection(Axis::Y)).compose(Isometry::reflection(Axis::Z)));

            // Triple axis
            group.push_back(p.compose(Isometry::reflection(Axis::X)).compose(Isometry::reflection(Axis::Y)).compose(Isometry::reflection(Axis::Z)));
        }

        return group;
    }
};

// Base unit the triangle (0, 0), (1, -1), (1, 1)
// Symmetry group ??? x Translations
// We mesh the half-space (x >= 0) since meshing within a prism is harder.
template<typename TOL>
struct Diagonal : public DoublyPeriodic<TOL>, SymmetryCRTP<Diagonal<TOL>> {
    typedef TOL Tolerance;

    // Disambiguate CRTP instances
    typedef SymmetryCRTP<Diagonal<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    using DoublyPeriodic<TOL>::tolerance;

    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() {
        return BBox<Point3<Real>>(Point3<Real>(-1, -1, 0),
                                  Point3<Real>(1, 1, 0));
    }

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        p = DoublyPeriodic<TOL>::mapToBaseUnit(p);
        if (p[0] < p[1]) std::swap(p[0], p[1]);
        if (p[1] < -p[0]) {
            Real tmp = p[1];
            p[1] = -p[0];
            p[0] = -tmp;
        }
        return p;
    }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        if (DoublyPeriodic<TOL>::inBaseUnit(p)) {
            return isPositive<TOL>(p[0]) && (p[0] + tolerance >= p[1]) && (p[0] + tolerance >= -p[1]);
        }
        return false;
    }

    template<typename Real>
    static bool inMeshingCell(const Point3<Real> &p) {
        return DoublyPeriodic<TOL>::inBaseUnit(p);
    }

    // Find the location of the independent vertex linked to p. For vertices in
    // the base cell's interior, this is just the vertex position itself. For
    // vertices on the period cell face(s), the diagonal symmetries impose that
    // the interface should have a reflective symmetry along the X and Y axes.
    template<typename Real>
    static Point3<Real> independentVertexPosition(Point3<Real> p) {
        assert(inBaseUnit(p));
        if (isZero<TOL>(std::abs(p[0] - 1.0))) {
            if (p[1] < 0) { p[1] = -p[1]; }
        }
        return p;
    }

    // We need augment DoublyPeriodic's symmetry group with the operations
    // taking region 1 into 2, 3, and 4 by reflections across the diagonals:
    // +---+
    // |\2/|
    // |3*1|
    // |/4\|
    // +---+
    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group;
        std::vector<Isometry> parentGroup = DoublyPeriodic<TOL>::symmetryGroup();
        for (const Isometry &p : parentGroup) {
            if (!p.affectsAxis(Axis::Z)) {
                group.push_back(p); // Identity (stay in region 1)

                group.push_back(p.compose(Isometry::permutation(Axis::X, Axis::Y))); // Region 1 to region 2 (swap x, y)
                group.push_back(p.compose(Isometry:: reflection(Axis::X))            // Region 1 to region 3: rotation by 180 degrees
                                 .compose(Isometry:: reflection(Axis::Y)));
                group.push_back(p.compose(Isometry:: reflection(Axis::X))            // Region 1 to region 4 (swap -x, y)
                                 .compose(Isometry::permutation(Axis::X, Axis::Y))   // (swap x, y in reflected space, transform back)
                                 .compose(Isometry:: reflection(Axis::X)));
            }
        }

        return group;
    }
};

// Base unit is one of the 6 tetrahedra in the positive octant:
// the tetrahedron with corners (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
// However, we still mesh the positive octant since meshing within a tetrahedron
// is harder.
// Symmetry group Oh x Translations
template<typename TOL>
struct Cubic : public Orthotropic<TOL>, SymmetryCRTP<Cubic<TOL>> {
    typedef TOL Tolerance;
    using Orthotropic<TOL>::representativeMeshCell;
    using TriplyPeriodic<TOL>::tolerance;

    // Disambiguate CRTP instances
    typedef SymmetryCRTP<Cubic<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        p = Orthotropic<TOL>::mapToBaseUnit(p);
        // Sort components descending: an optimal algorithm would still need 3
        // comparisons in the worst case (though 2 in the best), so this bubble
        // sort isn't too bad.
        if (p[0] < p[1]) std::swap(p[0], p[1]);
        if (p[1] < p[2]) std::swap(p[1], p[2]);
        if (p[0] < p[1]) std::swap(p[0], p[1]);
        return p;
    }

    // p is in the canonical base tetrahedron if it's in the positive unit cube
    // and its components are in non-ascending order.
    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        return Orthotropic<TOL>::inBaseUnit(p) &&
                (p[0] + tolerance >= p[1]) &&
                (p[1] + tolerance >= p[2]);
    }

    template<typename Real>
    static bool inMeshingCell(const Point3<Real> &p) {
        return Orthotropic<TOL>::inMeshingCell(p);
    }

    // All vertices in the cubic base unit are independent.
    template<typename Real>
    static Point3<Real> independentVertexPosition(Point3<Real> p) {
        assert(inBaseUnit(p));
        return p;
    }

    // Octahedral group symmetries are the reflections of the Orthotropic class'
    // group combined with all 6 axis permutations
    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group;
        std::vector<Isometry> parentGroup = Orthotropic<TOL>::symmetryGroup();
        for (const Isometry &p : parentGroup) {
            group.push_back(p); // X Y Z

            group.push_back(p.compose(Isometry::permutation(Axis::X, Axis::Y))); // Y X Z
            group.push_back(p.compose(Isometry::permutation(Axis::X, Axis::Z))); // Z Y X
            group.push_back(p.compose(Isometry::permutation(Axis::Y, Axis::Z))); // X Z Y

            group.push_back(p.compose(Isometry::permutation(Axis::Y, Axis::Z)).compose(Isometry::permutation(Axis::X, Axis::Y))); // Y Z X
            group.push_back(p.compose(Isometry::permutation(Axis::X, Axis::Y)).compose(Isometry::permutation(Axis::Y, Axis::Z))); // Z X Y
        }
        return group;
    }
};

// Square symmetry: implemented as cubic symmetry with the unwanted symmetry
// operations (those permuting out of the z = 0 plane) filtered out.
// Symmetry group D8 x Translations
template<typename TOL>
struct Square : public Orthotropic<TOL>, SymmetryCRTP<Square<TOL>> {
    typedef TOL Tolerance;
    using Orthotropic<TOL>::representativeMeshCell;
    using TriplyPeriodic<TOL>::tolerance;

    // Disambiguate CRTP instances
    typedef SymmetryCRTP<Square<TOL>> CRTP;
    using CRTP::nodePositioner;
    using CRTP::nodeType;

    template<typename Real>
    static Point3<Real> mapToBaseUnit(Point3<Real> p) {
        return Cubic<TOL>::mapToBaseUnit(p);
    }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) {
        return Cubic<TOL>::inBaseUnit(p);
    }

    template<typename Real>
    static bool inMeshingCell(const Point3<Real> &p) {
        return Orthotropic<TOL>::inMeshingCell(p);
    }

    template<typename Real>
    static Point3<Real> independentVertexPosition(Point3<Real> p) {
        return Orthotropic<TOL>::independentVertexPosition(p);
    }

    static std::vector<Isometry> symmetryGroup() {
        std::vector<Isometry> group;
        for (const Isometry &iso : Cubic<TOL>::symmetryGroup()) {
            if (!iso.affectsAxis(Axis::Z))
                group.push_back(iso);
        }
        return group;
    }
};

// Patterns with no symmetries
template<typename TOL>
struct Null {
    using Tolerance = TOL;
    static constexpr double tolerance = double(TOL::num) / double(TOL::den);

    // We still probably want to mesh only the "period cell"
    // (e.g., for periodic tilings, we mesh one cell of the tiling at a time for better efficiency.)
    template<typename Real>
    static BBox<Point3<Real>> representativeMeshCell() { return TriplyPeriodic<TOL>::template representativeMeshCell<Real>(); }
    template<typename Real>
    static bool inMeshingCell(const Point3<Real> &p) { return TriplyPeriodic<TOL>::inBaseUnit(p); }

    template<typename Real>
    static Point3<Real> mapToBaseUnit(const Point3<Real> &p) { return p; }

    template<typename Real>
    static bool inBaseUnit(const Point3<Real> &p) { return TriplyPeriodic<TOL>::inBaseUnit(p); }

    // No symmetries
    static std::vector<Isometry> symmetryGroup() { return std::vector<Isometry>(); }
};

} // end of namespace Symmetry

#endif /* end of include guard: SYMMETRY_HH */
