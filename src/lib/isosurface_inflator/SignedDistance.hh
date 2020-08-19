////////////////////////////////////////////////////////////////////////////////
// SignedDistance.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Signed distance function utilities. These functions are templated to
//      support auto-differentiation via adept.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/23/2015 14:59:39
////////////////////////////////////////////////////////////////////////////////
#ifndef SIGNEDDISTANCE_HH
#define SIGNEDDISTANCE_HH

#include "InflatorTypes.hh"
#include "TriangleClosestPoint.hh"
#include "AutomaticDifferentiation.hh"
#include <vector>

namespace SignedDistance {

template<typename Real>
Real clamp(Real x, Real a, Real b) { return std::min<Real>(std::max<Real>(x, a), b); }
template<typename Real>
Real mix(Real x, Real y, Real a) { return x * (1 - a) + y * a; }

// Min of a collection of values.
template<typename Real>
Real min(const std::vector<Real> &values) {
    Real res = std::numeric_limits<Real>::max();
    for (Real v : values)
        res = std::min(res, v);
    return res;
}

// exponential smooth min of two values (k = 32);
template<typename Real>
Real exp_smin(Real a, Real b, Real k = 32)
{
    Real res = exp(-k * a) + exp(-k * b);
    return -log(res) / k;
}

// exponential smooth min one or more values (k = 32);
template<typename Real>
Real exp_smin(const std::vector<Real> &values, Real k = 32)
{
    Real res = 0;
    for (Real v : values) {
        res += exp(-k * v);
    }
    return -log(res) / k;
}

// exponential smooth min of two values (s = 1/32);
// reparametrized version: shape velocities are better behaved if we
// parametrize by k = 1/s
template<typename Real>
Real exp_smin_reparam(Real a, Real b, Real s = 1.0/32)
{
    Real res = exp(-a / s) + exp(-b / s);
    return -log(res) * s;
}

// We can rewrite the smin of two numbers in a more numerically stable way by
// decomposing into the average and difference of the numbers:
// (a + b) / 2 - s log(exp((a - b)/(2s)) + exp((b - a)/(2s)))
// Then when an overflow occurs with exp, we know that we can safely clamp to
// one or the other value.
// We can predict the overflow in either the single or double precision case:
// In single precision, expf(x) overflows when x >= log(2** 128) ~= 88.7228
// In double precision, exp (x) overflows when x >= log(2**1024) ~= 709.7827.
// For now, we assume double precision...
template<typename Real>
Real exp_smin_reparam_accurate(Real a, Real b, Real s = 1.0/32) {
    if (s == 0.0) return std::min(a, b);
    Real d_div_s = (a - b) / (2.0 * s);
    const double maxDifference = 709.7827;
    if (d_div_s >  maxDifference) return b;
    if (d_div_s < -maxDifference) return a;
    // return (a + b) / 2.0 - s * log(exp(d_div_s) + exp(-d_div_s));
    // Note: the following sidesteps an overflow issue with the derivative wrt s:
    return (a + b) / 2.0 - s * (log_cosh(d_div_s) + log(2.0));
}

// We can rewrite the smin of N numbers in a more numerically stable way by
// decomposing into the average and difference of the numbers:
// (a + b + ...) / N - s log(exp((avg - a)/(2s)) + exp((avg - b)/(2s)) + ...)
// Then when an overflow occurs with exp, we know that we can safely clamp to
// the smallest number.
// We can predict the overflow in either the single or double precision case:
// In single precision, N * expf(x) overflows when x >= log(2** 128 / N) = log(2** 128) - log(N) ~= 88.7228
// In double precision, N * exp (x) overflows when x >= log(2**1024 / N) = log(2**1024) - log(N) ~= 709.7827.
// So, if we choose the cutoff for x at 700 (double precision), we should avoid
// overflow for N = exp(9.7827...) ~= 17724.
template<typename Real>
Real exp_smin_reparam_accurate(const std::vector<Real> &values, Real s) {
    Real minVal = values[0];
    Real avg = 0;
    Real n = values.size();
    for (const Real &val : values) {
        minVal = std::min(minVal, val);
        avg += val;
    }
    avg /= n;

    if (s == 0.0) return minVal;

    // TODO: make this computation more stable for autodiff partial derivatives
    // with small s (like the log_cosh version for 2 value). For now, we use a
    // more conservative threshold: const double maxDifference = 700;
    const double maxDifference = 600;
    if ((avg - minVal) > maxDifference * s) return minVal;

    Real k = 1.0 / s;
    Real res = 0;
    for (auto val : values)
        res += exp((avg - val) * k);
    return avg - s * log(res);
}

// exponential smooth min one or more values (s = 1/32);
template<typename Real>
Real exp_smin_reparam(const std::vector<Real> &values, Real s = 1.0/32) {
    Real res = 0;
    Real k = 1.0 / s;
    for (Real v : values)
        res += exp(-k * v);
    return -log(res) / k;
}

// Minimum max(a, b) - min(a, b) such that
// min(a, b) - exp_smin(a, b) < tol
// (Note: exp_smin always under-estimates, so this is a bound on the absolute
// error of the min approximation).
template<typename Real>
Real exp_smin_radius(Real k = 20, Real tol = 1e-3) {
    return -log(exp(k * tol) - 1.0) / k;
}

template<typename Real>
Real poly_smin(Real a, Real b, Real k = 0.2) {
    Real h = clamp(0.5 + 0.5 * (b - a) / k, Real(0.0), Real(1.0));
    return mix(b, a, h) - k * h * (1.0 - h);
}

// power smooth min (k = 8);
template<typename Real>
Real pow_smin( Real a, Real b, Real k = 5 )
{
    a = pow( a, k ); b = pow( b, k );
    return pow( (a*b)/(a+b), 1.0/k );
}

template<typename Real>
Real smax(Real a, Real b, Real k = 0.2) {
    return -poly_smin(-a, -b, k);
}

////////////////////////////////////////////////////////////////////////////////
// Primitives
////////////////////////////////////////////////////////////////////////////////
namespace Primitives {

template<typename Real>
struct Sphere {
    Sphere() { }
    Sphere(const Point3<Real> &c, Real r) { set(c, r); }
    void set(const Point3<Real> &c, Real r) { m_c = c; m_r = r; }

    Real signedDistance(const Point3<Real> &p) const {
        return sqrt((p - m_c).squaredNorm()) - m_r;
    }
private:
    Point3<Real> m_c;
    Real m_r;
};

template<typename Real>
struct Box {
    Box() { }
    Box(const BBox<Point3<Real>> &box) : m_box(box) { }
    void set(const BBox<Point3<Real>> &box) { m_box = box; }
    Real signedDistance(const Point3<Real> &p) const {
        // Box is intersection of 6 halfspaces
        Vector3<Real> d = (m_box.minCorner - p).cwiseMax(p - m_box.maxCorner);
        return std::min(d.maxCoeff(), Real(0)) + d.cwiseMax(Vector3<Real>::Zero()).norm();
    }
    bool isInside(Real p0, Real p1, Real p2) const {
        return (p0 >= m_box.minCorner[0]) &&
               (p1 >= m_box.minCorner[1]) &&
               (p2 >= m_box.minCorner[2]) &&
               (p0 <= m_box.maxCorner[0]) &&
               (p1 <= m_box.maxCorner[1]) &&
               (p2 <= m_box.maxCorner[2]);
    }
private:
    BBox<Point3<Real>> m_box;
};

template<typename Real>
struct ConicalFrustum {
    ConicalFrustum() { }
    ConicalFrustum(const Point3<Real> &a, const Point3<Real> &b, Real ra, Real rb) { set(a, b, ra, rb); }

    void set(const Point3<Real> &a, const Point3<Real> &b, Real ra, Real rb) {
        m_a = a;
        m_ra = ra; m_rb = rb;
        m_axis = b - a;
        m_axisLenSq = m_axis.squaredNorm();
        m_axisLen = sqrt(m_axisLenSq);
    }

    // WARNING: NOT A TRUE SIGNED DISTANCE FUNCTION
    Real signedDistance(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_a);
        Real frac_along = v.dot(m_axis) / m_axisLenSq;
        Real dist_axial = m_axisLen * (fabs(frac_along - 0.5) - 0.5);

        Vector3<Real> v_perp = v - frac_along * m_axis;
        Real r = m_ra + clamp(frac_along, Real(0.0), Real(1.0))  * (m_rb - m_ra);
        Real dist_perp = sqrt(v_perp.squaredNorm()) - r;

        Vector2<Real> posDists(std::max(dist_perp, Real(0.0)), std::max(dist_axial, Real(0.0)));
        return sqrt(posDists.squaredNorm()) + std::min(std::max(dist_perp, dist_axial), Real(0.0));
    }

private:
    Point3<Real> m_a;
    Vector3<Real> m_axis;
    Real m_axisLen, m_axisLenSq, m_ra, m_rb;
};

template<typename Real>
struct Cylinder {
    Cylinder() { }
    Cylinder(const Point3<Real> &a, const Point3<Real> &b, Real r) { set(a, b, r); }

    void set(const Point3<Real> &a, const Point3<Real> &b, Real r) {
        m_a = a;
        m_axis = b - a;
        m_r = r;
        m_axisLenSq = m_axis.squaredNorm();
        m_axisLen = sqrt(m_axisLenSq);
    }

    Real signedDistance(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_a);
        Real frac_along = v.dot(m_axis) / m_axisLenSq;
        Real dist_axial = m_axisLen * (std::abs(frac_along - 0.5) - 0.5);

        Vector3<Real> v_perp = v - frac_along * m_axis;
        Real dist_perp = sqrt(v_perp.squaredNorm()) - m_r;

        Vector2<Real> posDists(std::max(dist_perp, 0.0), std::max(dist_axial, 0.0));
        return sqrt(posDists.squaredNorm()) + std::min(std::max(dist_perp, dist_axial), 0.0);
    }

private:
    Point3<Real> m_a;
    Vector3<Real> m_axis;
    Real m_axisLen, m_axisLenSq, m_r;
};

// InflatedEdge "Primitive": edge geometry for an inflated wire mesh.
// The edge geometry consists of a "cylinder" with sphere endpoints. To support
// per-vertex thickness, the spheres can be of different radii, in which case
// the cylinder part is really a conical frustum (linearly interpolated
// radius). The conical frustum is created so that it joins smoothly with the
// sphere endpoints (tangent at the intersection--continuous normal).
template<typename Real>
class InflatedEdge {
public:
    InflatedEdge() { }
    InflatedEdge(const Point3<Real> &p1, const Point3<Real> &p2, Real r1, Real r2)
    { set(p1, p2, r1, r2); }

    void set(const Point3<Real> &p1, const Point3<Real> &p2, const Real r1, const Real r2) {
        m_c1 = p1, m_c2 = p2;
        m_r1 = r1, m_r2 = r2;

        if (p2 == p1) throw std::runtime_error("Degenerate inflated edge");
        m_axisUnit = p2 - p1;
        m_axisLength = sqrt(m_axisUnit.squaredNorm());
        m_axisUnit /= m_axisLength;

        // Note: this actually computes a clockwise angle... switch to ccw
        // later? (But all class methods use the clockwise angle properly.)
        m_sinTheta = (r1 - r2) / m_axisLength;
        m_sinTheta = std::max(Real(-1.0), std::min(Real(1.0), m_sinTheta)); // Avoid catastrophic failure in degenerate configurations
        m_theta = asin(m_sinTheta);
        m_cosTheta = sqrt(1 - m_sinTheta * m_sinTheta);

        m_edgeLength = m_axisLength * m_cosTheta;
    }

    // Additional real type to support automatic differentiation wrt. p
    // even when the class's real type doesn't support autodiff.
    template<typename Real2, bool DebugOutput = false>
    Real2 signedDistance(const Point3<Real2> &p) const {
        // Note: we must cast the vector-type member variables to Real2 types
        // for the automatic differentation case. There's a chance this will
        // make the Real2 == Real case less efficient; we can add some
        // SFINAE voodoo to optimize this if it's an issue.

        Vector3<Real2> v(p - m_c1.template cast<Real2>());
        Real2 v_normSq = v.squaredNorm();
        Real2 v_parallelComponent = v.dot(m_axisUnit.template cast<Real2>());
        // Max is to avoid numerical issues... (should be fine since we never
        // need to auto-diff where v_perpComponent = 0, i.e. deep inside
        // object).
        // Note: With the sphere hull joint blending, we *do* sometimes care
        // about derivatives at the midline; this max should set the derivative
        // to 0 rather than NaN, which should be fine: the blending modulation
        // function should be design so that the blending derivative is zero at
        // the hull medial axis. We just need to keep NaNs from propagating.
        Real2 v_perpComponent = sqrt(std::max(Real2(v_normSq - v_parallelComponent * v_parallelComponent), Real2(1e-16)));

        // Rotate so that conical frustum surface is horizontal
        Real2 x = m_cosTheta * v_parallelComponent - m_sinTheta * v_perpComponent;

        if (DebugOutput) {
            std::cerr << "InflatedEdge::signedDistance derivative info:" << std::endl;
            std::cerr << "v_normSq ( " << v_normSq << "):"; reportDerivatives(std::cerr, v_normSq); std::cerr << std::endl;
            std::cerr << "v_parallelComponent ( " << v_parallelComponent << "):"; reportDerivatives(std::cerr, v_parallelComponent); std::cerr << std::endl;
            std::cerr << "v_perpComponent ( " << v_perpComponent << "):"; reportDerivatives(std::cerr, v_perpComponent); std::cerr << std::endl;
            std::cerr << "x ( " << x << "):"; reportDerivatives(std::cerr, x); std::cerr << std::endl;
            std::cerr << std::endl;
        }

        // Closest surface is sphere 1
        if (x < 0) {
            Real2 result = sqrt(v_normSq); // work around make_coherent bug in Eigen autodiff
            return result - m_r1;
        }

        // Closest surface is the conical frustum part (the closest edge of
        // which is horizontal in rotated (x, y) coordinates).
        if (x < m_edgeLength) {
            Real2 y = m_sinTheta * v_parallelComponent + m_cosTheta * v_perpComponent;
            return y - m_r1;
        }

        // Closest surface is sphere 2
        return (p - m_c2.template cast<Real2>()).norm() - m_r2;
    }

    Point3<Real> closestMedialAxisPoint(const Point3<Real> &p) const {
        Real parallelDist = (p - m_c1).dot(m_axisUnit);
        // Clamp 'parallelDist' to stay in the line segment (0, m_axisLength)
        return m_c1 + clamp(parallelDist, Real(0), m_axisLength) * m_axisUnit;
    }

    // Closest point to "p" lying on the intersection of the conical frustum and
    // the (hemi)sphere endcap "eci" This is useful for determining the amount
    // of joint smoothing. (eci is 0 or 1)
    // TODO: autodiff-compatible types
    Point3<Real> closestFrustumBorderPoint(const Point3<Real> &p, int eci) const {
        Vector3<Real> v(p - m_c1);
        Real v_parallelComponent = v.dot(m_axisUnit);
        Vector3<Real> v_perp = v - v_parallelComponent * m_axisUnit;
        Real v_perpComponent = v_perp.norm();

        // Clamp to specified frustum boundary
        Real x = (eci == 0) ? 0.0 : m_edgeLength;
        Real y = m_r1;

        // Rotate back to get the result.
        Point3<Real> result;
        result = m_axisUnit *  ( m_cosTheta * x + m_sinTheta * y)
               + v_perp     * ((-m_sinTheta * x + m_cosTheta * y) / v_perpComponent); // (v_perp not normalized)
        return result;
    }

    Point3<Real> closestPoint(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_c1);
        Real v_parallelComponent = v.dot(m_axisUnit);
        Vector3<Real> v_perp = v - v_parallelComponent * m_axisUnit;
        Real v_perpComponent = v_perp.norm();

        // Rotate so that conical frustum surface is horizontal
        Real x = m_cosTheta * v_parallelComponent - m_sinTheta * v_perpComponent;

        // Closest point on sphere 1
        if (x < 0) { return m_c1 + (m_r1 / v.norm()) * v; }
        // Closest point on conical frustum
        if (x < m_edgeLength) {
            // Rotate back and add to c1 to get the result.
            return m_c1 +
                   m_axisUnit *  ( m_cosTheta * x + m_sinTheta * m_r1)
                 + v_perp     * ((-m_sinTheta * x + m_cosTheta * m_r1) / v_perpComponent); // (v_perp not normalized)
        }

        // Closest point on sphere 2
        Vector3<Real> dir = (p - m_c2);
        return m_c2 + (m_r2 / dir.norm()) * dir;
    }

    // Get the normal for the surface point closest to p.
    Vector3<Real> normal(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_c1);
        Real v_parallelComponent = v.dot(m_axisUnit);
        Vector3<Real> v_perp = v - v_parallelComponent * m_axisUnit;
        Real v_perpComponent = v_perp.norm();

        // Rotate so that conical frustum surface is horizontal
        Real x = m_cosTheta * v_parallelComponent - m_sinTheta * v_perpComponent;

        // Closest point on sphere 1; normal is just v normalized
        if (x < 0) return v / v.norm();

        // Closest surface is the conical frustum part; the normal is (0, 1)
        // in the rotated coordinates.
        if (x < m_edgeLength) {
            // Rotate normal back: R^T (0 1)
            return m_axisUnit * m_sinTheta
                 + v_perp * (m_cosTheta / v_perpComponent); // (v_perp not normalized)
        }

        // Closest point on sphere 2
        Vector3<Real> result = (p - m_c2);
        result /= result.norm();
        return result;
    }

    // Tangent vector along the conical frustum at the frustum surface point
    // closest to p
    Vector3<Real> frustumAxialTangent(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_c1);
        Real v_parallelComponent = v.dot(m_axisUnit);
        Vector3<Real> v_perp = v - v_parallelComponent * m_axisUnit;
        Real v_perpComponent = v_perp.norm();

        // The tangent is (1, 0) in the rotated coordinates.
        // Rotate tangent back: R^T (1 0)
        return m_axisUnit * m_cosTheta - v_perp * (m_sinTheta / v_perpComponent); // (v_perp not normalized)
    }

    // Closest surface frustum point to p
    Point3<Real> closestFrustumPoint(const Point3<Real> &p) const {
        Vector3<Real> v(p - m_c1);
        Real v_parallelComponent = v.dot(m_axisUnit);
        Vector3<Real> v_perp = v - v_parallelComponent * m_axisUnit;
        Real v_perpComponent = v_perp.norm();

        // Clamp to frustum
        Real x = std::max(std::min<Real>(x, m_edgeLength), 1.0);
        Real y = m_r1;

        // Rotate back to get the result.
        Point3<Real> result;
        result = m_axisUnit *  ( m_cosTheta * x + m_sinTheta * y)
               + v_perp     * ((-m_sinTheta * x + m_cosTheta * y) / v_perpComponent); // (v_perp not normalized)
        return result;
    }

    // Additional real type to support automatic differentiation wrt. p
    // even when the class's real type doesn't support autodiff.
    // Type: 0 means m_c1 is closest, 1 means m_c2 is closest, 2 = edge is
    // closest.
    template<typename Real2>
    Real2 signedDistanceAndClosestType(const Point3<Real2> &p, size_t &type) const {
        // Note: we must cast the vector-type member variables to Real2 types
        // for the automatic differentation case. There's a chance this will
        // make the Real2 == Real case less efficient; we can add some
        // SFINAE vodoo to optimize this if it's an issue.

        Vector3<Real2> v(p - m_c1.template cast<Real2>());
        Real2 v_normSq = v.squaredNorm();
        Real2 v_parallelComponent = v.dot(m_axisUnit.template cast<Real2>());
        // Max is to avoid numerical issues... (should be fine since we never
        // need to auto-diff where v_perpComponent = 0, i.e. deep inside
        // object).
        Real2 v_perpComponent = sqrt(std::max(Real2(v_normSq - v_parallelComponent * v_parallelComponent), Real2(0.0)));

        // Rotate so that conical frustum surface is horizontal
        Real2 x = m_cosTheta * v_parallelComponent - m_sinTheta * v_perpComponent;

        // Closest surface is sphere 1
        if (x < 0) {
            type = 0;
            return sqrt(v_normSq) - m_r1;
        }

        // Closest surface is the conical frustum part (the closest edge of
        // which is horizontal in rotated (x, y) coordinates).
        if (x < m_edgeLength) {
            type = 2;
            Real2 y = m_sinTheta * v_parallelComponent + m_cosTheta * v_perpComponent;
            return y - m_r1;
        }

        // Closest surface is sphere 2
        type = 1;
        return sqrt((p - m_c2.template cast<Real2>()).squaredNorm()) - m_r2;
    }

    // Signs are based on whether the overlap with the enpoint sphere decreases
    // (positive) or increases (negative)
    Real angleAtP1() const { return  m_theta; }
    Real angleAtP2() const { return -m_theta; }

    Point3<Real> c1() const { return m_c1; }
    Point3<Real> c2() const { return m_c2; }
    Real         r1() const { return m_r1; }
    Real         r2() const { return m_r2; }

private:
    Point3<Real> m_c1, m_c2;
    Real m_r1, m_r2;

    Real m_theta, m_sinTheta, m_cosTheta, m_axisLength, m_edgeLength;
    Vector3<Real> m_axisUnit;
};

// Convex hull of three spheres, **except in degenerate cases**, where it always
// becomes the convex hull of spheres 2 and 3.
// This ensures that the correct blending region is always computed for pair of
// edges ((1, 2), (2, 3)): if the true convex
// hull actually consists of the convex hulls of spheres (1, 2) and (2, 3), no
// concave edges are formed, meaning no blending is necessary. Thus, the region
// excluded by instead computing region hull((2, 3)) is irrelevant.
template<typename Real>
struct InflatedTriangle {
    InflatedTriangle() { }

    InflatedTriangle(const Point3<Real> &p1,
                     const Point3<Real> &p2,
                     const Point3<Real> &p3,
                     Real r1, Real r2, Real r3)
    { set(p1, p2, p3, r1, r2, r3); }

    void set(const Point3<Real> &p1,
             const Point3<Real> &p2,
             const Point3<Real> &p3,
             Real r1, Real r2, Real r3) {
        m_degenerate = false;

        m_c.col(0) = p1, m_c.col(1) = p2, m_c.col(2) = p3;
        m_r[0] = r1, m_r[1] = r2, m_r[2] = r3;

        Vector3<Real> e1 = p2 - p1,
                      e2 = p3 - p1;
        Real l1 = e1.norm(),
             l2 = e2.norm();
        e1 /= l1;
        e2 /= l2;

        m_midplaneNormal = e1.cross(e2);
        m_midplaneNormal /= m_midplaneNormal.norm();
        Vector3<Real> e1perp = m_midplaneNormal.cross(e1);
        // TODO: There are degenerate cases where this will fail: detect them!
        assert(std::abs(e1perp.norm() - 1.0) < 1e-10);

        // Determine the (positive) tangent plane for the three spheres.
        // We compute the tangent plane's normal in terms of its components
        // along e1, e1^perp, and midplane normal n_mp:
        //  n = alpha e1 + beta e1perp + gamma n_mp
        // Here e1, e2, e1perp, and n_mp are unit vectors
        Real alpha = (r1 - r2) / l1;
        Real beta = ((r1 - r3) / l2 - alpha * e2.dot(e1)) / e2.dot(e1perp);
        Real aSqbSq = alpha * alpha + beta * beta;
        if (aSqbSq <= 1) {
            Real gamma = sqrt(1 - aSqbSq);
            Vector3<Real> n = alpha * e1 + beta * e1perp + gamma * m_midplaneNormal;
            m_hullTangentTriangle.template set<Real>(p1 + r1 * n, p2 + r2 * n, p3 + r3 * n);

            // p3
            // ^--
            // |  --
            // |    --
            // p1 -----> p2
            m_prismBorderNormals[0] = (p3 - p2).cross(n);
            m_prismBorderNormals[1] = (p1 - p3).cross(n);
            m_prismBorderNormals[2] = (p2 - p1).cross(n);

            m_prismBorderNormals[0] /= m_prismBorderNormals[0].norm();
            m_prismBorderNormals[1] /= m_prismBorderNormals[1].norm();
            m_prismBorderNormals[2] /= m_prismBorderNormals[2].norm();

            // edge[i] is opposite sphere i (p_{i + 1})
            m_edge[0].set(p2, p3, r2, r3);
            m_edge[1].set(p3, p1, r3, r1);
            m_edge[2].set(p1, p2, r1, r2);
        }
        else {
            m_degenerate = true;
            m_edgeForDegenerateCase.set(p2, p3, r2, r3);
        }
    }

    template<typename Real2>
    Real2 signedDistance(const Point3<Real2> &p) const {
        if (m_degenerate)
            return m_edgeForDegenerateCase.signedDistance(p);

        // Reflect point to the positive side of the midplane
        Point3<Real2> p_pos = p - std::min<Real2>(2 * (p - m_c.col(0)).dot(m_midplaneNormal), 0) * m_midplaneNormal;

        // Find barycentric coordinates of closest point on convex hull tangent
        // triangle.
        auto lambda = m_hullTangentTriangle.closestBaryCoords(p_pos);

        int numZero = 0;
        int zeroCoords[3];
        if (lambda[0] == 0.0) { zeroCoords[numZero++] = 0; }
        if (lambda[1] == 0.0) { zeroCoords[numZero++] = 1; }
        if (lambda[2] == 0.0) { zeroCoords[numZero++] = 2; }
        assert(numZero != 3);

        // If exactly one barycentric coordinate is zero (i.e. we're on the hull
        // tangent triangle's edge), the closest surface point to p lies on the
        // convex hull of the two endpoint spheres; delegate to InflatedEdge.
        if (numZero == 1) {
            return m_edge[zeroCoords[0]].signedDistance(p);
        }
        if (numZero == 2) {
            // Closest point is on sphere corresponding to vertex with
            // nonzero (actually, 1) barycentric coordinate
            int nonzeroCoord = 3 - zeroCoords[0] - zeroCoords[1];

            // If the closest point on the hull tangent triangle is a vertex, we
            // we need an additional check to determine which of the incident
            // edges we are closest to. We do this by checking which side of the
            // plane spanned by m_tangentTriNormal and one of the incident
            // center edges
            // TODO: check if this heuristic computes the correct signed
            // distance.
            Vector3<Real2> v = p_pos - m_c.col(nonzeroCoord);
            if (v.dot(m_prismBorderNormals[zeroCoords[0]]) >= v.dot(m_prismBorderNormals[zeroCoords[1]]))
                return m_edge[zeroCoords[0]].signedDistance(p);
            else
                return m_edge[zeroCoords[1]].signedDistance(p);
        }

        // Now the closest surface point lies on a sphere whose center and radii
        // are given by interpolating the three spheres' centers and radii with
        // these barycentric coordinates.

        Point3<Real2> c = m_c * lambda;
        Real2 r = m_r.dot(lambda);

        return (p_pos - c).norm() - r;
    }

    // For debugging purposes
    const TriangleClosestPoint<Real> &hullTangentTriangle() const { return m_hullTangentTriangle; }

private:
    bool m_degenerate;
    // Normal vectors for plane containing sphere centers.
    // Tangent vector 1 is along edge m_c2 - m_c1
    Vector3<Real> m_midplaneNormal;

    InflatedEdge<Real> m_edgeForDegenerateCase;

    // Data structure representing the triangle tangent to the three spheres on
    // the positive side of the midplane and implementing fast nearest point
    // computations.
    TriangleClosestPoint<Real> m_hullTangentTriangle;
    // Plane normals for the side faces of the prism-like region bounded above
    // and below by the hull triangles. Indexed like m_edge below.
    Vector3<Real> m_prismBorderNormals[3];

    // When the closest point lies on one of the edges' inflated geometry, we
    // must delegate to the InflatedEdge class; the linearly interpolated sphere
    // we use to get the signed distance to the hull tangent triangle does not
    // compute the correct distance in this case.
    // Note: m_edge[i] is opposite sphere i
    InflatedEdge<Real> m_edge[3];

    // Sphere center positions (position i in column i).
    Eigen::Matrix<Real, 3, 3> m_c;
    // Sphere radii
    Vector3<Real> m_r;
};

} // End of namespace: Primitives

}

#endif /* end of include guard: SIGNEDDISTANCE_HH */
