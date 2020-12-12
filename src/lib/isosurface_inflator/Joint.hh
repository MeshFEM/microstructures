////////////////////////////////////////////////////////////////////////////////
// Joint.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Represents and computes signed distance to the smoothed joint geometry
//      at a particular vertex.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/08/2016 23:55:55
////////////////////////////////////////////////////////////////////////////////
#ifndef JOINT_HH
#define JOINT_HH

#include "SphereConvexHull.hh"
#include "AutomaticDifferentiation.hh"
#include <MeshFEM/Future.hh>
#include <stdexcept>

enum class JointBlendMode { FULL, HULL, HULL_HALF_EDGE };
enum class JointBlendFunction { EXPONENTIAL };
namespace SD = SignedDistance;

template<typename Real>
class Joint {
public:
    Joint(const std::vector<Point3<Real>> &centers,
          const std::vector<Real>         &radii,
          Real blendingAmt, JointBlendMode mode,
          JointBlendFunction blendFunction,
          const std::vector<Real> &blendingPolyCoeffs)
    { setParameters(centers, radii, blendingAmt, mode, blendFunction, blendingPolyCoeffs); }

    Joint(const std::vector<Point3<Real>> &centers,
          const std::vector<Real>         &radii,
          Real blendingAmt, JointBlendMode mode,
          JointBlendFunction blendFunction = JointBlendFunction::EXPONENTIAL)
    {
        std::vector<Real> emptyBlendingPolyCoeffs;
        setParameters(centers, radii, blendingAmt, mode, blendFunction, emptyBlendingPolyCoeffs);
    }

    // Assumes the first (center, radius) specifies the joint sphere
    void setParameters(std::vector<Point3<Real>> centers, // modified inside
                       std::vector<Real>         radii,   // modified inside
                       Real blendingAmt,
                       JointBlendMode blendMode,
                       JointBlendFunction blendFunction,
                       std::vector<Real> blendingPolyCoeffs) {
        m_mode = blendMode;
        m_blendFunction = blendFunction;
        assert(centers.size() == radii.size());
        if (centers.size() < 3) {
            {
                std::cout << "Joint constructor called on centers, radii:" << std::endl;
                for (const auto &pt : centers) {
                    std::cout << "{"
                        << pt[0] << ", "
                        << pt[1] << ", "
                        << pt[2] << "}, ";
                }
                std::cout << std::endl;
                std::cout << "{";
                for (const Real &r : radii)
                    std::cout << r << ", ";
                std::cout << "}" << std::endl;
            }
            throw std::runtime_error("Joint must comprise at least 3 spheres");
        }
        m_minRadius = *std::min_element(radii.begin(), radii.end());
        m_r1 = radii[0];
        m_c1 = centers[0];
        if (m_mode == JointBlendMode::HULL_HALF_EDGE) {
            // In HULL_HALF_EDGE mode, shrink the hull region to extend only
            // halfway through each edge so that neighboring joint blending
            // regions do not overlap. This should prevent creases from forming
            // when joint geometries are unioned together.
            for (size_t i = 1; i < centers.size(); ++i) {
                centers[i] += centers[0];
                radii[i]   += radii[0];
                centers[i] *= 0.5;
                radii[i]   *= 0.5;
            }
        }
        m_blendingAmt = blendingAmt;
        if (m_mode != JointBlendMode::FULL) {
            try {
                m_blendingHull = Future::make_unique<SD::Primitives::SphereConvexHull<Real>>(centers, radii);
            }
            catch (std::exception &e) {
                // Don't let robustness issues in the sphere convex hull crash
                // optimization. Simply disable blending in those cases.
                std::cerr << "Caught error while constructing joint hull: " << e.what() << std::endl;
                std::cerr << "Disabling joint blending" << std::endl;
                m_blendingHull = nullptr;
            }
        }

        m_blendingPolyCoeffs = blendingPolyCoeffs;
    }

    template<typename Real2, bool DebugOutput = false>
    Real2 smoothingAmt(const Point3<Real2> &p) const {
        Real2 result = m_blendingAmt;
        if (m_mode == JointBlendMode::FULL) { return m_blendingAmt; }
        if (!m_blendingHull) { return m_blendingAmt; }

        Real2 modulation = 1.0;
        try {
            Real2 hullDist = m_blendingHull->template signedDistance<Real2, DebugOutput>(p);
            // Real2 z = 1.0 + (hullDist / m_r1); // from 0 at "center" to 1 at outside
            // Note: hullDist is non-differentiable at the blending hull medial
            // axis, so we must scale the modulation region so that essentially no
            // modulation is done at the medial axis.
            Real2 z = 1.0 + (hullDist / m_minRadius); // from 0 at "center" to 1 at outside
            z = std::max<Real2>(z, 0.0);
            // Real2 modulation = 1.0 - tanh(pow(z, 8.0));
            z *= 1.025;
            modulation = 1.0 - tanh(pow(z, 10.0));

            result = modulation*result;
            if (DebugOutput || hasInvalidDerivatives(result)) {
                std::cerr << "smoothingAmt derivatives:" << std::endl;
                std::cerr << "     hullDist (" <<      hullDist << "):"; reportDerivatives(std::cerr,      hullDist); std::cerr << std::endl;
                std::cerr << "            z (" <<             z << "):"; reportDerivatives(std::cerr,             z); std::cerr << std::endl;
                std::cerr << "   modulation (" <<    modulation << "):"; reportDerivatives(std::cerr,    modulation); std::cerr << std::endl;
                std::cerr << "m_blendingAmt (" << m_blendingAmt << "):"; reportDerivatives(std::cerr, m_blendingAmt); std::cerr << std::endl;
                std::cerr << std::endl << std::endl;

                std::cerr << "m_r1, m_minRadius values: " << m_r1 << ", " << m_minRadius << std::endl;
            }
        }
        catch (std::exception &e) {
            std::cerr << "Caught error while evaluating joint hull: " << e.what() << std::endl;
            std::cerr << "Disabling joint blendig for eval pt" << std::endl;
        }

        // Original: p = 8, z *= 1.0;
        // Try p = 10, z *= 1.03 ==> 1/10th smoothing at hull, falling quickly
        // to zero

        return result;
    }

    Real blendParam() const { return m_blendingAmt; }

    template<typename Real2>
    std::vector<Real2> smoothingPolyCoeffs() {
        std::vector<Real2> result(m_blendingPolyCoeffs.size());

        for (size_t i = 0; i < m_blendingPolyCoeffs.size(); i++) {
            result[i] = m_blendingPolyCoeffs[i];
        }

        return result;
    }

    Real         r1() const { return m_r1; }
    Point3<Real> c1() const { return m_c1; }

    SD::Primitives::SphereConvexHull<Real> blendingHull() const { assert(m_blendingHull); return *m_blendingHull; }

    JointBlendFunction m_blendFunction;

private:
    std::unique_ptr<SD::Primitives::SphereConvexHull<Real>> m_blendingHull;
    Real m_blendingAmt;
    std::vector<Real> m_blendingPolyCoeffs;
    // Center and radius of the joint sphere.
    Point3<Real> m_c1;
    Real         m_r1;
    Real         m_minRadius;
    JointBlendMode m_mode;
};

#endif /* end of include guard: JOINT_HH */
