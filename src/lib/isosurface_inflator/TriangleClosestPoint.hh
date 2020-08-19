#ifndef TRIANGLECLOSESTPOINT_HH
#define TRIANGLECLOSESTPOINT_HH

#include <Eigen/Dense>

template<typename Real>
class TriangleClosestPoint {
public:
    using Vec = Eigen::Matrix<Real, 3, 1>;

    TriangleClosestPoint() { }

    template<typename Real2>
    TriangleClosestPoint(const Eigen::Matrix<Real2, 3, 1> &p0,
                         const Eigen::Matrix<Real2, 3, 1> &p1,
                         const Eigen::Matrix<Real2, 3, 1> &p2) { set(p0, p1, p2); }

    template<typename Real2>
    void set(const Eigen::Matrix<Real2, 3, 1> &p0,
             const Eigen::Matrix<Real2, 3, 1> &p1,
             const Eigen::Matrix<Real2, 3, 1> &p2) {
        m_p0 = p0; m_p1 = p1; m_p2 = p2;

        Vec e0 = p2 - p1,
            e1 = p0 - p2,
            e2 = p1 - p0;
        m_normal = e0.cross(e1);
        Real doubleA = m_normal.norm();
        m_normal /= doubleA;

        // grad lambda_i = e_i^perp / 2A
        m_jacobianLambda.row(0) = m_normal.cross(e0);
        m_jacobianLambda.row(1) = m_normal.cross(e1);
        m_jacobianLambda.row(2) = m_normal.cross(e2);
        m_jacobianLambda *= 1.0 / doubleA;

        // Compute M_ij:
        //      Change in lambda_i when moving 1 unit in lambda_j direction.
        m_M(0, 1) = m_M(1, 0) = e0.dot(e1);
        m_M(0, 2) = m_M(2, 0) = e0.dot(e2);
        m_M(1, 2) = m_M(2, 1) = e1.dot(e2);
        m_M.col(0) /= e0.squaredNorm();
        m_M.col(1) /= e1.squaredNorm();
        m_M.col(2) /= e2.squaredNorm();
        m_M(0, 0) = m_M(1, 1) = m_M(2, 2) = 1.0;
    }

    template<typename Real2>
    Eigen::Matrix<Real2, 3, 1> pointAtBarycoords(const Eigen::Matrix<Real2, 3, 1> &lambda) const {
        return m_p0.template cast<Real2>() * lambda[0] +
               m_p1.template cast<Real2>() * lambda[1] +
               m_p2.template cast<Real2>() * lambda[2];
    }

    // Closest triangle point to p
    template<typename Real2>
    Eigen::Matrix<Real2, 3, 1> operator()(const Eigen::Matrix<Real2, 3, 1> &p) const {
        return pointAtBarycoords(closestBaryCoords(p));
    }

    // Alias for operator(): find closest point on the triangle to p
    template<typename Real2>
    Eigen::Matrix<Real2, 3, 1> closestPoint(const Eigen::Matrix<Real2, 3, 1> &p) const {
        return (*this)(p);
    }

    // Get barycentric coordinates of point in triangle's plane closest to p.
    // Note: may be outside triangle; use closestBaryCoords to snap to triangle.
    template<typename Real2>
    Eigen::Matrix<Real2, 3, 1>
    baryCoords(const Eigen::Matrix<Real2, 3, 1> &p) const {
        Eigen::Matrix<Real2, 3, 1> lambda = m_jacobianLambda.template cast<Real2>() * (p - m_p0.template cast<Real2>());
        lambda[0] += 1.0;
        return lambda;
    }

    template<typename Real2>
    Eigen::Matrix<Real2, 3, 1>
	closestInternalBarycoordsToBaryCoords(Eigen::Matrix<Real2, 3, 1> lambda) const {
        using BaryC = Eigen::Matrix<Real2, 3, 1>;
        bool seenNegative = false;
        int prevNegIdx = -1;

        for (size_t i = 0; i < 3; ++i) {
            if (lambda[i] >= 0) continue;

            int j = (i == 2) ? 0 : i + 1;
            int k = (j == 2) ? 0 : j + 1;

            BaryC lambdaSnap = lambda - m_M.col(i).template cast<Real2>() * lambda[i];

            if (lambdaSnap[j] >= 0) {
                if (lambdaSnap[k] >= 0) return lambdaSnap; // snapped to closest pt in triangle
                else {
                    if (lambda[k] >= 0) { // initially positive but went negative
                        lambda[i] = 0; lambda[k] = 0;
                        lambda[j] = 1.0;
                        return lambda;
                    }
                }
            }
            else if (lambdaSnap[k] >= 0) {
                // we know lambdaSnap[j] < 0
                if (lambda[j] >= 0) { // initially positive but went negative
                    lambda[i] = 0; lambda[j] = 0;
                    lambda[k] = 1.0;
                    return lambda;
                }
            }

            if (seenNegative) {
                // There were two negative barycentric coordinates and both
                // snaps were unsuccessful. This means the closest point is the
                // triangle vertex whose corresponding barycentric coordinate is
                // positive.
                lambda.fill(1.0);
                lambda[i] = 0.0; lambda[prevNegIdx] = 0.0;
                return lambda;
            }
            seenNegative = true;
            prevNegIdx = i;
        }

        return lambda;
	}

    // Get barycentric coordinates of closest triangle point to p
    template<typename Real2>
    Eigen::Matrix<Real2, 3, 1>
    closestBaryCoords(const Eigen::Matrix<Real2, 3, 1> &p) const {
        return closestInternalBarycoordsToBaryCoords(baryCoords(p));
    }

    const Eigen::Matrix<Real, 3, 1> &p0() const { return m_p0; }
    const Eigen::Matrix<Real, 3, 1> &p1() const { return m_p1; }
    const Eigen::Matrix<Real, 3, 1> &p2() const { return m_p2; }

    const Eigen::Matrix<Real, 3, 1> &normal() const { return m_normal; }

private:
    Eigen::Matrix<Real, 3, 1> m_p0, m_p1, m_p2;
    Eigen::Matrix<Real, 3, 3> m_jacobianLambda;
    Eigen::Matrix<Real, 3, 3> m_M;
    Eigen::Matrix<Real, 3, 1> m_normal;
};

#endif /* end of include guard: TRIANGLECLOSESTPOINT_HH */
