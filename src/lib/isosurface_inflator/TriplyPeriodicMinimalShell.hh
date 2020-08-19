////////////////////////////////////////////////////////////////////////////////
// TriplyPeriodicMinimalShell.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  The extrusion of a triply-periodic minimal surface as defined in
//  [Cvijovi and Klinowski 1992] by the implicit surface function:
//      phi(x) = sum_{k = 1}^K A_k cos(2 * Pi (h_k . x) / lambda_k + P_k)
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  07/10/2018 14:14:56
////////////////////////////////////////////////////////////////////////////////
#ifndef TRIPLYPERIODICMINIMALSHELL_HH
#define TRIPLYPERIODICMINIMALSHELL_HH

#include "SignedDistance.hh"
#include <stdexcept>

struct TriplyPeriodicMinimalShell : public SignedDistanceRegion<3> {
    using Real = double;

    TriplyPeriodicMinimalShell(const std::vector<Real> &A,
                               const std::vector<Vector3d> &h,
                               const std::vector<Real> &lambda,
                               const std::vector<Real> &P,
                               const Real c) {
        setParameters(A, h, lambda, P, c);
    }

    void setParameters(const std::vector<Real> &A,
                       const std::vector<Vector3d> &h,
                       const std::vector<Real> &lambda,
                       const std::vector<Real> &P,
                       const Real c) {
        m_K = A.size();
        if ((m_K != h.size()) || (m_K != lambda.size()) || (m_K != P.size()))
            throw std::runtime_error("Parameter size mismatch");
        m_A = A;
        m_h = h;
        m_lambda = lambda;
        m_P = P;
        m_c = c;
    }

    virtual void boundingSphere(Point3D &c, double &r) const override {
        c.setZero();
        r = 3.0;
    }

    virtual const BBox<PointNd<3>> & boundingBox() const override {
        static auto box = BBox<Point3D>(Point3D(0, 0, 0), Point3D(1, 1, 1));
        return box;
    }

    virtual Real signedDistance(const Point3D &p) const override {
        Real phi = 0;
        for (size_t k = 0; k < m_K; ++k)
            phi += m_A[k] * cos(2 * M_PI * m_h[k].dot(p) / m_lambda[k] + m_P[k]);
        return std::abs(phi) - m_c;
    }

    virtual bool isInside(const Point3D &p) const override {
        return signedDistance(p) <= 0;
    }

private:
    size_t m_K;
    std::vector<Real> m_A;
    std::vector<Vector3d> m_h;
    std::vector<Real> m_lambda;
    std::vector<Real> m_P;
    Real m_c;
};

#endif /* end of include guard: TRIPLYPERIODICMINIMALSHELL_HH */
