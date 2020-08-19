////////////////////////////////////////////////////////////////////////////////
// BilinearMap.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Simple class for representing a bilinear map from a reference square [-1,1]²
//  to an arbitrary quadrilateral.
*/
//  Author:  Jérémie Dumas (jdumas), jeremie.dumas@ens-lyon.org
//  Created:  08/08/2018 11:00:00
////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "InflatorTypes.hh"

// By default it maps the square [-1,1]² to itself, and the jacobian is constant
// and equal to the identity.
class BilinearMap {

public:
    BilinearMap()
        : a(-1, -1, 0)
        , b( 1, -1, 0)
        , c( 1,  1, 0)
        , d(-1,  1, 0)
    { }

    template<typename T>
    BilinearMap(T pts) {
        a[0] = pts[0][0]; a[1] = pts[0][1];
        b[0] = pts[1][0]; b[1] = pts[1][1];
        c[0] = pts[2][0]; c[1] = pts[2][1];
        d[0] = pts[3][0]; d[1] = pts[3][1];
        a[2] = b[2] = c[2] = d[2] = 0.0;
    }

    template<typename Real>
    Point3<Real> apply(Real u, Real v) const {
        return (a*(u*v - u - v + 1) + b*(-u*v + u - v + 1) + c*(u*v + u + v + 1) + d*(-u*v - u + v + 1))/4;
    }

    Eigen::Matrix3d jacobian(double u, double v) const {
        //Point3d dfdu = -a/2 + b/2 + (v/2 + 1.0/2.0)*(a - b + c - d)/2;
        //Point3d dfdv = -a/2 + d/2 + (u/2 + 1.0/2.0)*(a - b + c - d)/2;
        Point3d dfdu  = a*(v - 1)/4 - b*(v - 1)/4 + c*(v + 1)/4 - d*(v + 1)/4;
        Point3d dfdv = a*(u - 1)/4 - b*(u + 1)/4 + c*(u + 1)/4 - d*(u - 1)/4;
        Eigen::Matrix3d jac;
        jac <<
            dfdu[0], dfdv[0], 0,
            dfdu[1], dfdv[1], 0,
            dfdu[2], dfdv[2], 1;
        return jac;
    }

private:
    Point3d a, b, c, d;
};
