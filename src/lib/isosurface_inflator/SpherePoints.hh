////////////////////////////////////////////////////////////////////////////////
// SpherePoints.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Generates points approximately evenly distributed on a unit sphere at
//      the origin using a Fibonacci/golden sector spiral. Simplified from:
//      http://stackoverflow.com/a/26127012/122710
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/01/2016 02:58:02
////////////////////////////////////////////////////////////////////////////////
#ifndef SPHEREPOINTS_HH
#define SPHEREPOINTS_HH
#include <cmath>

// Append numPts points sampled on a unit sphere to "pts"
template<class OutputCollection, typename Real = double>
void generateSpherePoints(size_t numPts, OutputCollection &pts) {
    // Spiral stretches over [-1, 1] in the z axis. We divide the interval into
    // numPts segments (numPts + 1 points), but then drop the point at z=1 and
    // shift by +1/2 segment--we don't want points at the singular points +/-1.
    Real dz = 2.0 / numPts;
    Real  z = -1 + 0.5 * dz;
    // Note: angular velocity (treating z as time) is proportional to numPts
    Real dtheta = M_PI * (sqrt(5.0) + 1.0);
    Real  theta = 0;
    for (size_t i = 0; i < numPts; ++i) {
        Real r = sqrt(1.0 - z * z);
        pts.emplace_back(r * cos(theta), r * sin(theta), z);
        z += dz;
        theta += dtheta;
    }
}

// Append to "pts" numPts sampled on a sphere of radius "r" and center "c"
template<class OutputCollection, typename Real, class Pt>
void generateSpherePoints(size_t numPts, OutputCollection &pts, Real r, const Pt &c) {
    size_t offset = pts.size();
    generateSpherePoints(numPts, pts);
    for (size_t i = offset; i < pts.size(); ++i) {
        pts[i] *= r;
        pts[i] += c;
    }
}

// Include the 6 vertices tounching the bounding box. This is is used to reduce
// the liklihood of degeneracies in SphereConvexHull, or at least make them more
// exact/detectable.
template<class OutputCollection, typename Real, class Pt>
void generateSpherePointsWithExtremeVertices(size_t numPts, OutputCollection
        &pts, Real r, const Pt &c) {
    generateSpherePoints(numPts - 6, pts, r, c);

    for (size_t d = 0; d < 3; ++d) {
        Pt offset = Pt::Zero();
        offset[d] = r;
        pts.emplace_back(c + offset);
        pts.emplace_back(c - offset);
    }
}

#endif /* end of include guard: SPHEREPOINTS_HH */
