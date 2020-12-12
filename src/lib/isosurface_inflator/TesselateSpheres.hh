////////////////////////////////////////////////////////////////////////////////
// TesselateSpheres.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Tesselates a collection of spheres using triangles.
//      Each sphere is tesselated by computing the convex hull of N points
//      sampled on the sphere.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/13/2016 11:55:34
////////////////////////////////////////////////////////////////////////////////
#ifndef TESSELATESPHERES_HH
#define TESSELATESPHERES_HH
#include "SpherePoints.hh"
#include "ConvexHullTriangulation.hh"

template<typename Real>
void tesselateSpheres(size_t nSpherePts,
        const std::vector<Point3<Real>> &centers,
        const std::vector<Real>         &radii,
        std::vector<MeshIO::IOVertex > &outVertices,
        std::vector<MeshIO::IOElement> &outElements,
        // Which sphere each vertex came from:
        std::vector<size_t> &sphereIndex) {
    const size_t numSpheres = centers.size();
    assert(radii.size() == numSpheres);
    // Generate unit sphere at the origin
    std::vector<Point3<double>> spherePts;
    generateSpherePoints(nSpherePts, spherePts);
    std::vector<MeshIO::IOVertex > unitSphereV;
    std::vector<MeshIO::IOElement> unitSphereE;
    convexHullFromTriangulation(spherePts, unitSphereV, unitSphereE);

    outVertices.clear(), outElements.clear(), sphereIndex.clear();
    outVertices.reserve(numSpheres * unitSphereV.size());
    outElements.reserve(numSpheres * unitSphereE.size());
    sphereIndex.reserve(outVertices.size());

    for (size_t i = 0; i < numSpheres; ++i) {
        size_t offset = outVertices.size();
        for (const auto &v : unitSphereV) {
            outVertices.push_back((v.point * radii[i] + centers[i]).eval());
            sphereIndex.push_back(i);
        }
        for (const auto &e : unitSphereE) {
            outElements.push_back(e);
            for (size_t &j : outElements.back()) j += offset;
        }
    }
}

#endif /* end of include guard: TESSELATESPHERES_HH */

