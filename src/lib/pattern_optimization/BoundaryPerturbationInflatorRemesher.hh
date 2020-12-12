////////////////////////////////////////////////////////////////////////////////
// BoundaryPerturbationInflatorRemesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Remeshes the mesh created by a BoundaryPerturbationInflator. In
//      particular, short boundary segments are merged, and the interior is
//      remeshed using triangle.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/28/2015 17:56:39
////////////////////////////////////////////////////////////////////////////////
#ifndef BOUNDARYPERTURBATIONINFLATORREMESHER_HH
#define BOUNDARYPERTURBATIONINFLATORREMESHER_HH

#include "PatternOptimizationConfig.hh"
#include "BoundaryPerturbationInflator.hh"
#include <MeshFEM/filters/CurveCleanup.hh>
#include <MeshFEM/filters/extract_hole_boundaries.hh>
#include <MeshFEM/Geometry.hh>
#include <MeshFEM/Triangulate.h>
#include <MeshFEM/utils.hh>
#include <iostream>
#include <vector>
#include <queue>

template<size_t N>
void remeshPerturbedShape(const BoundaryPerturbationInflator<N> &m,
                          Real maxVolume,
                          std::vector<MeshIO::IOVertex>  &outVertices,
                          std::vector<MeshIO::IOElement> &outElements);

template<>
void remeshPerturbedShape(const BoundaryPerturbationInflator<2> &m,
                          Real maxVolume,
                          std::vector<MeshIO::IOVertex>  &outVertices,
                          std::vector<MeshIO::IOElement> &outElements) {
    const auto &mesh = m.mesh();
    std::vector<Real> bdryLengths;
    size_t numBE = mesh.numBoundaryElements();
    bdryLengths.reserve(numBE);
    for (auto be : mesh.boundaryElements()) {
        auto p0 = be.vertex(0).volumeVertex().node()->p;
        auto p1 = be.vertex(1).volumeVertex().node()->p;
        bdryLengths.push_back((p1 - p0).norm());
    }
    std::vector<size_t> order;
    sortPermutation(bdryLengths, order);
    Real medianLen = bdryLengths.at(order.at(numBE / 2));

    // First, determine points in each hole. Assumes that hole boundary point
    // centroid lies within the hole, which should usually be true...)
    // TODO: better approach (based on winding number/shooting rays?)
    //      Don't tell triangle about holes; remove them ourselves
    //      note: this will allow triangle to insert Steiner points on the hole
    //      boundary. The downside is unneccessary subdivision may occur on the
    //      hole boundaries to mesh the void.
    std::vector<std::vector<size_t>> holeBoundaries;
    extract_hole_boundaries(mesh, holeBoundaries);
    std::vector<Point2D> holes;
    for (const auto &hb : holeBoundaries) {
        Point2D centroid(Point2D::Zero());
        for (size_t bei : hb)
            centroid += mesh.boundaryElement(bei).vertex(0).volumeVertex().node()->p;
        centroid /= hb.size();
        holes.push_back(centroid);
    }

    const auto &config = PatternOptimization::Config::get();
    Real minLength = medianLen * config.remeshMergeThreshold;
    Real maxLength = medianLen * config.remeshSplitThreshold;

    // Copy over existing boundary
    outVertices.clear(), outVertices.reserve(mesh.numBoundaryVertices());
    outElements.clear(), outElements.reserve(mesh.numBoundaryElements());
    for (auto bv : mesh.boundaryVertices()) outVertices.emplace_back(bv.volumeVertex().node()->p);
    for (auto be : mesh.boundaryElements()) outElements.emplace_back(be.vertex(0).index(), be.vertex(1).index());

    curveCleanup(outVertices, outElements, outVertices, outElements,
                 minLength, maxLength, config.remeshFeatureAngleThreshold, true);

    // Remesh the interior.
    // Q: quiet, Y: do not remesh the boundary
    triangulatePSLG(outVertices, outElements, holes, outVertices, outElements,
                    maxVolume, "QY");
}

template<>
void remeshPerturbedShape(const BoundaryPerturbationInflator<3> &/* m */,
                          Real /* maxVolume */,
                          std::vector<MeshIO::IOVertex>  &/* outVertices */,
                          std::vector<MeshIO::IOElement> &/* outElements */) {
    // 3D is a bit harder...
    throw std::runtime_error("remeshPerturbedShape<3> unimplemented.");
}

#endif /* end of include guard: BOUNDARYPERTURBATIONINFLATORREMESHER_HH */
