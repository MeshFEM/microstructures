////////////////////////////////////////////////////////////////////////////////
// ConvexHullTriangulation.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes the convex hull using the 3D triangulations package of CGAL.
//      The advantage of this approach over using convex_hull_3 is that it is
//      easier to customize the per-vertex data. Instead of needing to create a
//      new kernel with a custom point type, we can just use
//      Triangulation_vertex_base_with_info_3
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/05/2016 17:20:32
////////////////////////////////////////////////////////////////////////////////
#ifndef CONVEXHULLTRIANGULATION_HH
#define CONVEXHULLTRIANGULATION_HH

#include <MeshFEM/MeshIO.hh>

#include <cassert>
#include <sstream>
#include <vector>

template<class PointCollection>
void convexHullFromTriangulation(const PointCollection &points,
                std::vector<MeshIO::IOVertex > &hullVertices,
                std::vector<MeshIO::IOElement> &hullElements,
                std::vector<size_t>            &originatingVertexIndices);

template<class PointCollection>
void convexHullFromTriangulation(const PointCollection &points,
                std::vector<MeshIO::IOVertex > &hullVertices,
                std::vector<MeshIO::IOElement> &hullElements) {
    std::vector<size_t> dummy;
    convexHullFromTriangulation<PointCollection>(points, hullVertices, hullElements, dummy);
}

#endif /* end of include guard: CONVEXHULLTRIANGULATION_HH */
