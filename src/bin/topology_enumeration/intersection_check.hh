#ifndef INTERSECTION_CHECK_HH
#define INTERSECTION_CHECK_HH

#include <isosurface_inflator/Symmetry.hh>
#include <MeshFEM/MeshIO.hh>
#include <vector>
#include <utility>

bool hasSelfIntersection(const std::vector<Point3d> &nodes,
                         const std::vector<std::pair<size_t, size_t>> &edges);

#endif /* end of include guard: INTERSECTION_CHECK_HH */
