////////////////////////////////////////////////////////////////////////////////
// rasterize.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Implements rasterization of a signed distance region to a grid
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  11/17/2017 20:44:07
////////////////////////////////////////////////////////////////////////////////
#include "rasterize.hh"
#include <MeshFEM/filters/gen_grid.hh>

// Rasterizes the z = 0 midplane if 2 sizes are passed,
// full bounding box if 3 sizes are passed.
void rasterize(const SignedDistanceRegion<3> &sdf,
        std::vector<size_t> sizes,
        std::vector<MeshIO::IOVertex> &vertices,
        std::vector<MeshIO::IOElement> &elements,
        ScalarField<Real> &indicator) {
    // Generate rasterization grid filling the sdf bounding box
    if      (sizes.size() == 2) gen_grid(sizes[0], sizes[1],           vertices, elements);
    else if (sizes.size() == 3) gen_grid(sizes[0], sizes[1], sizes[2], vertices, elements);
    else throw std::runtime_error("Expected 2 or 3 rasterization dimensions.");

    const auto &bb = sdf.boundingBox();
    auto scale = bb.dimensions();
    for (size_t i = 0; i < sizes.size(); ++i) scale[i] /= sizes[i];
    for (auto &v : vertices) { v.point = (v.point.array() * scale.array()).matrix() + bb.minCorner; }

    // Run inside/outside check for each quad/hex center to generate indicator field.
    indicator.resizeDomain(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
        Point3D center(Point3D::Zero());
        for (size_t c : elements[i])
            center += vertices[c].point;
        center *= 1.0 / elements[i].size();

        indicator[i] = sdf.isInside(center) ? 1.0 : 0.0;
    }
}
