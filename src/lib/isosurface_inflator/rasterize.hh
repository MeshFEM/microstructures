////////////////////////////////////////////////////////////////////////////////
// rasterize.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Implements rasterization of a signed distance region to a grid
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  11/17/2017 20:42:14
////////////////////////////////////////////////////////////////////////////////
#ifndef RASTERIZE_HH
#define RASTERIZE_HH

#include "SignedDistanceRegion.hh"
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/Fields.hh>


// Rasterizes the z = 0 midplane if 2 sizes are passed,
// full bounding box if 3 sizes are passed.
void rasterize(const SignedDistanceRegion<3> &sdf,
        std::vector<size_t> sizes,
        std::vector<MeshIO::IOVertex> &vertices,
        std::vector<MeshIO::IOElement> &elements,
        ScalarField<Real> &indicator);

#endif /* end of include guard: RASTERIZE_HH */
