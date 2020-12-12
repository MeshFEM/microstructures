////////////////////////////////////////////////////////////////////////////////
// PostProcess.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Postprocess/clean up the geometry from the mesher and then compute
//  normals and shape velocities.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  08/23/2017 16:49:49
////////////////////////////////////////////////////////////////////////////////
#ifndef POSTPROCESS_HH
#define POSTPROCESS_HH

#include "IsosurfaceInflator.hh"
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/Geometry.hh>
#include <vector>

// Postprocess: Snap to base cell and then reflect if necessary Compute vertex normals and normal
// shape velocities
//
// @param[int,out] vertices               { Input mesh vertices positions }
// @param[int,out] elements               { Input mesh elements (triangles or tets) }
// @param[out]     normalShapeVelocities  { The normal shape velocities }
// @param[out]     vertexNormals          { The vertex normals }
// @param[in]      inflator               { The inflator }
// @param[in]      meshedFullPeriodCell   { ??? }
// @param[in]      requestFullPeriodCell  { ??? }
// @param[in]      meshCell               { Bounding box }
// @param[in]      opts                   { Meshing options }
// @param[in]      cheapPostProcessing    { Enable cheap post processing (???) }
// @param[in]      nonPeriodicity         { Whether the cell is not periodic }
//
// @tparam         N                      { Dimension }
// @tparam         Point                  { Point type }
//
template<size_t N, class Point>
void postProcess(std::vector<MeshIO::IOVertex>  &vertices,
                 std::vector<MeshIO::IOElement> &elements,
                 std::vector<std::vector<Real>> &normalShapeVelocities,
                 std::vector<Point>             &vertexNormals,
                 const IsosurfaceInflator::Impl &inflator,
                 bool                      meshedFullPeriodCell,
                 bool                      requestFullPeriodCell,
                 const BBox<Point>         &meshCell,
                 const MeshingOptions      &opts,
                 bool cheapPostProcessing,
                 bool nonPeriodicity);

#endif /* end of include guard: POSTPROCESS_HH */
