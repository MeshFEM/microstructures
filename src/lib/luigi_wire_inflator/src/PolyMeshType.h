#ifndef POLYMESHTYPE_H
#define POLYMESHTYPE_H

#include <vcg/complex/complex.h>

// include the support for polygon meshes (function to convert from/to trimesh)
#include <vcg/complex/algorithms/polygon_support.h>

// include the support for polygon meshes (the component for the face )
#include <vcg/simplex/face/component_polygon.h>

// Definition of a mesh of polygons
class PolyFace;
class PolyVertex;
struct PolyUsedTypes : public vcg::UsedTypes<vcg::Use<PolyVertex>::AsVertexType,
                                             vcg::Use<PolyFace>::AsFaceType> {};

class PolyVertex : public vcg::Vertex<
          PolyUsedTypes
        , vcg::vertex::Coord3d
        , vcg::vertex::Normal3d
        , vcg::vertex::BitFlags
        , vcg::vertex::Mark
        > {};

class PolyFace : public vcg::Face<
          PolyUsedTypes
        , vcg::face::PolyInfo // this is necessary if you use component in vcg/simplex/face/component_polygon.h
                              // It says "this class is a polygon and the memory for its components (e.g. pointer to its vertices
                              // will be allocated dynamically")
        , vcg::face::PFVAdj   // Pointer to the vertices (just like FVAdj )
        , vcg::face::Mark
        , vcg::face::PFFAdj   // Pointer to edge-adjacent face (just like FFAdj )
        , vcg::face::BitFlags // bit flags
        , vcg::face::Normal3d // normal
        > {};

class PolyMesh: public vcg::tri::TriMesh<
          std::vector<PolyVertex> // the vector of vertices
        , std::vector<PolyFace>   // the vector of faces
        >{};

#endif // POLYMESHTYPE_H
