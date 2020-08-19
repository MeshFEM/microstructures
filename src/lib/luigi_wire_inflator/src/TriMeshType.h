#ifndef TRIMESHTYPE_H
#define TRIMESHTYPE_H

#include <vcg/complex/complex.h>

class TFace;
class TVertex;

struct TUsedTypes : public vcg::UsedTypes<vcg::Use<TVertex>::AsVertexType, vcg::Use<TFace>::AsFaceType> {};

class TVertex : public vcg::Vertex<
          TUsedTypes
        , vcg::vertex::Coord3d
        , vcg::vertex::Normal3d
        , vcg::vertex::Qualityd
        , vcg::vertex::VFAdj
        , vcg::vertex::BitFlags
        , vcg::vertex::Mark>{};
class TFace   : public vcg::Face<
          TUsedTypes
        , vcg::face::VertexRef
        , vcg::face::Normal3d
        , vcg::face::BitFlags
        , vcg::face::VFAdj
        , vcg::face::FFAdj> {};
class TMesh   : public vcg::tri::TriMesh< std::vector<TVertex>, std::vector<TFace> > {};

#endif // TRIMESHTYPE_H
