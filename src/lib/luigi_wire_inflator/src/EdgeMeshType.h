#ifndef EDGEMESHTYPE_H
#define EDGEMESHTYPE_H

#include <vcg/complex/complex.h>

class EmEdgeType;
class EmVertexType;

struct EUsedTypes : public vcg::UsedTypes< vcg::Use<EmVertexType>::AsVertexType,
                                           vcg::Use<EmEdgeType>::AsEdgeType> {};

class EmVertexType : public vcg::Vertex< EUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::Qualityd, vcg::vertex::BitFlags, vcg::vertex::VEAdj> {};
class EmEdgeType   : public vcg::Edge< EUsedTypes, vcg::edge::VertexRef, vcg::edge::Qualityd, vcg::edge::BitFlags, vcg::edge::EEAdj, vcg::edge::VEAdj> {};
class EMesh       : public vcg::tri::TriMesh< std::vector<EmVertexType>, std::vector<EmEdgeType> > {};

#endif // EDGEMESHTYPE_H
