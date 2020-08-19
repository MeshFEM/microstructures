#ifndef OUTMESH_H
#define OUTMESH_H

#include <array>
#include <vector>
#include <unordered_map>
#include "edge_hash.h"

template <size_t Dim = 3, size_t NodesPerElement = 3>
struct OutMesh
{
	typedef int                                    IndexType;
	typedef double                                 ScalarType;
	typedef std::array<ScalarType, Dim>            CoordType;
	typedef std::array<ScalarType, Dim>            NormalType;
	typedef std::array<IndexType, NodesPerElement> FaceType;
	typedef std::pair<IndexType, IndexType>        EdgeType;

	typedef std::vector<CoordType>  NodeVector;
	typedef std::vector<FaceType>   ElementVector;

	typedef std::vector<ScalarType>                         Fields;
	typedef std::unordered_map<EdgeType, Fields, edge_hash> EdgeFields;

	NodeVector    nodes;
	ElementVector elements;
	EdgeFields    edge_fields;

	// Julian: per-vertex vector velocities instead of the old per-edge average
	// normal velocities in "edge_fields"
	using VField = std::vector<CoordType>;
	using VertexVelocities = std::vector<VField>;
	VertexVelocities vertex_velocities;
};

#endif // OUTMESH_H
