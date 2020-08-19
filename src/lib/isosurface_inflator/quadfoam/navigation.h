#ifdef HAS_LIBIGL

#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <vector>
#include <array>
////////////////////////////////////////////////////////////////////////////////

namespace micro_quadfoam {

struct NavigationIndex {
	int vertex;
	int edge;
	int face;
	int face_corner;
};

struct NavigationData {
	int num_vertex_per_face_;
	std::vector<int> corner_to_vertex_;
	std::vector<int> corner_to_edge_;
	std::vector<int> corner_to_adjacent_face_;
	std::vector<std::array<int, 2>> edge_to_vertex_;

	NavigationData() = default;
	NavigationData(const Eigen::MatrixXi &F);

	int corner_vertex(int c) const { return corner_to_vertex_[c]; }
	int corner_edge(int c) const { return corner_to_edge_[c]; }
	int corner_adjacent_face(int c) const { return corner_to_adjacent_face_[c]; }

	int edge_vertex(int e, int lv) const{ return edge_to_vertex_[e][lv]; }
	int face_corner(int f, int lv) const { return f * num_vertex_per_face_ + lv; }

	int face_corners_begin(int f) const { return f * num_vertex_per_face_; }
	int face_corners_end(int f) const { return (f + 1) * num_vertex_per_face_; }

	int next_corner_around_face(int f, int c) const { return face_corners_begin(f) + (c + 1) % num_vertex_per_face_; }
	int prev_corner_around_face(int f, int c) const { return face_corners_begin(f) + (c + num_vertex_per_face_ - 1) % num_vertex_per_face_; }
};

// Retrieves the index (v,e,f) of one vertex incident to the given face
NavigationIndex index_from_face(const Eigen::MatrixXi &F, const NavigationData &data, int f, int lv);

// Navigation in a surface mesh
NavigationIndex switch_vertex(const NavigationData &data, NavigationIndex idx);
NavigationIndex switch_edge(const NavigationData &data, NavigationIndex idx);
NavigationIndex switch_face(const NavigationData &data, NavigationIndex idx);

// Iterate in a mesh
inline NavigationIndex next_around_face(const NavigationData &data, NavigationIndex idx) { return switch_edge(data, switch_vertex(data, idx)); }
inline NavigationIndex next_around_edge(const NavigationData &data, NavigationIndex idx) { return switch_vertex(data, switch_face(data, idx)); } // TODO: Fixme
inline NavigationIndex next_around_vertex(const NavigationData &data, NavigationIndex idx) { return switch_face(data, switch_edge(data, idx)); } // TODO: Fixme

} // namespace micro_quadfoam

#endif
