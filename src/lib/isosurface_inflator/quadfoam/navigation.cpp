#ifdef HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#include "navigation.h"
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace micro_quadfoam {

// Build the navigation data for a given mesh
NavigationData::NavigationData(const Eigen::MatrixXi &F)
	: num_vertex_per_face_((int) F.cols())
	, corner_to_vertex_(F.rows() * num_vertex_per_face_, -1)
	, corner_to_edge_(F.rows() * num_vertex_per_face_, -1)
	, corner_to_adjacent_face_(F.rows() * num_vertex_per_face_, -1)
{
	// Map edges (vi, vj) to a unique edge id ek
	typedef std::pair<int, int> EdgeVertices;
	std::vector<std::pair<EdgeVertices, int>> list_of_edges;

	for (int f = 0; f < F.rows(); ++f) {
		for (int lv = 0; lv < F.cols(); ++lv) {
			int c = face_corner(f, lv);
			corner_to_vertex_[c] = F(f, lv);

			int v1 = F(f, lv);
			int v2 = F(f, (lv+1)%F.cols());
			EdgeVertices e(std::min(v1, v2), std::max(v1, v2));
			list_of_edges.emplace_back(e, c);
		}
	}

	// Sort edges and determine unique edges
	std::sort(list_of_edges.begin(), list_of_edges.end());

	int num_edges = 0;
	std::vector<int> corners_around_edge;
	for (size_t i = 0; i < list_of_edges.size(); ++i) {
		// Accumulate for current edge
		corners_around_edge.push_back(list_of_edges[i].second);

		if (i + 1 == list_of_edges.size() || list_of_edges[i + 1].first != list_of_edges[i].first) {
			// Done with the current edge
			int e = num_edges++;
			for (int c : corners_around_edge) {
				corner_to_edge_[c] = e;
			}
			assert(corners_around_edge.size() <= 2 && "Non manifold edges detected!");
			if (corners_around_edge.size() == 2) {
				int c1 = corners_around_edge[0];
				int c2 = corners_around_edge[1];
				int f1 = c1 / num_vertex_per_face_;
				int f2 = c2 / num_vertex_per_face_;
				corner_to_adjacent_face_[c1] = f2;
				corner_to_adjacent_face_[c2] = f1;
				// std::cout << c1 << ' ' << c2 << " -> " << f1 << ' ' << f2 << std::endl;
			}
			corners_around_edge.clear();
			auto uv = list_of_edges[i].first;
			edge_to_vertex_.push_back({{ uv.first, uv.second }});
			// std::cout << "e" << e << ": " << uv.first << "-" << uv.second << std::endl;
		}
	}
}

// Retrieves the index (v,e,f) of one vertex incident to the given face
NavigationIndex index_from_face(const Eigen::MatrixXi &F, const NavigationData &data, int f, int lv) {
	NavigationIndex idx;
	idx.face_corner = data.face_corner(f, lv);
	idx.vertex = data.corner_vertex(idx.face_corner);
	idx.face = f;
	idx.edge = data.corner_edge(idx.face_corner);
	int v2 = F(f, (lv+1) % F.cols());
	if (switch_vertex(data, idx).vertex != v2) {
		assert(false);
		return switch_edge(data, idx);
	}
	return idx;
}

// -----------------------------------------------------------------------------

// Navigation in a surface mesh
NavigationIndex switch_vertex(const NavigationData &data, NavigationIndex idx) {
	int c1 = data.next_corner_around_face(idx.face, idx.face_corner);
	int v1 = data.corner_vertex(c1);
	if (v1 == data.edge_vertex(idx.edge, 0) || v1 == data.edge_vertex(idx.edge, 1)) {
		idx.face_corner = c1;
		idx.vertex = v1;
		return idx;
	} else {
		int c2 = data.prev_corner_around_face(idx.face, idx.face_corner);
		int v2 = data.corner_vertex(c2);
		assert(v2 == data.edge_vertex(idx.edge, 0) || v2 == data.edge_vertex(idx.edge, 1));
		idx.face_corner = c2;
		idx.vertex = v2;
		return idx;
	}
}

NavigationIndex switch_edge(const NavigationData &data, NavigationIndex idx) {
	int v2 = data.edge_vertex(idx.edge, 0);
	if (v2 == idx.vertex) {
		v2 = data.edge_vertex(idx.edge, 1);
	}
	int c1 = data.next_corner_around_face(idx.face, idx.face_corner);
	int v1 = data.corner_vertex(c1);
	if (v1 == v2) {
		int c2 = data.prev_corner_around_face(idx.face, idx.face_corner);
		idx.edge = data.corner_edge(c2);
		return idx;
	} else {
		idx.edge = data.corner_edge(idx.face_corner);
		return idx;
	}
}

NavigationIndex switch_face(const NavigationData &data, NavigationIndex idx) {
	int c1 = idx.face_corner;
	if (data.corner_edge(c1) != idx.edge) {
		c1 = data.prev_corner_around_face(idx.face, c1);
	}
	int f2 = data.corner_adjacent_face(c1);
	if (f2 < 0) {
		idx.face = -1;
		return idx;
	} else {
		// Iterate over all corners of the new face until we find the vertex we came from.
		// Not ideal but for now it will do the job.
		for (int c2 = data.face_corners_begin(f2); c2 < data.face_corners_end(f2); ++c2) {
			int v2 = data.corner_vertex(c2);
			if (v2 == (int) idx.vertex) {
				idx.face = f2;
				idx.face_corner = c2;
				return idx;
			}
		}
		std::cout << "Not found" << std::endl;
		assert(false); // This should not happen
		return idx;
	}
}

} // namespace micro_quadfoam

#endif
