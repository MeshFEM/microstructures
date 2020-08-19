#ifdef HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#include "instantiate.h"
#include "navigation.h"
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_unreferenced.h>
#include <igl/is_boundary_edge.h>
#include <igl/edges.h>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace micro_quadfoam {

////////////////////////////////////////////////////////////////////////////////

namespace {

enum class PatternType {
	NOT_PERIODIC,    // Opposite borders do not match
	NOT_SYMMETRIC,   // At least one border is not symmetric along its bisector
	SIMPLE_PERIODIC, // Opposite borders match individually
	DOUBLE_PERIODIC, // All four borders match between each others
};

// -----------------------------------------------------------------------------

// Determine the type of periodicity of the given pattern. This function rely
// only on the input vertex positions, but does not check if the underlying
// topology is consistent.
//
// @param[in]  V                { #V x 3 matrix of vertex positions }
// @param[out] border_vertices  { List of vertices along each side }
//
// @return     { Type of periodicity of the pattern. }
//
PatternType compute_pattern_type(const Eigen::MatrixXd &V,
	std::array<Eigen::VectorXi, 4> &border_vertices)
{
	Eigen::Vector2d lower = V.colwise().minCoeff().head<2>();
	Eigen::Vector2d upper = V.colwise().maxCoeff().head<2>();

	std::array<Eigen::MatrixXd, 2> border;

	// Compute vertices along the border. Sides are numbered as follows:
	//
	//        2
	//   v3 ----- v2
	//   |        |
	// 3 |        | 1
	//   |        |
	//   v0 ----- v1
	//        0
	// ↑y
	// └─→x
	for (int d = 0; d < 2; ++d) {
		std::vector<std::pair<double, int>> low, upp;
		for (int i = 0; i < V.rows(); ++i) {
			if (V(i, 1-d) == lower[1-d]) { low.emplace_back(V(i, d), i); }
			if (V(i, 1-d) == upper[1-d]) { upp.emplace_back(V(i, d), i); }
		}
		std::sort(low.begin(), low.end());
		std::sort(upp.begin(), upp.end());
		border_vertices[d?3:0].resize(low.size());
		for (size_t i = 0; i < low.size(); ++i) {
			border_vertices[d?3:0][i] = low[i].second;
		}
		border_vertices[d?1:2].resize(upp.size());
		for (size_t i = 0; i < upp.size(); ++i) {
			border_vertices[d?1:2][i] = upp[i].second;
		}
	}
	for (int i : {2, 3}) {
		border_vertices[i].reverseInPlace();
	}

	// Check if borders have the same topology and geometry
	for (int d = 0; d < 2; ++d) {
		std::vector<double> low, upp;
		for (int i = 0; i < V.rows(); ++i) {
			if (V(i, d) == lower[d]) { low.push_back(V(i, 1-d)); }
			if (V(i, d) == upper[d]) { upp.push_back(V(i, 1-d)); }
		}
		std::sort(low.begin(), low.end());
		std::sort(upp.begin(), upp.end());
		if (low.size() != upp.size()) {
			return PatternType::NOT_PERIODIC;
		}
		size_t n = low.size();
		border[d].resize(n, 2);
		for (size_t i = 0; i < n; ++i) {
			border[d](i, 0) = low[i];
			border[d](i, 1) = upp[i];
		}
		if (!border[d].col(0).isApprox(border[d].col(1))) {
			return PatternType::NOT_PERIODIC;
		}
	}

	// Check if borders are symmetric
	for (int d = 0; d < 2; ++d) {
		if (!(border[d].col(0).array() - lower[d]).isApprox(
			upper[d] - border[d].col(0).array().reverse()))
		{
			return PatternType::NOT_SYMMETRIC;
		}
	}

	// Check if horizontal and vertical borders are matching
	if (border[0].size() == border[1].size()) {
		if (!border[0].col(0).isApprox(border[1].col(0))) {
			std::cerr << "Warning: pattern boundaries have the same number of vertices, but their position differ slighly." << std::endl;
			Eigen::MatrixXd X(border[0].size(), 2);
			X.col(0) = border[0].col(0);
			X.col(1) = border[1].col(0);
			std::cout << X << std::endl << std::endl;;
		}
		return PatternType::DOUBLE_PERIODIC;
	} else {
		return PatternType::SIMPLE_PERIODIC;
	}

	return PatternType::SIMPLE_PERIODIC;
}

// -----------------------------------------------------------------------------

Eigen::VectorXi vertex_degree(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
	Eigen::VectorXi count = Eigen::VectorXi::Zero(V.rows());

	for (unsigned i=0; i<F.rows();++i) {
		for (unsigned j=0; j<F.cols();++j) {
			// Avoid duplicate edges
			if (F(i,j) < F(i,(j+1)%F.cols())) {
				count(F(i,j  )) += 1;
				count(F(i,(j+1)%F.cols())) += 1;
			}
		}
	}

	return count;
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

bool instantiate_pattern(
	const Eigen::MatrixXd &IV,
	const Eigen::MatrixXi &IQ,
	const std::vector<Eigen::MatrixXd> &PV,
	const std::vector<Eigen::MatrixXi> &PF,
	Eigen::MatrixXd &OV,
	Eigen::MatrixXi &OF,
	bool remap_duplicated_vertices,
	double tolerance,
	Eigen::VectorXi *remap_vertices,
	Eigen::VectorXi *parent_face,
	Eigen::MatrixXi *boundary_edges)
{
	assert(PV.size() == (size_t) IQ.rows());
	assert(PF.empty() || PF.size() == (size_t) IQ.rows());
	size_t num_quads = IQ.rows();

	// List of vertices along each border (from lv to (lv+1)%4)
	std::vector<std::array<Eigen::VectorXi, 4>> border_vertices(num_quads);
	for (size_t q = 0; q < num_quads; ++q) {
		compute_pattern_type(PV[q], border_vertices[q]);
		// std::cout << PV[q] << std::endl;
		// for (auto x : border_vertices[q]) {
		// 	std::cout << x << std::endl;
		// }
		// std::cout << std::endl;
	}

	std::function<const Eigen::MatrixXi *(int)> PF_func;
	if (!PF.empty()) { PF_func = [&](int q) { return &PF[q]; }; }

	// Call generic function
	return instantiate_pattern_aux(
		IV, IQ,
		[&](int q) { return &PV[q]; },
		PF_func,
		[&](int q) { return border_vertices[q]; },
		OV, OF,
		remap_duplicated_vertices,
		tolerance,
		remap_vertices,
		parent_face,
		boundary_edges
	);
}

////////////////////////////////////////////////////////////////////////////////

bool instantiate_pattern_aux(
	const Eigen::MatrixXd &IV,
	const Eigen::MatrixXi &IF,
	std::function<const Eigen::MatrixXd *(int)> PV,
	std::function<const Eigen::MatrixXi *(int)> PF,
	std::function<std::array<Eigen::VectorXi, 4>(int)> border_vertices,
	Eigen::MatrixXd &OV,
	Eigen::MatrixXi &OF,
	bool remap_duplicated_vertices,
	double tolerance,
	Eigen::VectorXi *vertex_id_map_ptr,
	Eigen::VectorXi *parent_face_ptr,
	Eigen::MatrixXi *boundary_edges_ptr)
{
	int num_quads = (int) IF.rows();

	auto valence = vertex_degree(IV, IF);

	// Instantiate (duplicating vertices)
	int num_vertices = 0;
	int num_facets = 0;
	int num_corners = -1;
	for (int q = 0; q < num_quads; ++q) {
		num_vertices += (int) PV(q)->rows();
		if (PF) {
			num_facets += (int) PF(q)->rows();
			int size = (int) PF(q)->cols();
			if (num_corners < 0) { num_corners = size; }
			else { assert(num_corners == size); }
		} else {
			num_corners = 0;
		}
	}
	Eigen::MatrixXd V(num_vertices, 3);
	Eigen::MatrixXi F(num_facets, num_corners);
	Eigen::VectorXi parent_face(num_facets);
	Eigen::VectorXi vertex_offset(num_quads);

	for (int q = 0, v0 = 0, f0 = 0; q < num_quads; ++q) {
		// Normalized coordinates (between 0 and 1 for barycentric coordinates)
		auto lower = PV(q)->colwise().minCoeff();
		auto upper = PV(q)->colwise().maxCoeff();
		Eigen::MatrixXd PVN = (PV(q)->rowwise() - lower).array().rowwise()
			/ (upper - lower).cwiseMax(1e-5).array();
		const auto & u = PVN.col(0).array();
		const auto & v = PVN.col(1).array();
		Eigen::RowVector3d a = IV.row(IF(q, 0));
		Eigen::RowVector3d b = IV.row(IF(q, 1));
		Eigen::RowVector3d c = IV.row(IF(q, 2));
		Eigen::RowVector3d d = IV.row(IF(q, 3));
		V.middleRows(v0, PVN.rows()) = ((1-u)*(1-v)).matrix()*a
			+ (u*(1-v)).matrix()*b
			+ (u*v).matrix()*c
			+ ((1-u)*v).matrix()*d;
		if (PVN.cols() > 2) {
			V.middleRows(v0, PVN.rows()).col(2) += PVN.col(2);
		}
		if (PF) {
			F.middleRows(f0, PF(q)->rows()) = PF(q)->array() + v0;
			parent_face.segment(f0, PF(q)->rows()).setConstant(q);
		}
		vertex_offset(q) = v0;
		v0 += (int) PV(q)->rows();
		if (PF) { f0 += (int) PF(q)->rows(); }
	}

	// Remapped vertex id, with union-find (after duplicate removal)
	std::vector<int> parent(V.rows());
	std::iota(parent.begin(), parent.end(), 0);
	std::function<int(int)> find = [&](int i) {
		if (parent[i] == i) { return i; }
		parent[i] = find(parent[i]);
		return parent[i];
	};

	auto merge = [&](int i, int j) {
		parent[find(i)] = find(j);
	};

	// Navigation data
	NavigationData data(IF);

	// [Helper] Retrieve vertices located along a given edge
	auto edge_vertices = [&](NavigationIndex idx0) {
		NavigationIndex idx = index_from_face(IF, data, idx0.face, 0);
		int lv;
		for (lv = 0; lv < 4; ++lv) {
			if (idx.edge == idx0.edge) { break; }
			idx = next_around_face(data, idx);
		}
		assert(idx.edge == idx0.edge);
		int q = idx0.face;
		Eigen::VectorXi side = border_vertices(q)[lv].array() + vertex_offset(q);
		if (idx.vertex != idx0.vertex) {
			side.reverseInPlace();
		}
		return side;
	};

	// [Helper] Stitch vertices from adjacent quads
	bool incompatible_num = false;
	bool incompatible_pos = false;
	auto remap_vertices = [&](NavigationIndex idx1, NavigationIndex idx2) {
		assert(idx1.edge == idx2.edge);
		if (idx1.vertex != idx2.vertex) { idx1 = switch_vertex(data, idx1); }
		auto side1 = edge_vertices(idx1);
		auto side2 = edge_vertices(idx2);
		if (side1.size() != side2.size()) {
			incompatible_num = true;
			return;
		}
		if (tolerance >= 0) {
			// Optionally check that positions are matching
			for (int i = 0; i < (int) side1.size(); ++i) {
				if ((V.row(side1[i]) - V.row(side2[i])).squaredNorm() > tolerance * tolerance) {
					incompatible_pos = true;
					return;
				}
			}
		}
		for (int i = 0; i < (int) side1.size(); ++i) {
			const int x1 = side1[i];
			const int x2 = side2[i];
			merge(x1, x2);
		}
	};

	// Iterate over quads and remap vertices on adjacent quads
	for (int q = 0; q < IF.rows(); ++q) {
		auto idx = index_from_face(IF, data, q, 0);
		for (int lv = 0; lv < 4; ++lv) {
			auto idx2 = switch_face(data, idx);
			if (idx2.face >= 0) {
				remap_vertices(idx, idx2);
			}
			idx = next_around_face(data, idx);
		}
	}

	// Assign ids to vertices based on their groups
	Eigen::VectorXi ids(V.rows());
	{
		int cnt = 0;
		ids.setConstant(-1);
		for (int v = 0; v < V.rows(); ++v) {
			if (find(v) == v) { ids(v) = cnt++; }
		}
		for (int v = 0; v < V.rows(); ++v) {
			ids(v) = ids(find(v));
		}
	}
	if (vertex_id_map_ptr) {
		if (incompatible_num) {
			std::cerr << "WARNING: Incompatible number of vertices between adjacent patterns!" << std::endl;
			return false;
		}
		if (incompatible_pos) {
			std::cerr << "WARNING: Mismatched vertex positions between adjacent patterns!" << std::endl;
			return false;
		}
	}

	// Retrieve edges on the boundary of the original mesh
	typedef std::pair<int, int> Edge;
	std::vector<Edge> boundary_edges;
	if (boundary_edges_ptr && PF) {
		for (int q = 0; q < num_quads; ++q) {
			Eigen::MatrixXi E;
			igl::edges(*PF(q), E);

			// Boundary of the pattern mesh
			Eigen::Matrix<bool, Eigen::Dynamic, 1> BE;
			igl::is_boundary_edge(E, *PF(q), BE);

			// Vertices on the boundary of the original mesh
			std::vector<bool> boundary_vertex(PV(q)->rows(), false);
			auto idx = index_from_face(IF, data, q, 0);
			for (int lv = 0; lv < 4; ++lv) {
				auto idx2 = switch_face(data, idx);
				if (idx2.face < 0) {
					auto side = edge_vertices(idx);
					side.array() -= vertex_offset(q);
					for (int i = 0; i < side.size(); ++i) {
						boundary_vertex[side[i]] = true;
					}
				}
				idx = next_around_face(data, idx);
			}

			// Keep only pattern boundary edges whose endpoints are on a
			// boundary edge of the original mesh
			for (int le = 0; le < E.rows(); ++le) {
				int v1 = E(le, 0);
				int v2 = E(le, 1);
				if (BE(le) && boundary_vertex[v1] && boundary_vertex[v2]) {
					boundary_edges.emplace_back(v1 + vertex_offset(q), v2 + vertex_offset(q));
				}
			}
		}
	}

	// Remap vertices according to 'remap'
	if (remap_duplicated_vertices) {
		int num_remapped = ids.maxCoeff() + 1;
		OV.resize(num_remapped, V.cols());
		Eigen::VectorXi count(num_remapped);
		OV.setZero();
		count.setZero();
		for (int v = 0; v < V.rows(); ++v) {
			OV.row(ids(v)) += V.row(v);
			count(ids(v))++;
		}
		OV = OV.array().colwise() / count.array().cast<double>();
		OF = F.unaryExpr([&](int v){ return ids(v); });
		for (auto &e : boundary_edges) {
			e.first = ids(e.first);
			e.second = ids(e.second);
		}
		std::sort(boundary_edges.begin(), boundary_edges.end());
		auto it = std::unique(boundary_edges.begin(), boundary_edges.end());
		boundary_edges.resize(std::distance(boundary_edges.begin(), it));
	} else {
		OV = V;
		OF = F;
	}

	// Assign optional output arguments
	if (vertex_id_map_ptr) {
		*vertex_id_map_ptr = ids;
	}
	if (parent_face_ptr) {
		*parent_face_ptr = parent_face;
	}
	if (boundary_edges_ptr) {
		Eigen::MatrixXi &E = *boundary_edges_ptr;
		E.resize(boundary_edges.size(), 2);
		for (int e = 0; e < E.rows(); ++e) {
			E.row(e) << boundary_edges[e].first, boundary_edges[e].second;
		}
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace micro_quadfoam

#endif
