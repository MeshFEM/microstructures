#ifdef HAS_LIBIGL

#pragma once
////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace micro_quadfoam {

///
/// Instantiate a periodic 2D pattern (triangle-mesh) on a given quad mesh
///
/// @param[in]  IV                         { #IV x 3 input quad mesh vertices }
/// @param[in]  IF                         { #IF x 4 input quad mesh facets }
/// @param[in]  PV                         { quad index -> #PV x (2|3) input pattern vertices in
///                                        [0,1]^2 }
/// @param[in]  PF                         { quad index -> #PF x (3|4) input pattern facets }
/// @param[in]  border_vertices            { quad index -> list of vertices along each local border
///                                        (from lv to (lv+1)%4) }
/// @param[out] OV                         { #OV x 3 output mesh vertices }
/// @param[out] OF                         { #OF x 3 output mesh facets }
/// @param[in]  remap_duplicated_vertices  { Should duplicated vertices be remapped in the output }
/// @param[in]  tolerance                  { Tolerance for merging vertices. Negative means ignored
///                                        (remapping is based on vertex indices only) }
/// @param[out] vertex_id_map              { Map from the concatenated vertex ids to the indices
///                                        after de-duplicatation }
/// @param[out] parent_face                { #OF list of parent face ids }
/// @param[out] boundary_edges             { m x 2 list of undirected edges that lie on the boundary
///                                        of the original mesh }
///
/// @return     { True in case of success }
///
bool instantiate_pattern_aux(
	const Eigen::MatrixXd &IV,
	const Eigen::MatrixXi &IF,
	std::function<const Eigen::MatrixXd *(int)> PV,
	std::function<const Eigen::MatrixXi *(int)> PF,
	std::function<std::array<Eigen::VectorXi, 4>(int)> border_vertices,
	Eigen::MatrixXd &OV,
	Eigen::MatrixXi &OF,
	bool remap_duplicated_vertices = true,
	double tolerance = -1.0,
	Eigen::VectorXi *vertex_id_map = nullptr,
	Eigen::VectorXi *parent_face = nullptr,
	Eigen::MatrixXi *boundary_edges = nullptr
);

// Each quad has its own pattern mesh (PV, PF)
bool instantiate_pattern(
	const Eigen::MatrixXd &IV,
	const Eigen::MatrixXi &IF,
	const std::vector<Eigen::MatrixXd> &PV,
	const std::vector<Eigen::MatrixXi> &PF,
	Eigen::MatrixXd &OV,
	Eigen::MatrixXi &OF,
	bool remap_duplicated_vertices = true,
	double tolerance = -1.0,
	Eigen::VectorXi *vertex_id_map = nullptr,
	Eigen::VectorXi *parent_face = nullptr,
	Eigen::MatrixXi *boundary_edges = nullptr
);

} // namespace micro_quadfoam

#endif
