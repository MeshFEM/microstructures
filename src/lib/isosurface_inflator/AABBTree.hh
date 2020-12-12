#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace micro {

// AABB tree for an edge graph
// Leaves corresponds to rounded cylinders
struct AABBTree {

	struct Node {
		Eigen::AlignedBox3d bbox;
		int parent; // Index of the parent node (-1 for root)
		int left; // Index of the left child (-1 for a leaf)
		int right; // Index of the right child (-1 for a leaf)
		int index; // Edge id for the leaf (-1 for internal nodes)

		bool isLeaf() const { return left < 0; }
	};

private:
	std::vector<Node> m_Nodes;
	int m_Root;

public:
	AABBTree() = default; // Default empty constructor

	///
	/// Construct an AABB over the given edge graph. Each edge is interpreted
	/// as a cone frustum with hemispherical endpoints. The radius of each
	/// endpoint (vertex) being defined by the user.
	///
	/// @param[in]  V     #V x 3 input vertex positions
	/// @param[in]  E     #E x 2 input edge vertices
	/// @param[in]  R     #V x 1 input vertex radii
	///
	AABBTree(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::VectorXd &R);

	///
	/// Compute the list of cylinders whose AABB tree overlaps with the query point
	///
	/// @param[in]  p          3 x 1 query point
	/// @param[out] cylinders  List of edge ids of the cylinder overlapping with the query point
	///
	void intersects(const Eigen::RowVector3d &p, std::vector<int> &cylinders) const;

	// Test whether the tree is empty
	bool empty() const { return m_Nodes.empty(); }
};

} // namespace micro
