////////////////////////////////////////////////////////////////////////////////
#include "AABBTree.hh"
#include <numeric>
#include <stack>
////////////////////////////////////////////////////////////////////////////////

namespace micro {

// -----------------------------------------------------------------------------

namespace {

// Build an orthogonal frame given a single vector
void orthogonal_frame(const Eigen::RowVector3d &x, Eigen::RowVector3d &y, Eigen::RowVector3d &z) {
	int imin;
	x.array().abs().minCoeff(&imin);
	Eigen::RowVector3d u;
	// std::cout << imin << std::endl;
	for (int i = 0, s = -1; i < 3; ++i) {
		if (i == imin) {
			u[i] = 0;
			// std::cout << "0, " << std::endl;
		} else {
			int j = (i+1)%3;
			if (j == imin) { j = (i+2)%3; }
			u[i] = s * x[j];
			// std::cout << (s < 0 ? "-" : "") << "v" << j << ", " << std::endl;
			s *= -1;
		}
	}
	// std::cout << u << std::endl;
	z = x.cross(u).normalized();
	y = z.cross(x).normalized();
}

// Bounding box of a sphere
Eigen::AlignedBox3d bbox_sphere(const Eigen::RowVector3d &a, double r) {
	return Eigen::AlignedBox3d(a.array() - r, a.array() + r);
}

// Bounding box of a rounder cone frustum
Eigen::AlignedBox3d bbox_cylinder(
	const Eigen::RowVector3d &a, const Eigen::RowVector3d &b, double ra, double rb)
{
	return bbox_sphere(a, ra).extend(bbox_sphere(b, rb));
}

// Sdf to a rounded cone frustum centered at the origin, aligned with the y axis
double sdf_cylinder_canonical(const Eigen::RowVector3d &p, double r1, double r2, double h) {
	typedef Eigen::RowVector2d Vec2;

    Vec2 q(Vec2(p.x(), p.z()).norm(), p.y());

    double b = (r1-r2)/h;
    double a = std::sqrt(1.0-b*b);
    double k = q.dot(Vec2(-b, a));

    if( k < 0.0 ) return q.norm() - r1;
    if( k > a*h ) return (q-Vec2(0.0,h)).norm() - r2;

    return q.dot(Vec2(a,b)) - r1;
}

} // anonymous namespace

// Sdf to a rounded cone frustum at an arbitrary position
double sdf_cylinder(const Eigen::RowVector3d &a, const Eigen::RowVector3d &b,
	double ra, double rb, const Eigen::RowVector3d &p)
{
	typedef Eigen::RowVector3d Vec3;
	Vec3 y = (b - a).normalized();
	Vec3 x, z;
	orthogonal_frame(y, z, x);
	Eigen::Matrix3d M;
	M.col(0) = x.transpose();
	M.col(1) = y.transpose();
	M.col(2) = z.transpose();
	Vec3 q = p - a;
	return sdf_cylinder_canonical((M.inverse() * q.transpose()).transpose(), ra, rb, (b - a).norm());
}

// -----------------------------------------------------------------------------

AABBTree::AABBTree(const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::VectorXd &R) {
	// Compute the centroids of all the edges in the input mesh
	Eigen::MatrixXd centroids(E.rows(), V.cols());
	centroids.setZero();
	for (int i = 0; i < E.rows(); ++i) {
		for (int k = 0; k < E.cols(); ++k) {
			centroids.row(i) += V.row(E(i, k));
		}
		centroids.row(i) /= E.cols();
	}

	// Top-down approach: split each set of primitives into 2 sets of roughly equal size,
	// based on sorting the centroids along one direction or another.
	std::vector<int> cylinder(E.rows());
	std::iota(cylinder.begin(), cylinder.end(), 0);

	std::function<int(int, int, int)> top_down = [&] (int i, int j, int parent) {
		// Scene is empty, so is the aabb tree
		if (j - i == 0) {
			return -1;
		}

		// If there is only 1 cylinder left, then we are at a leaf
		if (j - i == 1) {
			Node node;
			int e = cylinder[i];
			Eigen::RowVector3d a = V.row(E(e, 0)); double ra = R(E(e, 0));
			Eigen::RowVector3d b = V.row(E(e, 1)); double rb = R(E(e, 1));
			node.bbox = bbox_cylinder(a, b, ra, rb);
			node.parent = parent;
			node.left = node.right = -1;
			node.index = e;
			m_Nodes.push_back(node);
			return (int) (m_Nodes.size() - 1);
		}

		// Otherwise, we need to sort centroids along the longest dimension, and split recursively
		Eigen::AlignedBox3d centroid_box;
		for (int k = i; k < j; ++k) {
			Eigen::Vector3d c = centroids.row(cylinder[k]).transpose();
			centroid_box.extend(c);
		}
		Eigen::Vector3d extent = centroid_box.diagonal();
		int longest_dim = 0;
		for (int dim = 1; dim < 3; ++dim) {
			if (extent(dim) > extent(longest_dim)) {
				longest_dim = dim;
			}
		}
		std::sort(cylinder.begin() + i, cylinder.begin() + j, [&] (int f1, int f2) {
			return centroids(f1, longest_dim) < centroids(f2, longest_dim);
		});

		// Then we can create a new internal node
		int current = m_Nodes.size();
		m_Nodes.resize(current + 1);
		int midpoint = (i + j) / 2;
		int left = top_down(i, midpoint, current);
		int right = top_down(midpoint, j, current);
		Node &node = m_Nodes[current];
		node.left = left;
		node.right = right;
		node.parent = parent;
		node.index = -1;
		node.bbox = m_Nodes[node.left].bbox.extend(m_Nodes[node.right].bbox);

		return current;
	};

	m_Root = top_down(0, cylinder.size(), -1);
}

// -----------------------------------------------------------------------------

void AABBTree::intersects(const Eigen::RowVector3d &query, std::vector<int> &cylinders) const {
	cylinders.clear();
	if (empty()) { return; }
	std::stack<int> q;
	q.push(m_Root);
	while (!q.empty()) {
		const AABBTree::Node &node = m_Nodes[q.top()];
		q.pop();
		if (node.bbox.contains(query.transpose())) {
			if (node.isLeaf()) {
				cylinders.push_back(node.index);
			} else {
				// Internal node, we need to test intersection with both children
				q.push(node.left);
				q.push(node.right);
			}
		}
	}
}

} // namespace micro
