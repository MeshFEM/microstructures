#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <nanoflann.hpp>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace micro {

template<typename Real, int Dim>
struct PointCloud
{
	using Point = Eigen::Matrix<Real, Dim, 1>;
	using TaggedPoint = std::pair<Point, size_t>;

	std::vector<TaggedPoint> pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the dim'th component of the idx'th point in the class
	inline Real kdtree_get_pt(const size_t idx, const size_t dim) const { return pts[idx].first[dim]; }

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

////////////////////////////////////////////////////////////////////////////////

template<typename Real, int Dim>
class DynamicKdTree {
public:

	using KdTreeType = nanoflann::KDTreeSingleIndexDynamicAdaptor<
		nanoflann::L2_Simple_Adaptor<Real, PointCloud<Real, Dim> >,
		PointCloud<Real, Dim>,
		Dim>;

	using Point = typename PointCloud<Real, Dim>::Point;

public:
	DynamicKdTree();

	// No copy/move allows
    DynamicKdTree(DynamicKdTree&&) = delete;
    DynamicKdTree& operator=(DynamicKdTree&&) = delete;
    DynamicKdTree(const DynamicKdTree&) = delete;
    DynamicKdTree& operator=(const DynamicKdTree&) = delete;

	std::pair<int, Point> getClosestPoint(const Point &p, double eps) const;

	void addPoint(const Point &p, size_t tag = 0);

protected:
	PointCloud<Real, Dim> m_PointCloud;
	KdTreeType m_Tree;
};

////////////////////////////////////////////////////////////////////////////////

using DynamicKdTree2d = DynamicKdTree<double, 2>;
using DynamicKdTree3d = DynamicKdTree<double, 3>;

} // namespace micro
