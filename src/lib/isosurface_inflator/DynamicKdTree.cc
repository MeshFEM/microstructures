#include "DynamicKdTree.hh"

namespace micro {

template<typename Real, int Dim>
DynamicKdTree<Real, Dim>::DynamicKdTree()
	: m_PointCloud()
	, m_Tree(Dim, m_PointCloud)
{ }

template<typename Real, int Dim>
auto DynamicKdTree<Real, Dim>::getClosestPoint(const Point &p, double eps) const -> std::pair<int, Point>
{
	const size_t numResults = 1;
	size_t retIndex;
	Real outDistSqr;
	nanoflann::KNNResultSet<Real> resultSet(numResults);
	resultSet.init(&retIndex, &outDistSqr);
	m_Tree.findNeighbors(resultSet, p.data(), nanoflann::SearchParams());
	if (outDistSqr <= eps * eps) {
		return std::make_pair(m_PointCloud.pts[retIndex].second, m_PointCloud.pts[retIndex].first);
	} else {
		return std::make_pair(-1, Point::Zero());
	}
}

template<typename Real, int Dim>
void DynamicKdTree<Real, Dim>::addPoint(const Point &p, size_t tag) {
	size_t idx = m_PointCloud.pts.size();
	m_PointCloud.pts.emplace_back(p, tag);
	m_Tree.addPoints(idx, idx);
}

template class DynamicKdTree<double, 2>;
template class DynamicKdTree<double, 3>;

} // namespace micro

