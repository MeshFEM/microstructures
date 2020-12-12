#include "ConvexHullTriangulation.hh"

#include <QuickHull.hpp>

template<class PointCollection>
void convexHullFromTriangulation(const PointCollection &points,
                std::vector<MeshIO::IOVertex > &hullVertices,
                std::vector<MeshIO::IOElement> &hullElements,
                std::vector<size_t>            &originatingVertexIndices)
{
    if (points.empty()) {
        hullVertices.clear();
        hullElements.clear();
        originatingVertexIndices.clear();
        return;
    }

    quickhull::QuickHull<Real> engine;
    std::vector<quickhull::Vector3<Real>> pointCloud;
    pointCloud.reserve(points.size());
    for (const auto &p : points) {
        pointCloud.emplace_back(p(0), p(1), p(2));
    }

    auto hull = engine.getConvexHull(pointCloud, false, true);
    const auto& vertexBuffer = hull.getVertexBuffer();
    const auto& indexBuffer = hull.getIndexBuffer();

    // Find out isolated vertices and reindex accordingly
    size_t numInHull = 0;
    const size_t invalidId = std::numeric_limits<size_t>::max();
    std::vector<size_t> newVertexId(points.size(), invalidId);
    originatingVertexIndices.clear();
    hullVertices.clear();
    for (auto v : indexBuffer) {
        if (newVertexId[v] == invalidId) {
            newVertexId[v] = numInHull++;
            originatingVertexIndices.push_back(v);
            hullVertices.emplace_back(vertexBuffer[v].x, vertexBuffer[v].y, vertexBuffer[v].z);
        }
    }

    hullElements.resize(indexBuffer.size() / 3);
    for (size_t i = 0; 3 * i < indexBuffer.size(); ++i) {
        hullElements[i].resize(3);
        for (size_t lv = 0; lv < 3; ++lv) {
            hullElements[i][lv] = newVertexId[indexBuffer[i * lv + lv]];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
template
void convexHullFromTriangulation<std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>>>(
    const std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>> &points,
    std::vector<MeshIO::IOVertex > &hullVertices,
    std::vector<MeshIO::IOElement> &hullElements,
    std::vector<size_t>            &originatingVertexIndices);
