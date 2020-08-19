////////////////////////////////////////////////////////////////////////////////
// MarchingSquaresStitch.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Implements an axis-aligned grid that can extract the polygon of an
//  object with marching cubes. Instead of marching counterclockwise around
//  the contour, this version generates all the points on a per-cell basis,
//  stitching them together accordingly.
//
//  This version works with objects that may extend outside the grid; in
//  that case, the resulting polygon is actually the intersection of the
//  object with the meshing box.
//
//      16 Cases of the corners:
//      .---.  .---.  .---.  .---.
//      | 0 |  | 1 |  | 2 |  | 3 |
//      '---'  *---'  '---*  *---*
//      .---*  .---*  .---*  .---*
//      | 4 |  | 5 |  | 6 |  | 7 |
//      '---'  *---'  '---*  *---*
//      *---.  *---.  *---.  *---.
//      | 8 |  | 9 |  |10 |  |11 |
//      '---'  *---'  '---*  *---*
//      *---*  *---*  *---*  *---*
//      |12 |  |13 |  |14 |  |15 |
//      '---'  *---'  '---*  *---*
//
//  The following node ordering convention is used:
//      3--2
//      |  |
//      0--1
//  (As implemented in Grid2D)
//
//  This version does *not* use a lookup table. Instead, it notes that all cases
//  are a rotation of the following:
//      .---.  .---.  .---.  .-u-*  *---*
//      | 0 |  v 1 |  v 3 u  v 7 |  |15 |
//      '---'  *-u-'  *---*  *---*  *---*
//  and the ambiguous case:
//      .---*
//      | 5 |
//      *---'
//  Points are generated for each cell edge with mixed endpoints. All
//  unambiguous cases have exactly two points. The orientation of the new
//  segment is also clear: it should point from the point *after* an
//  inside-object corner (in counter-clockwise order) to the point after an
//  outside-object corner. In other words, it points from the "falling edge"
//  (of the domain indicator function) vertex to the "rising edge" vertex. The
//  new points are marked as u and v above, and the segment would be u->v.
//  In the ambiguous case, four points are created. Sampling the center
//  determines which points must connect (u->v, p->q)
//      .---*       .-p-*    .-p-*
//      | 5 |  ==>  q * v    v O q
//      *---'       *-u-'    *-u-'
//                   (a)      (b)
//  a) 1st falling to 1st rising, 2nd falling to 2nd rising
//  b) 1st falling to 2nd rising, 2nd falling to 1st rising
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/03/2015 15:46:48
////////////////////////////////////////////////////////////////////////////////
#ifndef MARCHINGSQUARESSTITCH_HH
#define MARCHINGSQUARESSTITCH_HH
#include <vector>
#include <map>
#include "Grid.hh"
#include "AdaptiveEvaluator.hh"
#include <algorithm>
#include <cassert>
#include <functional>

#include <string>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/GlobalBenchmark.hh>

#include <queue>

class MarchingSquaresGrid : public Grid2D {
public:
    MarchingSquaresGrid(size_t Nx, size_t Ny, size_t clevels = 0) : Grid2D(Nx, Ny, BBox<Vector2D>()) {
        setCoarsening(clevels);
    }

    typedef std::pair<size_t, size_t> Edge;
    enum class SegmentType { Interior, Border };
    typedef std::pair<SegmentType, std::vector<Edge>> Segment;
    // CCW ordering of CCW boundary segments.
    // Distinguishes between:
    //      interior segment: one or both endpoints are non-border points
    //      border   segment: both endpoints are border points
    struct MarchingSquaresResult {
        // Separate edges into CCW ordered segments of border/interior edges.
        // Assumptions (which hold for marching squares edge soup):
        // 1)  Input is manifold
        // 2)  All polygons are closed
        MarchingSquaresResult(std::vector<Vector2d> &&p, std::vector<Edge> &&edges,
                              const std::vector<int> &borderMarkers);

        size_t numEdges() const {
            size_t count = 0;
            for (const auto &s : segments) count += s.second.size();
            return count;
        }

        std::vector<Vector2d> points;
        std::vector<Segment>  segments;
    };

    static MeshIO::IOVertex canonicalEmbedding(const Vector2d &p) {
        return MeshIO::IOVertex(p);
    }

    using EmbeddingFunction = std::function<MeshIO::IOVertex(const Vector2d &)>;
    template<typename Domain>
    void outputSignedDistanceField(const std::string &path, const Domain &domain,
                                   const EmbeddingFunction &embedder = canonicalEmbedding) {
        m_bbox = domain.boundingBox();
        std::vector<Real> signedDistance(numVertices());
        for (size_t v = 0; v < numVertices(); ++v)
            signedDistance[v] = domain.signedDistance(vertexPosition(v));
        outputSignedDistanceField(path, signedDistance, embedder);
    }

    void outputSignedDistanceField(const std::string &path, std::vector<Real> &sdVec,
                                   const EmbeddingFunction &embedder = canonicalEmbedding) {
        assert(sdVec.size() == numVertices());
        std::vector<MeshIO::IOVertex>  vertices;
        std::vector<MeshIO::IOElement> quads;
        ScalarField<Real> sd(sdVec.size());
        for (size_t i = 0; i < numVertices(); ++i) {
            vertices.emplace_back(embedder(vertexPosition(i)));
            sd[i] = sdVec[i];
        }
        Grid2D::AdjacencyVec corners;
        for (size_t ci = 0; ci < numCells(); ++ci) {
            cellVertices(ci, corners);
            quads.emplace_back(corners[0], corners[1], corners[2], corners[3]);
        }
        MSHFieldWriter writer(path, vertices, quads, MeshIO::MESH_QUAD);
        writer.addField("signed distance", sd);
    }

    void outputContours(const std::string &path, std::vector<Vector2d> &points, std::vector<Edge> &edges)
    {
        std::vector<MeshIO::IOVertex>  vertices;
        std::vector<MeshIO::IOElement> elems;
        for (size_t i = 0; i < points.size(); ++i) {
            vertices.emplace_back(canonicalEmbedding(points[i]));
        }
        for (size_t e = 0; e < edges.size(); ++e) {
            elems.emplace_back(edges[e].first, edges[e].second);
        }
        save(path, vertices, elems);
    }

    // Configure the number of coarsening levels used in adaptive sampling.
    // Initially the signed distance will be sampled at a grid 2^levels times
    // coarser than this grid. Then the grid is refined adaptively in regions
    // where the contour is detected.
    void setCoarsening(size_t levels) { m_coarseningLevels = levels; }

    template<typename Domain>
    MarchingSquaresResult extractBoundaryPolygons(
            const Domain &domain, typename Domain::Real /*mergeThreshold*/ = 0.10);
private:
    // Use adapative evaluation: initially sample at a grid
    // 2^m_coarseningLevels times coarser than this grid, then refine around
    // the detected contour.
    size_t m_coarseningLevels;
    template<typename Real>
    size_t m_getLerpPoint(size_t a, size_t b, Real sda, Real sdb,
                          std::vector<Vector2d> &points,
                          std::map<Edge, size_t> &edgePointMap);
    bool m_closeBoundary(std::vector<Vector2d> &points, std::vector<Edge> &edges,
                         std::vector<int> &borderMarkers);
};


template<typename Domain>
MarchingSquaresGrid::MarchingSquaresResult
MarchingSquaresGrid::extractBoundaryPolygons(
        const Domain &domain, typename Domain::Real /*mergeThreshold*/)
{
    typedef typename Domain::Real Real;

    std::vector<Vector2d> points;
    std::vector<Edge> edges;

    BENCHMARK_START_TIMER("Eval signed distance");
    m_bbox = domain.boundingBox();
    AdaptiveEvaluator sdEval(domain, *this, m_coarseningLevels);
    BENCHMARK_STOP_TIMER("Eval signed distance");

#if DEBUG_MARCHING_SQUARES
    std::vector<Real> signedDistance(numVertices());
    for (size_t v = 0; v < numVertices(); ++v)
        signedDistance[v] = domain.signedDistance(vertexPosition(v));
    static size_t _meshRun = 0;
    outputSignedDistanceField("sd_" + std::to_string(_meshRun++) + ".msh", signedDistance);
#endif // DEBUG_MARCHING_SQUARES

    Grid2D::AdjacencyVec corners;
    // Indices of cell edge vertices, in ccw order.
    // Guaranteed to alternate falling-rising, but may start on a rising edge.
    std::vector<size_t> edgeVertices;
    std::map<Edge, size_t> edgePointMap;
    for (size_t ci = 0; ci < numCells(); ++ci) {
        if (!sdEval.cellOverlapsContour(ci)) continue; // Ignore cells away from contour
        edgeVertices.clear();
        cellVertices(ci, corners);
        // Offset in edgeVertices of the cell's first falling-edge vertex.
        size_t firstFalling = 4;
        for (size_t cv = 0; cv < 4; ++cv) {
            size_t a = corners[ cv      % 4],
                   b = corners[(cv + 1) % 4];
            Real sda = sdEval.signedDistance(a), sdb = sdEval.signedDistance(b);
            bool aInside = sda <= 0, bInside = sdb <= 0;
            if (aInside != bInside) {
                // Detect the cell's first falling edge.
                if (aInside && (firstFalling == 4))
                    firstFalling = edgeVertices.size();
                edgeVertices.push_back(m_getLerpPoint(a, b, sda, sdb, points, edgePointMap));
            }
        }
        assert(edgeVertices.size() % 2 == 0);
        assert((edgeVertices.size() == 0) ||
               (firstFalling < edgeVertices.size()));
        if (edgeVertices.size() == 2) {
            // Unambiguous case: falling to rising
            edges.push_back({edgeVertices[firstFalling],
                             edgeVertices[(firstFalling + 1) % 2]});
        }
        if (edgeVertices.size() == 4) {
            // Ambiguous case: sample center
            //      .---*       .-p-*    .-p-*
            //      | 5 |  ==>  q * v    v O q
            //      *---'       *-u-'    *-u-'
            //                   (a)      (b)
            //  Center inside  (a): connect falling to next rising.
            //  Center outisde (b): connect falling to prev rising.
            if (domain.isInside(cellMidpointPosition(ci))) {
                edges.push_back({edgeVertices[ firstFalling         ], edgeVertices[(firstFalling + 1) % 4]});
                edges.push_back({edgeVertices[(firstFalling + 2) % 4], edgeVertices[(firstFalling + 3) % 4]});
            }
            else {
                edges.push_back({edgeVertices[ firstFalling         ], edgeVertices[(firstFalling + 3) % 4]});
                edges.push_back({edgeVertices[(firstFalling + 2) % 4], edgeVertices[(firstFalling + 1) % 4]});
            }
        }
    }

    std::vector<int> borderMarkers;
    bool hadOpenPolygons = m_closeBoundary(points, edges, borderMarkers);
    if (!hadOpenPolygons) {
        // We must handle the case where the entire grid border lies within
        // the object. This happens iff a single grid corner lies inside the
        // object and no open polygons were created by marching squares.
        std::vector<Vector2d> gridCorner = {
            Vector2d(m_bbox.minCorner[0], m_bbox.minCorner[1]),
            Vector2d(m_bbox.maxCorner[0], m_bbox.minCorner[1]),
            Vector2d(m_bbox.maxCorner[0], m_bbox.maxCorner[1]),
            Vector2d(m_bbox.minCorner[0], m_bbox.maxCorner[1])
        };
        if (domain.signedDistance(gridCorner[0]) < 0) {
            std::vector<size_t> newPoints;
            for (const auto &c : gridCorner) {
                // All corners better be inside if any one is...
                assert(domain.signedDistance(c) < 0);
                newPoints.push_back(points.size());
                points.push_back(c);
            }
            edges.push_back({newPoints[0], newPoints[1]});
            edges.push_back({newPoints[1], newPoints[2]});
            edges.push_back({newPoints[2], newPoints[3]});
            edges.push_back({newPoints[3], newPoints[0]});
            borderMarkers.push_back(1);
            borderMarkers.push_back(2);
            borderMarkers.push_back(3);
            borderMarkers.push_back(4);
        }
    }

    return MarchingSquaresResult(std::move(points), std::move(edges), borderMarkers);
}

// Get the point linearly interpolating between corner vertices a and b based on
// the corresponding signed distances sda and sdb.
// This method assumes sign(sda) != sign(sdb).
// A point is only created once per (undirected) cell edge (a, b), so vertices
// will be merged appropriately.
template<typename Real>
size_t MarchingSquaresGrid::
m_getLerpPoint(size_t a, size_t b, Real sda, Real sdb,
               std::vector<Vector2d> &points,
               std::map<Edge, size_t> &edgePointMap)
{
    // Relabel so that sda < sdb. This both implements unordered cell edge
    // lookup in edgePointMap and ensures that we get the *identical* floating
    // point result on the shared edge of two grids on adjacent faces of an
    // axis-aligned cube. The problem this solves is:
    //      Swapping a and b replaces "alpha" with "1 - alpha" below. While
    //      this represents an identical interpolated point in exact arithmetic,
    //      roundoff error in alpha/1-alpha results in a different point.
    if (sda > sdb) {
        std::swap(sda, sdb);
        std::swap(a, b);
    }

    // Look up unordered cell edges...
    auto key = std::make_pair(a, b);
    auto it = edgePointMap.find(key);
    if (it != edgePointMap.end()) return it->second;

    auto va = vertexPosition(a);
    auto vb = vertexPosition(b);

    assert(sda != sdb);
    Vector2d newPt;
    // Explicitly handle cases where grid points lie on the contour.
    bool needsLerp = true;
    if (sda == 0) { newPt = va; needsLerp = false; }
    if (sdb == 0) { newPt = vb; needsLerp = false; }

    if (needsLerp) {
        assert(sda * sdb < 0);
        // Linear approximation of zero crossing...
        // f(x) = [(x - a) / (a - b)] * sdb + [1 - (x - a) / (a - b)] * sda
        // f(x) = 0 ==> x = (sdb * a - sda * b) / (sdb - sda)
        //                = alpha * a + (1 - alpha) * b
        // where alpha = sdb / (sdb - sda)
        // Assuming sdb and sda differ in sign, 0 <= alpha <= 1.
        Real alpha = sdb / (sdb - sda);
        newPt = alpha  * va + (1 - alpha) * vb;

        // Because neither gridpoint is on the contour (sda, sdb != 0), if
        // "newPoint" coincides with a grid point, there was a roundoff error
        // in the linear interpolation that has caused a topology error.
        // Correct this topology error by constructing the closest possible
        // double-precision-quantized point to va on segment (va, vb).
        bool needsPerturbation = false;
        Vector2d limitPoint, candidatePoint;
        if (newPt == va) { limitPoint = va; candidatePoint = (va + vb) / 2; needsPerturbation = true; }
        if (newPt == vb) { limitPoint = vb; candidatePoint = (va + vb) / 2; needsPerturbation = true; }
        if (needsPerturbation) {
            // Approach the limit point along (va, vb), stopping just before reaching it.
            while (true) {
                Vector2d newCandidate = (limitPoint + candidatePoint) / 2;
                if (newCandidate == limitPoint) break;
                candidatePoint = newCandidate;
            }
            newPt = candidatePoint;
        }

    }
    // else {
    //     std::cout << " exact gridpoint crossing at " << va.transpose() << std::endl;
    // }

    size_t newPointIdx = points.size();
    edgePointMap.emplace(key, newPointIdx);
    points.emplace_back(newPt);

    return newPointIdx;
}

// If the object extends outside the marching squares grid, there might be open
// curves incident on the grid border. Close these by inserting edges along the
// grid border. Note that we may need to add new vertices at the grid corners if
// the object overlaps them.
//
// Note that even if there are no such open curves incident on the grid border,
// the object may still extend outside the grid: the entire grid border may lie
// inside the object. This case is handled by the marching squares algorithm
// after we return, since it requires access to the signed distance function
// (and really isn't boundary closing...)
//
// The orientation of the incident curve tells us which direction the new edge
// should head (out or in). We consider three types of points:
//      source intersection vertex: new edge will leave these
//      sink   intersection vertex: new edge will enter these
//      grid corner:                0 or 2 (entering, leaving) new edges
// The grid corner is only created if new incident edges are needed.
// Algorithm:
//      Sort all curve-border intersection points corners in counter-clockwise
//      order, inserting grid corners.
//      Discard grid corners if they are unneeded
//          needed only if a source comes before and a sink comes after
//      Verify that vertices now alternate (source, sink), ignoring grid
//      corners.
//      starting at a source vertex, insert the source->(grid point)->sink edges
// @param[in]    points         points created by marching squares
// @param[inout] edges          marching squares edges (no border edges yet)
//                              Border edges are appended by this method.
// @param[out]   borderMarkers  Zero for all original (interior) edges,
//                              nonzero for border edges. Encodes which border
//                              the edge lies on.
//                              WARNING: these markers are only able to indicate
//                              whether or not two edges in a segment are on the
//                              same border; the grid borders' markers can be
//                              any cyclic permutation of {1, 2, 3, 4}.
// @return      didCloseBorder  report if there were actually any open polygons
//                              incident on the border... if not, it's possible
//                              the entire grid border is within the object,
//                              which must be handled by the caller.
inline
bool MarchingSquaresGrid::m_closeBoundary(std::vector<Vector2d> &points,
                                          std::vector<Edge> &edges,
                                          std::vector<int> &borderMarkers)
{
    borderMarkers.assign(edges.size(), 0);

    std::vector<int> flux(points.size(), 0);
    // -1 for sink, +1 for source
    for (const auto &e : edges) {
        --flux[e.first];
        ++flux[e.second];
    }

#if DEBUG_MARCHING_SQUARES
    outputContours("mc_contours.obj", points, edges);
#endif

    // Point sets for 1d ordering: bottom, right, top, left
    std::vector<std::vector<std::pair<size_t, Real>>> edgeVertices(4);
    for (size_t i = 0; i < points.size(); ++i) {
        const auto &p = points[i];
        if (flux[i] != 0) {
            // Use a small tolerance to check that a point lies on the bounding box
            // When the bbox coordinates are not a power of 2, the lerp interpolation
            // can introduce roundoff errors even when interpolating between the
            // same floating point. More specifically, due to roundoff errors,
            // `alpha * x + (1 - alpha) * x` is not always equal to `x` (unless
            // the code is compiled with -ffast-math). For coordinates between
            // [-10^3, 10^3], a quick stochastic experiment shows that the error
            // stays below 10^-12, which is what we use here.
            const double tol = 1e-12;
            bool onBottom = std::abs(p[1] - m_bbox.minCorner[1]) < tol;
            bool onRight  = std::abs(p[0] - m_bbox.maxCorner[0]) < tol;
            bool onTop    = std::abs(p[1] - m_bbox.maxCorner[1]) < tol;
            bool onLeft   = std::abs(p[0] - m_bbox.minCorner[0]) < tol;
            assert((onBottom + onRight + onTop + onLeft == 1) &&
                   "Valence 1 vertices must lie on exactly one grid border");

            if      (onBottom) edgeVertices[0].push_back({i, p[0]}); // bottom
            else if (onRight ) edgeVertices[1].push_back({i, p[1]}); // right
            else if (onTop   ) edgeVertices[2].push_back({i, p[0]}); // top
            else if (onLeft  ) edgeVertices[3].push_back({i, p[1]}); // left
        }
    }
    auto  ascending = [](const std::pair<size_t, Real> &a, const std::pair<size_t, Real> &b) { return a.second < b.second; };
    auto descending = [](const std::pair<size_t, Real> &a, const std::pair<size_t, Real> &b) { return a.second > b.second; };
    std::sort(edgeVertices[0].begin(), edgeVertices[0].end(),  ascending);
    std::sort(edgeVertices[1].begin(), edgeVertices[1].end(),  ascending);
    std::sort(edgeVertices[2].begin(), edgeVertices[2].end(), descending);
    std::sort(edgeVertices[3].begin(), edgeVertices[3].end(), descending);

    // If all polygons are entirely inside the grid, there is nothing to close.
    size_t numEndpoints = edgeVertices[0].size() +
                          edgeVertices[1].size() +
                          edgeVertices[2].size() +
                          edgeVertices[3].size();
    if (numEndpoints == 0) return false;

    // Type and index of each border vertex
    enum class BoundaryVertexType { Source, Sink, GridCorner };
    typedef std::pair<BoundaryVertexType, size_t> BVertex;
    std::vector<BVertex> sortedVertices;
    // (GridCorner vertices get special indices:
    //      3--2
    //      |  |
    //      0--1)
    for (size_t ei = 0; ei < 4; ++ei) {
        sortedVertices.push_back({BoundaryVertexType::GridCorner, ei});
        const auto &edge = edgeVertices[ei];
        for (size_t vi = 0; vi < edge.size(); ++vi) {
            size_t vidx = edge[vi].first;
            auto type = flux[vidx] == -1 ? BoundaryVertexType::Sink
                                         : BoundaryVertexType::Source;
            sortedVertices.push_back({type, vidx});
        }
    }

    // We want to start traversal at the first source vertex.
    // We are guaranteed at least one source and sink in the MS result.
    auto it = std::find_if(sortedVertices.begin(), sortedVertices.end(),
        [](const BVertex &v) { return v.first == BoundaryVertexType::Source; });
    assert(it != sortedVertices.end());
    size_t startOffset = std::distance(sortedVertices.begin(), it);
    size_t numBorderVertices = sortedVertices.size();
    auto wrappedIndex = [=](size_t i) { return i % numBorderVertices; };
    auto advance = [=](size_t &curr, size_t &next, size_t &i) { curr = next; next = wrappedIndex(next + 1); ++i; };
    std::vector<Vector2d> gridCorner = {
        Vector2d(m_bbox.minCorner[0], m_bbox.minCorner[1]),
        Vector2d(m_bbox.maxCorner[0], m_bbox.minCorner[1]),
        Vector2d(m_bbox.maxCorner[0], m_bbox.maxCorner[1]),
        Vector2d(m_bbox.minCorner[0], m_bbox.maxCorner[1])
    };
    // Get a vertex's index in "points," inserting the grid corner vertices on
    // demand. Note: once a grid corner vertex is inserted, it becomes a source
    // vertex.
    auto vtxIndex = [&](size_t svidx) {
        BVertex &v = sortedVertices[svidx];
        // Convert GridCorner vertices into sources by creating them.
        if (v.first == BoundaryVertexType::GridCorner) {
            size_t idx = points.size();
            points.push_back(gridCorner.at(v.second));
            v.first = BoundaryVertexType::Source;
            v.second = idx;
        }
        return v.second;
    };
    // Traverse border vertices in counterclockwise order, forming the border
    // edges by joining sources to sinks.
    // During traversal, we must always alternate between source and sink
    // vertices. We discard grid corners if they come after a sink vertex and
    // keep them if they come after a source vertex.
    // Note: grid corners are delimiters for the borders. Since we start at the
    //       *first* source vertex, we can assign it border marker 1 and
    //       increment ever time we see a grid corner.
    int borderMarker = 1;
    auto addBorderEdge = [&](size_t u, size_t v) {
        borderMarkers.push_back(borderMarker);
        edges.push_back({vtxIndex(u), vtxIndex(v)});
    };
    size_t u = startOffset, v = wrappedIndex(startOffset + 1);
    for (size_t i = 0; i < numBorderVertices; /* advanced inside */ ) {
        assert(sortedVertices[u].first == BoundaryVertexType::Source);
        // Add any grid corners between source and sink
        while (sortedVertices[v].first == BoundaryVertexType::GridCorner) {
            addBorderEdge(u, v);
            advance(u, v, i); // u = new source, v = next grid corner or sink
            ++borderMarker;
        }
        assert(sortedVertices[v].first == BoundaryVertexType::Sink);
        addBorderEdge(u, v);
        advance(u, v, i); // u = sink, v = next vertex
        // Skip any grid corner vertices following a sink
        while (sortedVertices[v].first == BoundaryVertexType::GridCorner) {
            advance(u, v, i); // v = next grid corner or source
            ++borderMarker;
        }
        advance(u, v, i); // u = next source
    }
    // Report that there were actually open polygons that we closed.
    return true;
}

inline
void mergeEdges(const std::vector<std::pair<size_t, size_t>> &edgesToMerge,
                std::vector<Vector2d> &p,
                std::vector<std::pair<size_t, size_t>> &edges) {
    // "edgesToMerge" represents a graph connecting vertices that should be
    // merged into a single vertex.
    // We extract the connected components of this graph
    const size_t numPts = p.size();
    std::vector<std::vector<size_t>> adj(numPts);
    for (const auto &e : edgesToMerge) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }

    const size_t none = std::numeric_limits<size_t>::max();
    std::vector<size_t> mergedIdxForPt(numPts, none);
    auto visited = [&](size_t i) { return mergedIdxForPt.at(i) != none; };
    std::queue<size_t> bfsQ;
    size_t mergedIdx = 0;
    for (size_t i = 0; i < numPts; ++i) {
        if (visited(i)) continue;
        mergedIdxForPt[i] = mergedIdx;
        bfsQ.push(i);
        while (!bfsQ.empty()) {
            size_t u = bfsQ.front(); bfsQ.pop();
            for (size_t v : adj[u]) {
                if (visited(v)) continue;
                mergedIdxForPt[v] = mergedIdx;
                bfsQ.push(v);
            }
        }
        ++mergedIdx;
    }

    std::vector<Vector2d> newPoints(mergedIdx);
    std::vector<bool> isSet(mergedIdx, false);
    for (size_t i = 0; i < numPts; ++i) {
        size_t mi = mergedIdxForPt[i];
        if (!isSet.at(mi)) {
            newPoints[mi] = p[i];
            isSet[mi] = true;
        }
        else assert(newPoints[mi] == p[i]);
    }
    for (bool b : isSet) assert(b && "Error: merged point missing.");

    p.swap(newPoints);

    // Re-link edges to the new points.
    for (auto &e : edges) {
        e.first  = mergedIdxForPt[e.first];
        e.second = mergedIdxForPt[e.second];
    }

    std::cerr << "WARNING: merged " << edgesToMerge.size() << " degenerate edges from marching squares output." << std::endl;
}

// Separate edges into CCW ordered segments of border/interior edges.
// Assumptions (which hold for marching squares edge soup):
// 1)  Input is manifold
// 2)  All polygons are closed
// A new segment is also started upon transitioning between grid borders
// (as indicated by the borderMarkers array).
inline
MarchingSquaresGrid::MarchingSquaresResult::MarchingSquaresResult(
    std::vector<Vector2d> &&p, std::vector<Edge> &&edges,
    const std::vector<int> &borderMarkers)
{
    size_t numEdges = edges.size();
    assert(numEdges == borderMarkers.size());
    assert(p.size() == numEdges);
    points.swap(p);

    // Mark and merge degenerate edges.
    std::vector<bool> isDegenerate(edges.size(), false);
    std::vector<Edge> edgesToMerge;
    for (size_t ei = 0; ei < edges.size(); ++ei) {
        const auto &e = edges[ei];
        if (points[e.first] == points[e.second]) {
            edgesToMerge.push_back(e);
            isDegenerate[ei] = true;
        }
    }

    if (edgesToMerge.size() > 0)
        mergeEdges(edgesToMerge, points, edges);
    size_t numPoints = points.size();

    const size_t none = std::numeric_limits<size_t>::max();
    std::vector<size_t>  exitingEdge(numPoints, none),
                        enteringEdge(numPoints, none); // Debugging

    // Fill out CCW edge adjacency
    std::vector<size_t> nextEdge(numEdges, none);
    {
        for (size_t ei = 0; ei < numEdges; ++ei) {
            if (isDegenerate[ei]) continue;
            const auto &e = edges[ei];
            assert(( exitingEdge[ e.first] == none) && "Error: nonmanifold point");
            assert((enteringEdge[e.second] == none) && "Error: nonmanifold point");
             exitingEdge.at( e.first) = ei;
            enteringEdge.at(e.second) = ei; // Debugging
        }
        for (size_t ee :  exitingEdge) assert((ee != none) && "Error: Dangling vertex (no exiting edge)");
        for (size_t ee : enteringEdge) assert((ee != none) && "Error: Dangling vertex (no entering edge)");

        for (size_t ei = 0; ei < numEdges; ++ei) {
            if (isDegenerate[ei])
                continue;
            nextEdge[ei] = exitingEdge.at(edges[ei].second);
        }
    }

    // Visit all polygons
    std::vector<bool> visited(numEdges, false);
    size_t visitedCount = 0;
    auto visit = [&](size_t i) { assert(!visited[i]); visited[i] = true; ++visitedCount; };
    for (size_t start = 0; start < numEdges; ++start) {
        if (visited[start] || isDegenerate[start]) continue;

        // start, nextEdge[start], ... is now a cycle of unvisited edges. Break
        // it into interior/border segments:
        // If the entire polygon consists of only interior or border edges, we
        // can start the segment anywhere. Otherwise, to ensure a full segment,
        // we must start on the first of a new run of interior or border edges.
        size_t curr = start;
        while (borderMarkers[curr] == borderMarkers[nextEdge[curr]]) {
            assert(nextEdge[curr] < numEdges);
            assert(!isDegenerate[nextEdge[curr]]);
            curr = nextEdge[curr];
            if (curr == start) break;
        }
        // Start segment on the first edge of the next run...
        // (Or anywhere if all edges are of same type).
        assert(nextEdge[curr] < numEdges);
        assert(!isDegenerate[nextEdge[curr]]);
        curr = nextEdge[curr];
        int segmentBorderMarker = borderMarkers[curr];

        std::vector<Edge> segmentEdges;
        while (!visited[curr]) {
            segmentEdges.push_back(edges[curr]);
            visit(curr);
            assert(nextEdge[curr] < numEdges);
            assert(!isDegenerate[nextEdge[curr]]);
            curr = nextEdge[curr];
            // Commit segment when we detect its end.
            if (visited[curr] || (borderMarkers[curr] != segmentBorderMarker)) {
                auto type = (segmentBorderMarker == 0) ? SegmentType::Interior
                                                       : SegmentType::Border;
                segments.push_back({type, std::vector<Edge>()});
                segments.back().second.swap(segmentEdges);
                segmentBorderMarker = borderMarkers[curr];
            }
        }
    }

    assert(visitedCount + edgesToMerge.size() == numEdges);
}

#endif /* end of include guard: MARCHINGSQUARESSTITCH_HH */
