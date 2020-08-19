////////////////////////////////////////////////////////////////////////////////
// BoxIntersection1DFeatures.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      To resolve the sharp edges created the boolean intersection of a smooth
//      object with a box, CGAL needs an explicit representation of the
//      intersection's sharp feature curves in its "polylines" format (list of
//      point lists). CGAL appears to constrain the output to contain the exact
//      endpoints of these polylines (but the interior of the polyline can be
//      deformed/merged).
//
//      To create this list of feature curves, we run marching squares on the
//      smooth signed distance function restricted to each face of the box.
//      The feature curves lying on the box's edges will be present in two
//      face's marching squares result, so we must be careful to resolve this
//      duplication.
//
//      CGAL is very particular about the positions of the polyline vertices.
//      If the polyline endpoints that are meant to coincide are not perfectly
//      equal, meshing will fail (CGAL tries to keep both points distinct in
//      the output mesh). This means that the marching squares implementation
//      must guarantee exactly consistent coordinates for points on the edges
//      of the box (when created for the incident faces).
//
//      Note: the boolean intersection box is read from the signed distance
//      function's "bounding box" member.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/23/2015 17:03:32
////////////////////////////////////////////////////////////////////////////////
#ifndef BOXINTERSECTION1DFEATURES_HH
#define BOXINTERSECTION1DFEATURES_HH
#include "InflatorTypes.hh"
#include "SignedDistanceRegion.hh"
#include "MarchingSquares/MarchingSquaresStitch.hh"
#include <MeshFEM/MSHFieldWriter.hh>
#include <list>

// Axis-aligned 2D slice of the 3D signed distance function meant for running
// marching squares on the bounding box faces.
// +/- 1: X,    +/- 2: Y,    +/- 2: Z
template<class VolumeSDF>
class BoundaryFaceSlice : public SignedDistanceRegion<2> {
public:
    BoundaryFaceSlice(const VolumeSDF &vsdf, const int faceIdx)
        : m_volRegion(vsdf)
    {
        assert((faceIdx != 0) && (faceIdx >= -3)  && (faceIdx <= 3));
        // Cylic ordering
        m_faceIdx = abs(faceIdx) - 1;
        m_firstIndex  = (m_faceIdx + 1) % 3;
        m_secondIndex = (m_faceIdx + 2) % 3;
        auto bb = vsdf.boundingBox();
        m_faceCoordinate = (faceIdx < 0) ? bb.minCorner[m_faceIdx] : bb.maxCorner[m_faceIdx];
        m_2DBBox = BBox<Point2d>(Point2d(bb.minCorner[m_firstIndex], bb.minCorner[m_secondIndex]),
                                 Point2d(bb.maxCorner[m_firstIndex], bb.maxCorner[m_secondIndex]));
    }

    Point3d volumePoint(const Point2d &facePoint) const {
        Point3d p;
        p[m_faceIdx]     = m_faceCoordinate;
        p[m_firstIndex]  = facePoint[0];
        p[m_secondIndex] = facePoint[1];
        return p;
    }

    virtual const BBox<Point2d> &boundingBox()      const override { return m_2DBBox; }
    virtual double signedDistance(const Point2d &p) const override { return m_volRegion.signedDistance(volumePoint(p)); }

    // There are two copies of each boundary segment: one for each incident box
    // face. We assign unique "ownership" of boundary segments to box faces to
    // avoid this duplication:
    // x face: y = min, max segments \.
    // y face: z = min, max segments  |-- m_firstIndex = min, max (cyclic)
    // z face: x = min, max segments /
    bool ownsSegment(const MarchingSquaresGrid::MarchingSquaresResult &r, size_t si) const {
        const auto &s = r.segments[si];
        if (s.first == MarchingSquaresGrid::SegmentType::Interior) return true;

        // Count how many of the segment's edges are on each border.
        // borders 0 and 1 are m_firstIndex min/max
        // borders 2 and 3 are the other min/max borders
        std::vector<size_t> edgesOnBorder(4, 0);
        size_t numEdges = s.second.size();
        for (const auto &e : s.second) {
            const auto p1 = r.points[e.first];
            const auto p2 = r.points[e.second];
            // Edge is on a border if both of its endpoints are
            edgesOnBorder[0] += (p1[0] == m_2DBBox.minCorner[0]) && (p2[0] == m_2DBBox.minCorner[0]);
            edgesOnBorder[1] += (p1[0] == m_2DBBox.maxCorner[0]) && (p2[0] == m_2DBBox.maxCorner[0]);
            edgesOnBorder[2] += (p1[1] == m_2DBBox.minCorner[1]) && (p2[1] == m_2DBBox.minCorner[1]);
            edgesOnBorder[3] += (p1[1] == m_2DBBox.maxCorner[1]) && (p2[1] == m_2DBBox.maxCorner[1]);
        }

        // Either all edges or no edges should be on a single border.
        // Also, border membership should be disjoint (could happen with
        // degeneracies.
        size_t total = 0;
        for (size_t count : edgesOnBorder) {
            assert((count == 0) || (count == numEdges));
            total += count;
        }
        // All edges of a border segment should actually be on a border, and
        // there should be no double-counting.
        if (total != numEdges) {
            std::cerr << "Invalid border edge memebership:";
            for (size_t count : edgesOnBorder) std::cerr << " " << count;
            std::cerr << std::endl << "segment: " << si << std::endl;
            std::cerr << "edges: " << std::endl;
            for (const auto &e : s.second) {
                const auto p1 = r.points[e.first];
                const auto p2 = r.points[e.second];
                std::cerr << p1[0] << ", " << p1[1] << " -> "
                          << p2[0] << ", " << p2[1] << std::endl;
            }
            assert(total == numEdges); // trigger the assertion
        }
        // Is segment on one of the two border edges we own?
        return (edgesOnBorder[0] + edgesOnBorder[1]);
    }

private:
    const SignedDistanceRegion<3> &m_volRegion;
    BBox<Point2d> m_2DBBox;
    size_t m_faceIdx, m_firstIndex, m_secondIndex;
    double m_faceCoordinate;
};

// Extract the (de-duplicated) 1D sharp features in an obj-like edge-soup
// format. Edges are organized into ccw-orderd "segments," which are chains of
// edges that either form a closed polygon or begin and terminate at a box
// edge. Each segment will become a polyline...
//
// TODO: FIGURE THIS OUT
// It might be possible to get a dangling point (grid corner) from
// de-duplication, but haven't thought it through...
// Each segment is a list of edges indexing into points.
template<typename VolumeSDF, typename _Point>
void boxIntersection1DFeatures(const VolumeSDF &sdf,
        size_t gridSize, size_t gridCoarsening,
        std::vector<_Point> &points,
        std::vector<std::vector<std::pair<size_t, size_t>>> &segmentEdges) {
    MarchingSquaresGrid msquares(gridSize, gridSize, gridCoarsening);
    for (int f = -3; f <= 3; ++f) {
        if (f == 0) continue;
        BoundaryFaceSlice<VolumeSDF> faceSlice(sdf, f);
        std::vector<std::pair<size_t, size_t>> newEdges;
        auto result = msquares.extractBoundaryPolygons(faceSlice, 0.0);

#if SDF_DEBUG_OUT
        {
            auto embedder = [&](const Point2<Real> &p) -> MeshIO::IOVertex { return MeshIO::IOVertex(faceSlice.volumePoint(p)); };
            msquares.outputSignedDistanceField("sdFace_" + std::to_string(f) + ".msh", faceSlice, embedder);
        }
#endif

        size_t offset = points.size();
        points.reserve(points.size() + result.points.size());
        for (const auto &p : result.points) points.push_back(faceSlice.volumePoint(p));

        // De-duplication: only copy over segments that the current face owns
        for (size_t si = 0; si < result.segments.size(); ++si) {
            if (!faceSlice.ownsSegment(result, si)) continue;

            segmentEdges.push_back(std::vector<std::pair<size_t, size_t>>());
            auto &s = segmentEdges.back();
            s.reserve(result.segments[si].second.size());
            for (const auto &e : result.segments[si].second)
                s.push_back({e.first + offset, e.second + offset});
        }
    }
}

// Extract the (de-duplicated) 1D sharp features. Same as above, but this time
// in a CGAL-like polylines format (though the points themselves may not be
// CGAL points).
template<typename VolumeSDF, typename _Point>
void boxIntersection1DFeatures(const VolumeSDF &sdf, size_t gridSize, size_t gridCoarsening,
        std::list<std::vector<_Point>> &polylines)
{
    std::vector<_Point> points;
    std::vector<std::vector<std::pair<size_t, size_t>>> segmentEdges;
    boxIntersection1DFeatures(sdf, gridSize, gridCoarsening, points, segmentEdges);
    polylines.clear();

    // Organize ccw-ordered edge segments into polyline vectors.
    // (CGAL's feature format).
    for (const auto &s : segmentEdges) {
        assert(s.size());
        polylines.emplace_back();
        auto &line = polylines.back();
        line.reserve(s.size());
        size_t i = s[0].first; // for continuity validation
        line.push_back(points.at(i));
        for (const auto &e : s) {
            assert(i == e.first);
            line.push_back(points.at(e.second));
            i = e.second;
        }
    }
}

#endif /* end of include guard: BOXINTERSECTION1DFEATURES_HH */
