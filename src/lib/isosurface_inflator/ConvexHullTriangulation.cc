#include "ConvexHullTriangulation.hh"

#include "DisableWarnings.hh"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>
#include "EnableWarnings.hh"

#include <iterator>

template<class PointCollection>
void convexHullFromTriangulation(const PointCollection &points,
                std::vector<MeshIO::IOVertex > &hullVertices,
                std::vector<MeshIO::IOElement> &hullElements,
                std::vector<size_t>            &originatingVertexIndices)
{
    using K             = CGAL::Exact_predicates_inexact_constructions_kernel;
    using PointInfo     = size_t;
    using VertexBase    = CGAL::Triangulation_vertex_base_with_info_3<PointInfo, K>;
    using TDS           = CGAL::Triangulation_data_structure_3<VertexBase>;
    using Triangulation = CGAL::Triangulation_3<K, TDS>;
    using Point_3       = Triangulation::Point;
    using Vertex_handle = Triangulation::Vertex_handle;
    using Cell_handle   = Triangulation::Cell_handle;

    Triangulation T;
    // Annoyingly, we have to insert points one at a time, because
    // Triangulation_3 doesn't provide a constructor taking a pair of vertices
    // and their attached info. Delaunay_triangulation_3 does, but is probably
    // slower.
    // We must do a spatial sort as in CGAL's Triangulation_3::insert for ranges
    // to make this efficient (otherwise it is over 40x slower than range
    // insert...)
    {
        using IndexedPt = std::pair<Point_3, size_t>;
        using SortTraits = CGAL::Spatial_sort_traits_adapter_3<K,
              CGAL::First_of_pair_property_map<IndexedPt>>;
        std::vector<IndexedPt> cgal_pts;
        cgal_pts.reserve(points.size());
        size_t i = 0;
        for (const auto &pt : points)
            cgal_pts.emplace_back(Point_3(pt[0], pt[1], pt[2]), i++);
        CGAL::spatial_sort(cgal_pts.begin(), cgal_pts.end(), SortTraits());
        Vertex_handle vh;
        for (const auto &ipt : cgal_pts) {
            vh = T.insert(ipt.first, vh); // use prev as hint to accelerate
            vh->info() = ipt.second;
        }
    }

    // Determine vertices on the convex hull
    std::vector<Vertex_handle> hull_vertex_handles;
    hull_vertex_handles.reserve(T.number_of_vertices());
    T.adjacent_vertices(T.infinite_vertex(),
                        std::back_inserter(hull_vertex_handles));

    // Get the originating vertex index for each hull vertex, and reindex for
    // hull extraction
    originatingVertexIndices.clear();
    originatingVertexIndices.reserve(hull_vertex_handles.size());
    {
        size_t i = 0;
        for (const auto &vh : hull_vertex_handles) {
            originatingVertexIndices.push_back(vh->info());
            vh->info() = i++;
        }
    }

    // The convex hull will consist of faces of the infinite cells opposite the
    // infinite vertex.
    std::vector<Cell_handle> hull_incident_cells;
    hull_incident_cells.reserve(T.number_of_cells());
    T.incident_cells(T.infinite_vertex(),
                     std::back_inserter(hull_incident_cells));

    // Exactly one of the cell vertices should be infinite; the rest should be
    // finite and form a triangle of the convex hull. Extract these faces in
    // MeshIO format.
    hullElements.clear();
    hullElements.reserve(hull_incident_cells.size());
    for (const Cell_handle &ch : hull_incident_cells) {
        const size_t NO_IVERTEX = 4;
        size_t infiniteVtx = NO_IVERTEX;
        MeshIO::IOElement e;
        for (size_t i = 0; i < 4; ++i) {
            const Vertex_handle &vh = ch->vertex(i);
            if (T.is_infinite(vh)) { assert(infiniteVtx == NO_IVERTEX); infiniteVtx = i; continue; }
            e.push_back(vh->info());
        }
        assert(infiniteVtx != NO_IVERTEX);
        // Reorient all faces to point outward. Note that CGAL indexes the cell
        // corners as follows:
        //       3
        //       *
        //      / \`
        //     /   \ `* 2
        //    / __--\ /
        //  0*-------* 1
        // Since the cell in question is an infinite cell *outside* the object,
        // we want the folllowing cell-inward orientations:
        //   3      2      1      0      <== opposite (infinite) vertex
        // 0-1-2, 0-3-1, 0-2-3, 1-3-2
        //        flip          flip
        // The two configurations marked "flip" indicate the ones where the
        // current, sorted, ordering needs to be permuted to correct the
        // orientation.
        if ((infiniteVtx == 0) || (infiniteVtx == 2))
            std::swap(e[0], e[1]);

        hullElements.emplace_back(e);
    }

    hullVertices.clear();
    hullVertices.reserve(hull_vertex_handles.size());
    for (const Vertex_handle &vh : hull_vertex_handles) {
        const auto &cgal_pt = vh->point();
        hullVertices.emplace_back(cgal_pt.x(), cgal_pt.y(), cgal_pt.z());
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
