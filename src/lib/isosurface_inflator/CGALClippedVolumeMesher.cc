#include "CGALClippedVolumeMesher.hh"

#if MICRO_WITH_TBB
#define CGAL_LINKED_WITH_TBB
#define CGAL_CONCURRENT_MESH_3
#endif

#include "BoxIntersection1DFeatures.hh"
#include <MeshFEM/GlobalBenchmark.hh>
#include "DisableWarnings.hh"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>
#include "EnableWarnings.hh"
#include <vector>

// avoid verbose function and named parameters call
using namespace CGAL::parameters;
namespace SD = SignedDistance;

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Domain
typedef K::FT FT;
typedef K::Point_3 Point;

struct CGALClippedVolumeMesher::
ClippedSignedDistanceFunction {
    ClippedSignedDistanceFunction(const SignedDistanceRegion<3> &sdf)
        : m_sdf(sdf), m_meshingBox(sdf.boundingBox()) { }

    FT operator()(const Point &p) const {
        // Note: CGAL actually thresholds the signed distance to a binary
        // inside/outside and runs dumb bisection, so don't waste time computing
        // signed distances accurately outside the meshing box.
        if (!m_meshingBox.isInside(p[0], p[1], p[2]))
            return 1.0;
        return m_sdf.isInside(Point3<Real>(p[0], p[1], p[2])) ? -1.0 : 1;
    }
private:
    const SignedDistanceRegion<3> &m_sdf;
    SD::Primitives::Box<Real> m_meshingBox;
};

void CGALClippedVolumeMesher::
mesh(const SignedDistanceRegion<3> &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements) const
{
    using Mesh_domain = CGAL::Mesh_domain_with_polyline_features_3<
                            CGAL::Implicit_mesh_domain_3<ClippedSignedDistanceFunction, K>>;

    // Polyline
    typedef std::vector<Point>    Polyline_3;
    typedef std::list<Polyline_3> Polylines;

    // Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
    using Tr = typename CGAL::Mesh_triangulation_3<Mesh_domain, K,
            CGAL::Parallel_tag // Tag to activate parallelism
        >::type;
#else
    using Tr = typename CGAL::Mesh_triangulation_3<Mesh_domain>::type;
#endif
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, typename Mesh_domain::Corner_index,
                                                    typename Mesh_domain::Curve_segment_index> C3t3;
    using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

    // std::cout << "Meshing 1D features" << std::endl;
    ClippedSignedDistanceFunction cgal_sdfunc(sdf);
    std::list<std::vector<MeshIO::IOVertex>> polylinesMeshIO;
    boxIntersection1DFeatures(sdf, meshingOptions.marchingSquaresGridSize,
                              meshingOptions.marchingSquaresCoarsening, polylinesMeshIO);
    // std::cout << "Checking 1D features..." << std::endl;
    // Check for near-intersection of polylines--this should only happen if
    // we failed to stitch up the boundary curves correctly.
    // This brute-force O(n^2) could easily be sped up...
    for (auto it1 = polylinesMeshIO.begin(); it1 != polylinesMeshIO.end(); ++it1) {
        for (auto it2 = polylinesMeshIO.begin(); it2 != polylinesMeshIO.end(); ++it2) {
            // Close points on the same curve are fine.
            if (it1 == it2) continue;
            for (const auto &v1 : *it1) {
                for (const auto &v2 : *it2) {
                    if (((v1.point - v2.point).norm() < 1e-8) &&
                            ((v1.point[0] != v2.point[0]) ||
                             (v1.point[1] != v2.point[1]) ||
                             (v1.point[2] != v2.point[2]))) {
                        throw std::runtime_error("Inexact intersection of polylines");
                    }
                }
            }
        }
    }
    // std::cout << "Done checking 1D features" << std::endl;

#if 0
    Real maxDist = 0;
    for (const auto &pl : polylinesMeshIO) {
        for (auto &v : pl)
            maxDist = std::max(maxDist, std::abs(sdf.signedDistance(v.point)));
    }
    std::cout << "Max feature dist to surface: " << maxDist << std::endl;

    {
        std::vector<MeshIO::IOVertex > vertices;
        std::vector<MeshIO::IOElement> elements;
        for (const auto &pl : polylinesMeshIO) {
            std::cout << "feature line curve of length: " << pl.size() << std::endl;
            assert(pl.size() >= 2);
            size_t offset = vertices.size();
            for (auto &v : pl) vertices.emplace_back(v);
            for (size_t i = offset; i < vertices.size() - 1; ++i)
                elements.emplace_back(i, i + 1);
        }
        MeshIO::save("cgal_features.msh", vertices, elements);
    }
#endif

    // Convert to CGAL's polyline format
    Polylines polylines;
    for (const auto &l : polylinesMeshIO) {
        polylines.push_back(Polyline_3());
        Polyline_3 &line = polylines.back();
        line.reserve(l.size());
        for (size_t i = 0; i < l.size(); ++i)
            line.emplace_back(l[i][0], l[i][1], l[i][2]);
    }

    Point3d c;
    double r;
    sdf.boundingSphere(c, r);
#if 0
    std::cout << "meshing bounding sphere, radius = " << r << ", c = " << c.transpose() << std::endl;
    std::cout << "signed distance at sphere center: " << sdf.signedDistance(c) << std::endl;
    std::cout << "bounding box: " << sdf.boundingBox() << std::endl;
#endif

    Mesh_domain domain(cgal_sdfunc,
            K::Sphere_3(Point(c[0], c[1], c[2]), r * r), meshingOptions.domainErrorBound);
    // std::cout << "Adding features..." << std::endl;
    domain.add_features(polylines.begin(), polylines.end());

    // Mesh criteria
    Mesh_criteria criteria(facet_angle            = meshingOptions.facetAngle,
                           facet_size             = meshingOptions.facetSize,
                           facet_distance         = meshingOptions.facetDistance,
                           cell_radius_edge_ratio = meshingOptions.cellRadiusEdgeRatio,
                           cell_size              = meshingOptions.cellSize,
                           edge_size              = meshingOptions.edgeSize);

    // Mesh generation
    // std::cout << "Making mesh..." << std::endl;
    BENCHMARK_START_TIMER("make_mesh_3");
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
    // CGAL sometimes returns an empty mesh for some patterns due to
    // insufficient initialization:
    // https://github.com/CGAL/cgal/issues/2416
    // For some reason this is considered "not a bug," and Laurent recommended
    // the following workaround:
    if (c3t3.number_of_facets() == 0) {
        CGAL::internal::Mesh_3::init_c3t3(c3t3, domain, criteria, 20);
        refine_mesh_3(c3t3, domain, criteria);
    }
    BENCHMARK_STOP_TIMER("make_mesh_3");

    // Access triangulation directly
    const Tr &tr = c3t3.triangulation();
    vertices.clear(), elements.clear();
    vertices.reserve(tr.number_of_vertices());
    std::map<typename C3t3::Vertex_handle, size_t> V;
    size_t i = 0;
    for (auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it) {
        V[it] = i++;
        auto p = it->point();
        vertices.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()),
                CGAL::to_double(p.z()));
    }

    elements.reserve(c3t3.number_of_cells_in_complex());
    if (c3t3.number_of_cells_in_complex() == 0) {
        std::cerr << "WARNING: no elements in CGAL::make_mesh_3 result." << std::endl;
    }
    for (auto it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it) {
        elements.emplace_back(V[it->vertex(0)],
                              V[it->vertex(1)],
                              V[it->vertex(2)],
                              V[it->vertex(3)]);
    }
}
