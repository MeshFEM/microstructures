#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#include "WireQuadMesh.hh"
#include "MeshingOptions.hh"
#include <igl/doublearea.h>
#include <igl/write_triangle_mesh.h>
#include "quadfoam/instantiate.h"
#include "quadfoam/jacobians.h"
#include <cmath>
#include <map>
////////////////////////////////////////////////////////////////////////////////

namespace {

std::string lowercase(std::string data) {
    std::transform(data.begin(), data.end(), data.begin(), ::tolower);
    return data;
}

#define TRY_SYMMETRY(s, x, v, e)                                    \
    if (lowercase(x) == lowercase(#s))                              \
    {                                                               \
        return std::make_shared<WireMesh<Symmetry::s<>>>((v), (e)); \
    }

#define TRY_KEY_VAL(s, a, x, v, e)                                  \
    if (lowercase(x) == lowercase(#a))                              \
    {                                                               \
        return std::make_shared<WireMesh<Symmetry::s<>>>((v), (e)); \
    }

std::shared_ptr<WireMeshBase> load_wire_mesh(
    const std::string &sym,
    const std::vector<MeshIO::IOVertex > &vtx,
    const std::vector<MeshIO::IOElement> &elm)
{
    TRY_SYMMETRY(Square, sym, vtx, elm);
    TRY_SYMMETRY(Cubic, sym, vtx, elm);
    TRY_SYMMETRY(Orthotropic, sym, vtx, elm);
    TRY_SYMMETRY(Diagonal, sym, vtx, elm);
    TRY_KEY_VAL(DoublyPeriodic, Doubly_Periodic, sym, vtx, elm);
    TRY_KEY_VAL(TriplyPeriodic, Triply_Periodic, sym, vtx, elm);
    TRY_KEY_VAL(NonPeriodic, Non_Periodic, sym, vtx, elm);
    throw std::runtime_error("Could not find a matching symmetry type to load for: " + sym);
}

// -----------------------------------------------------------------------------

template<typename T>
void sort_unique(std::vector<T> &x) {
    std::sort(x.begin(), x.end());
    auto it = std::unique(x.begin(), x.end());
    x.resize(std::distance(x.begin(), it));
}

Eigen::MatrixXd to_eigen_matrix(const std::vector<WireQuadMesh::Point> &verts) {
    Eigen::MatrixXd V(verts.size(), 3);
    for (int i = 0; i < (int) verts.size(); ++i) {
        V.row(i) = verts[i].transpose();
    }
    return V;
}

std::vector<WireQuadMesh::Point> to_std_vector(const Eigen::MatrixXd &V) {
    std::vector<WireQuadMesh::Point> verts(V.rows());
    for (int i = 0; i < (int) verts.size(); ++i) {
        verts[i] = V.row(i).transpose();
    }
    return verts;
}

// -----------------------------------------------------------------------------

void build_stitched_graph(
    const Eigen::MatrixXd & reduced_positions,
    const Eigen::VectorXi & full_to_reduced,
    const std::vector<WireQuadMesh::Edge> & all_edges,
    std::vector<std::vector<size_t>> & stitchedVertices,
    std::vector<WireQuadMesh::Edge> & stitchedEdges)
{
    stitchedVertices.resize(reduced_positions.rows(), {});
    for (size_t i = 0; i < (size_t) full_to_reduced.size(); ++i) {
        int j = full_to_reduced[i];
        assert(j < reduced_positions.size());
        stitchedVertices[j].push_back(i);
    }
    for (auto &e : all_edges) {
        stitchedEdges.emplace_back(full_to_reduced[e.first], full_to_reduced[e.second]);
    }
    sort_unique(stitchedEdges);
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

WireQuadMesh::WireQuadMesh(
    const std::vector<MeshIO::IOVertex> &V,
    const std::vector<MeshIO::IOElement> &F,
    const nlohmann::json &params)
{
    // Read input parameters
    size_t numQuads = F.size();
    m_allTopologies.clear(); m_allTopologies.resize(numQuads);
    m_allParameters.clear(); m_allParameters.resize(numQuads);
    m_allJacobians.clear(); m_allJacobians.resize(numQuads);
    bool need_compute_jacobians = false;
    typedef std::pair<std::vector<MeshIO::IOVertex>, std::vector<MeshIO::IOElement>> Mesh;
    std::map<std::string, Mesh> cache;
    for (size_t i = 0; i < numQuads; ++i) {
        auto entry = params[i];
        std::string wirePath = entry["pattern"];
        if (!cache.count(wirePath)) {
            std::vector<MeshIO::IOVertex> inVertices;
            std::vector<MeshIO::IOElement> inElements;
            MeshIO::load(wirePath, inVertices, inElements);
            cache[wirePath] = std::make_pair(inVertices, inElements);
        }
        const auto &mesh = cache[wirePath];
        m_allParameters[i] = entry["params"].get<std::vector<double>>();
        m_allTopologies[i] = load_wire_mesh(entry["symmetry"], mesh.first, mesh.second);
        assert(m_allTopologies[i]->thicknessType() == m_thicknessType);
        assert(m_allTopologies[i]->numParams() == m_allParameters[i].size());
        if (entry.count("jacobian")) {
            m_allJacobians[i] = MeshingOptions::read_jacobian(entry["jacobian"]);
        } else {
            need_compute_jacobians = true;
        }
    }

    // Copy background quad mesh
    m_V.resize(V.size(), 3);
    m_F.resize(F.size(), 4);
    for (int v = 0; v < (int) V.size(); ++v) {
        for (int c = 0; c < 3; ++c) {
            m_V(v, c) = V[v][c];
        }
    }
    for (int f = 0; f < (int) F.size(); ++f) {
        assert(F[f].size() == 4);
        for (int lv = 0; lv < 4; ++lv) {
            m_F(f, lv) = F[f][lv];
        }
    }

    if (need_compute_jacobians) {
        std::cerr << "Missing Jacobian information, will compute values based on input mesh." << std::endl;
        compute_jacobians();
    }

    // Compute vertex id offset in concatenated graph
    std::vector<Eigen::MatrixXd> allVertices;
    std::vector<Edge> allEdges;
    for (size_t f = 0, vertex_offset = 0; f < numQuads; ++f) {
        std::vector<Point> verts;
        std::vector<Edge>  edges;
        m_allTopologies[f]->periodCellGraph(verts, edges);
        allVertices.push_back(to_eigen_matrix(verts));
        for (const auto &e : edges) { allEdges.emplace_back(e.first + vertex_offset, e.second + vertex_offset); }
        vertex_offset += verts.size();
        // _OutputGraph("test_inflation_graph.obj", verts, edges);
    };

    // Stitch input graph and remove duplicate vertices
    Eigen::MatrixXd reducedPositions;
    Eigen::MatrixXi reducedFaces;
    Eigen::VectorXi fullToReduced;
    bool success = micro_quadfoam::instantiate_pattern(m_V, m_F, allVertices, {},
        reducedPositions, reducedFaces, true, PatternSymmetry::tolerance, &fullToReduced);
    if (!success) {
        throw std::runtime_error("Could not stitch the wiremesh together!");
    }

    // Build inverse map (list of full ids mapped to a single id)
    std::vector<std::vector<size_t>> stitchedVertices;
    std::vector<WireQuadMesh::Edge> stitchedEdges;
    build_stitched_graph(reducedPositions, fullToReduced, allEdges,
        stitchedVertices, stitchedEdges);

    // [Debug]
    // std::vector<Point3d> allPts;
    // for (int i = 0; i < reducedPositions.rows(); ++i) {
    //     Eigen::RowVector3d x = reducedPositions.row(i);
    //     allPts.push_back(Point3d(x[0], x[1], x[2]));
    // }
    // _OutputGraph("test_inflation_graph.msh", allPts, stitchedEdges);

    // Set active quad to -1 (initialize bilinear map and bounding box)
    setActiveQuad(-1);
}

////////////////////////////////////////////////////////////////////////////////

// Fill m_allJacobians based on the background quad mesh (m_V, m_F)
void WireQuadMesh::compute_jacobians() {
    using json = nlohmann::json;
    Eigen::MatrixXd Q, P;
    micro_quadfoam::jacobians(m_V, m_F, Q, P);
    assert(m_F.rows() == (int) m_allJacobians.size());
    for (int i = 0; i < m_F.rows(); ++i) {
        json obj = json::array({Q(i, 0), Q(i, 1), Q(i, 2), Q(i, 3)});
        m_allJacobians[i] = MeshingOptions::read_jacobian(obj);
    }
}

BilinearMap WireQuadMesh::getBilinearMap(int i) const {
    Point2d pts[4];
    for (int lv = 0; lv < 4; ++lv) {
        int v = m_F(i, lv);
        pts[lv][0] = m_V(v, 0);
        pts[lv][1] = m_V(v, 1);
    }
    return BilinearMap(pts);
}

// Compute a global scaling factor between the reference domain [-1, 1]^2 to the
// actual quad on the background mesh.
double WireQuadMesh::getScalingFactor(int i) const {
    Point2d pts[4];
    double len = 0;
    for (int lv = 0; lv < 4; ++lv) {
        int v0 = m_F(i, lv);
        int v1 = m_F(i, (lv + 1)%4);
        len += (m_V.row(v0) - m_V.row(v1)).stableNorm();
    }
    return len / 8.0;
}

////////////////////////////////////////////////////////////////////////////////

// Set currently active quad
void WireQuadMesh::setActiveQuad(int idx) {
    if (idx < 0) {
        // Meshing domain = the whole mesh
        m_bilinearMap = BilinearMap(); // identity
        Point3d minV = m_V.colwise().minCoeff().transpose();
        Point3d maxV = m_V.colwise().maxCoeff().transpose();
        m_bbox = BBox<Point3d>(minV, maxV);
    } else {
        // Meshing domain restricted to the reference square [-1,1]²
        m_activeQuad = idx;
        m_bilinearMap = getBilinearMap(idx);
        m_bbox = BBox<Point3d>(Point3d(-1, -1, 0), Point3d(1, 1, 0));
    }
}

// Return wiremesh associated to the currently active quad
const WireMeshBase &WireQuadMesh::activeWireMesh() const {
    if (m_activeQuad < 0) {
        throw std::runtime_error("No active quad has been selected!");
    } else {
        return *m_allTopologies[m_activeQuad];
    }
}

// Inflation parameters the whole graph (simple concatenation)
std::vector<double> WireQuadMesh::params() const {
    // Otherwise simply concatenate all parameters together
    std::vector<double> params;
    for (const auto &p : m_allParameters) {
        params.insert(params.end(), p.begin(), p.end());
    }
    return params;
}

Eigen::VectorXd WireQuadMesh::areas() const {
    Eigen::VectorXd A;
    igl::doublearea(m_V, m_F, A);
    return A;
}

// -----------------------------------------------------------------------------

// Build the inflation graph for the active quad mesh, stitching together adjacent
// nodes (averaging stitched points' locations, thicknesses, and blending params).
//
// @param[in]  allParams               { Inflation parameters for the whole mesh }
// @param[out] stitchedPoints          { Graph vertices positions }
// @param[out] stitchedEdges           { Graph edges indices }
// @param[out] stitchedThicknesses     { Graph edge thicknesses }
// @param[out] stitchedBlendingParams  { Graph vertices blending parameters }
//
void WireQuadMesh::inflationGraph(const std::vector<double> &allParams,
    std::vector<WireQuadMesh::Point> &stitchedPoints,
    std::vector<WireQuadMesh::Edge> &stitchedEdges,
    std::vector<double> &stitchedThicknesses,
    std::vector<double> &stitchedBlendingParams,
    std::vector<std::vector<double>> &stitchedBlendingPolyParams
) const
{
    std::vector<Eigen::MatrixXd> allVertices;
    std::vector<Edge> allEdges;
    std::vector<double> allThicknesses;
    std::vector<double> allBlendingParams;

    // Build period cell graph for each wire mesh separately
    int paramOffset = 0;
    for (int i = 0, vertex_offset = 0; i < m_F.rows(); ++i) {
        const auto wmesh = m_allTopologies[i];
        const size_t np = wmesh->numParams();
        std::vector<Real> wmparams(&allParams[paramOffset], &allParams[paramOffset] + np);
        paramOffset += np;
        std::vector<Point> points;
        std::vector<Edge> edges;
        std::vector<double> thicknesses;
        std::vector<double> blendingParams;
        wmesh->periodCellGraph(wmparams, points, edges, thicknesses, blendingParams);

        // Graph parameters (radii & blending) for each cell are expressed in
        // the parallelogram space (i.e. after mapping the reference square [-1,1]²
        // to the target parallelogram).
        //
        // Whether we are meshing a single cell or the whole mesh, the evaluation
        // of the SDF is performed in world space (the i.e. after mapping to the
        // actual quadrilateral on the mesh).
        //
        // The scaling we thus need to apply is the Jacobian of the transformation
        // that goes from the pattern's target parallelogram to the actual physical
        // quadrilateral of the cell (which may not be a parallelogram).

        double jac_det = m_allJacobians[i].determinant();
        if (jac_det <= 0) {
            std::cerr << "Warning: Jacobian of the mapped parallelogram has negative volume for quad " << i << std::endl;
        }
        assert(jac_det > 0);
        double pre_scaling = getScalingFactor(i) / std::sqrt(std::abs(jac_det));
        auto map = getBilinearMap(i);

        allVertices.push_back(to_eigen_matrix(points));
        for (const auto &e : edges) { allEdges.emplace_back(e.first + vertex_offset, e.second + vertex_offset); }

        int numVertices = (int) points.size();
        assert(numVertices == (int) thicknesses.size());
        assert(numVertices == (int) blendingParams.size());
        for (int lv = 0; lv < numVertices; ++lv) {
            double jdet = map.jacobian(points[lv][0], points[lv][1]).determinant();
            if (jdet <= 0 || std::isnan(jdet)) {
                std::cerr << "Warning: Mapping local vertex [" << lv << "] for quad [" << i
                    << "] onto the reference square [-1,1]^2 has negative volume: ("
                    << points[lv][0] << "," << points[lv][1] << ")" << std::endl;
            }
            assert(jdet > 0);
            double scaling = pre_scaling * std::sqrt(std::abs(jdet));
            assert(!std::isnan(scaling));
            allThicknesses.push_back(scaling * thicknesses[lv]);
            allBlendingParams.push_back(scaling * blendingParams[lv]);
        }
        vertex_offset += points.size();
    };

    // Stitch input graph and remove duplicate vertices
    Eigen::MatrixXd reducedPositions;
    Eigen::MatrixXi reducedFaces;
    Eigen::VectorXi fullToReduced;
    bool success = micro_quadfoam::instantiate_pattern(m_V, m_F, allVertices, {},
        reducedPositions, reducedFaces, true, -1, &fullToReduced);
    if (!success) {
        throw std::runtime_error("Could not stitch the wiremesh together!");
    }

    // Build inverse map (list of full ids mapped to a single id)
    std::vector<std::vector<size_t>> stitchedVertices;
    build_stitched_graph(reducedPositions, fullToReduced, allEdges,
        stitchedVertices, stitchedEdges);

    // Extracted stitched graph from the replicated graph by averaging
    // merged vertices' values.
    stitchedPoints = to_std_vector(reducedPositions);
    stitchedThicknesses.clear();
    stitchedBlendingParams.clear();
    const size_t nsv = stitchedVertices.size();
    stitchedThicknesses.reserve(nsv);
    stitchedBlendingParams.reserve(nsv);
    for (const auto &sv : stitchedVertices) {
        // Take mean of location, thickness, and blending of all vertices
        // merged into this stitched output vertex.
        Point pt(Point::Zero());
        double t = 0, b = 0.0;
        for (size_t v : sv) {
            // pt += allVerts.at(v);
            t += allThicknesses.at(v);
            b += allBlendingParams.at(v);
        }

        assert(!std::isnan(t));
        assert(sv.size() > 0);
        // stitchedPoints.push_back(pt / sv.size());
        stitchedThicknesses.push_back(t / sv.size());
        stitchedBlendingParams.push_back(b / sv.size());
    }

    stitchedBlendingPolyParams.clear();
    stitchedBlendingPolyParams.resize(stitchedVertices.size());

    // Now, for each vertices, scale back the (radius, blending) parameters based
    // on the jacobian of the bilinear map (that maps the ref square [-1,1]² to
    // the active quad).

    // for (int i = 0; i < (int) stitchedPoints.size(); ++i) {
    //     Point3d p = stitchedPoints[i];
    //     auto jac = m_bilinearMap.jacobian(p[0], p[1]);
    //     double scaling = 1.0 / std::sqrt(jac.determinant());
    //     std::cout << scaling << std::endl;

    //     // stitchedThicknesses[i] *= scaling;
    //     // stitchedBlendingParams[i] *= scaling;
    // }

    // _OutputGraph("test_inflation_graph.obj", stitchedPoints, stitchedEdges);
}

#endif
