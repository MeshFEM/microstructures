#include "DynamicKdTree.hh"
#include <MeshFEM/Utilities/apply.hh>
#include <MeshFEM/CollisionGrid.hh>
#include <nanoflann.hpp>
#include <algorithm>
#include <iterator>
#include <type_traits>

// Set from embedded graph
template<class Sym>
void WireMesh<Sym>::
set(const std::vector<MeshIO::IOVertex > &inVertices,
    const std::vector<MeshIO::IOElement> &inElements)
{
    for (const auto &e : inElements)
        if (e.size() != 2) throw std::runtime_error("Expected line element mesh.");

    // Clear old state.
    m_fullVertices.clear(), m_fullEdges.clear();
    m_baseVertices.clear(), m_baseEdges.clear();
    m_baseVertexVarOffsets.clear();
    m_baseVertexPositioners.clear();

    m_fullVertices.reserve(inVertices.size());
    m_fullEdges.reserve(inElements.size());
    for (const auto &e : inElements) {
        assert(e.size() == 2);
        m_fullEdges.push_back({e[0], e[1]});
    }

    // Convert vertices to Point type, scaling into [-1, 1]
    BBox<Point> bb(inVertices);
    auto dim = bb.dimensions();
    if ((std::abs(dim[0]) < 1e-6) || (std::abs(dim[1]) < 1e-6))
        throw std::runtime_error("Degenerate pattern");
    const bool is2D = (std::abs(dim[2]) <= 1e-6);
    if (std::is_same<Symmetry::NonPeriodic<typename PatternSymmetry::Tolerance, 2>, Sym>::value ||
        std::is_same<Symmetry::NonPeriodic<typename PatternSymmetry::Tolerance, 3>, Sym>::value) {

        for (const auto &v : inVertices) {
            Point p(Point::Zero());
            // Hack for 2D case: leave z coords at 0
            for (size_t i = 0; i < (is2D ? 2 : 3); ++i)
                p[i] = v[i];
            m_fullVertices.push_back(p);
        }
    }
    else {
        for (const auto &v : inVertices) {
            // transform graph to [-1, 1]
            Point p(Point::Zero());
            // Hack for 2D case: leave z coords at 0
            for (size_t i = 0; i < (is2D ? 2 : 3); ++i)
                p[i] = (v[i] - bb.minCorner[i]) * (2.0 / dim[i]) - 1.0;
            m_fullVertices.push_back(p);
        }
    }

    // Determine the vertices in the symmetry base unit subgraph
    std::vector<int> baseVertexIndex(m_fullVertices.size(), -1);
    for (size_t i = 0; i < m_fullVertices.size(); ++i) {
        if (PatternSymmetry::inBaseUnit(m_fullVertices[i])) {
            baseVertexIndex[i] = m_baseVertices.size();
            m_baseVertices.push_back(m_fullVertices[i]);
        }
    }

    // Compute edges in the induced subgraph
    for (const auto &e : m_fullEdges) {
        int u = baseVertexIndex.at(e.first), v = baseVertexIndex.at(e.second);
        if ((u >= 0) && (v >= 0))
            m_baseEdges.push_back({u, v});
    }

    // Create positioners for each base vertex
    m_baseVertexPositioners.reserve(m_baseVertices.size());
    for (const auto &p : m_baseVertices) {
        m_baseVertexPositioners.push_back(PatternSymmetry::nodePositioner(p));
    }

    // Determine the "independent" and "dependent" base vertices
    // First, assume there are no dependent vertices
    m_numIndepBaseVertices = numBaseVertices();
    m_numDepBaseVertices = 0;
    m_indepVtxForBaseVtx.resize(numBaseVertices());
    auto isIndep = [&](size_t u) { return m_indepVtxForBaseVtx.at(u) == u; };
    std::iota(m_indepVtxForBaseVtx.begin(), m_indepVtxForBaseVtx.end(), 0);
    // Next, link dependent vertices to their corresponding independent vertex
    for (size_t u = 0; u < numBaseVertices(); ++u) {
        const Point &up = m_baseVertices.at(u);
        const Point &vp = PatternSymmetry::independentVertexPosition(up);
        // Most vertices will be independent (coincide with their independent
        // vertex position), so check this first for efficiency
        if ((vp - up).norm() < PatternSymmetry::tolerance)
            continue;

        // For dependent vertices, we must search for corresponding indep vtx
        m_indepVtxForBaseVtx.at(u) = m_findBaseVertex(vp);
        --m_numIndepBaseVertices;
        ++m_numDepBaseVertices;
    }

    // Enumerate variables
    m_baseVertexVarOffsets.resize(numBaseVertices());
    for (size_t i = 0; i < numBaseVertices(); ++i)
        m_baseVertexVarOffsets[i].blendingPoly.resize(m_blendingPolySize);


    // Create variables for each independent vertex.
    {
        size_t varOffset = 0;
        // Create position vars
        for (size_t i = 0; i < numBaseVertices(); ++i) {
            if (!isIndep(i)) continue;
            m_baseVertexVarOffsets[i].position = varOffset;
            varOffset += m_baseVertexPositioners[i].numDoFs();
        }
        m_numPositionParams = varOffset;

        // Create thickness vars
        for (size_t i = 0; i < numBaseVertices(); ++i) {
            if (!isIndep(i)) continue;
            m_baseVertexVarOffsets[i].thickness = varOffset++;
        }

        // Create blending vars
        for (size_t i = 0; i < numBaseVertices(); ++i) {
            if (!isIndep(i)) continue;
            m_baseVertexVarOffsets[i].blending = varOffset++;
        }

        // Create blending poly vars
        for (size_t d = 0; d < m_blendingPolySize; d++) {
            for (size_t i = 0; i < numBaseVertices(); ++i) {
                if (!isIndep(i)) continue;
                m_baseVertexVarOffsets[i].blendingPoly[d] = varOffset++;
            }
        }
    }

    // Link dependent vertices to the independent vertices' variables.
    for (size_t u = 0; u < numBaseVertices(); ++u)
        m_baseVertexVarOffsets[u] = m_baseVertexVarOffsets[m_indepVtxForBaseVtx[u]];

    ////////////////////////////////////////////////////////////////////////////
    // Construct inflation graph
    ////////////////////////////////////////////////////////////////////////////
    // Determine vertices/edges in the inflation graph
    std::vector<TransformedVertex> rVertices;
    std::vector<TransformedEdge>   rEdges;
    replicatedGraph(symmetryGroup(), rVertices, rEdges);

    // We want to include all vertices either inside the symmetry base unit or
    // m_inflationNeighborhoodEdgeDist edges away. First, we compute the edge
    // distance to the symmetry base unit.
    std::vector<size_t> vtxEdgeDist;
    std::queue<size_t> bfsQ;
    for (size_t i = 0; i < rVertices.size(); ++i) {
        size_t dist = std::numeric_limits<size_t>::max();
        if (PatternSymmetry::inBaseUnit(rVertices[i].pt)) {
            dist = 0;
            bfsQ.push(i);
        }
        vtxEdgeDist.push_back(dist);
    }

    std::vector<std::vector<size_t>> adj(rVertices.size());
    for (const auto &re : rEdges) {
        adj.at(re.e[0]).push_back(re.e[1]);
        adj.at(re.e[1]).push_back(re.e[0]);
    }

    while (!bfsQ.empty()) {
        size_t u = bfsQ.front();
        bfsQ.pop();
        size_t uDist = vtxEdgeDist[u];
        assert(uDist < std::numeric_limits<size_t>::max());
        for (size_t v : adj[u]) {
            if (uDist + 1 < vtxEdgeDist[v]) {
                vtxEdgeDist[v] = uDist + 1;
                bfsQ.push(v);
            }
        }
    }

    auto keepVtx = MeshFEM::apply(vtxEdgeDist, [=](size_t dist) { return dist <= m_inflationNeighborhoodEdgeDist; });

    // Reindex kept vertices
    m_inflVtx.clear();
    std::vector<size_t> vertexRenumber(rVertices.size(), std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < rVertices.size(); ++i) {
        if (!keepVtx[i]) continue;
        vertexRenumber[i] = m_inflVtx.size();
        m_inflVtx.push_back(rVertices[i]);
    }

    // Copy over edges m_inflationNeighborhoodEdgeDist away from cell
    // (i.e. the induced graph of kept vertices, but excluding the edges
    // between the frontier vertices)
    m_inflEdge.clear(), m_inflEdgeOrigin.clear();
    for (const auto &re : rEdges) {
        bool keepEdge = vtxEdgeDist[re.e[0]] < m_inflationNeighborhoodEdgeDist
                     || vtxEdgeDist[re.e[1]] < m_inflationNeighborhoodEdgeDist;
        // Always keep edges inside the inflation graph.
        keepEdge |= (vtxEdgeDist[re.e[0]] == 0) && (vtxEdgeDist[re.e[1]] == 0);
        if (!keepEdge) continue;
        size_t u = vertexRenumber.at(re.e[0]),
               v = vertexRenumber.at(re.e[1]);
        assert((u < m_inflVtx.size()) && (v < m_inflVtx.size()));
        m_inflEdge.push_back({u, v});
        m_inflEdgeOrigin.push_back(re.origEdge);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Construct full period cell graph
    ////////////////////////////////////////////////////////////////////////////
    {
        std::vector<Isometry> pcellIsometries;
        for (const auto &iso : symmetryGroup()) {
            if (iso.hasTranslation()) continue;
            pcellIsometries.push_back(iso);
        }
        std::vector<TransformedEdge> pEdges;
        replicatedGraph(pcellIsometries, m_periodCellVtx, pEdges);
        for (const auto &pe : pEdges)
            m_periodCellEdge.push_back({pe.e[0], pe.e[1]});
    }

    ////////////////////////////////////////////////////////////////////////////
    // Construct printability graph
    ////////////////////////////////////////////////////////////////////////////
    // Reflect base graph into the "printability column." This is done by
    // applying all permutation isometries, but no translation and only the
    // z-axis reflection
    std::vector<Isometry> printabiltyIsometries;
    for (const auto &iso : symmetryGroup()) {
        if (iso.hasTranslation() ||
            iso.hasReflection(Symmetry::Axis::X) ||
            iso.hasReflection(Symmetry::Axis::Y)) {
            continue;
        }
        printabiltyIsometries.push_back(iso);
    }

    std::vector<TransformedEdge> pEdges;
    replicatedGraph(printabiltyIsometries,
                    m_printGraphVtx, pEdges);
    m_printGraphEdge.clear(), m_printGraphEdge.reserve(pEdges.size());
    for (const auto &pe : pEdges)
        m_printGraphEdge.push_back({pe.e[0], pe.e[1]});
}

// Construct (stitched) replicated graph along with the maps from parameters
// to vertex positions/thicknesses/blending parameters
// Note: these maps operate on "homogeneous parameters" (i.e. the vector of
// parameters with 1 appended).
template<class Sym>
void WireMesh<Sym>::
replicatedGraph(const std::vector<Isometry> &isometries,
                std::vector<TransformedVertex> &outVertices,
                std::vector<TransformedEdge  > &outEdges) const
{
    std::vector<TransformedVertex> rVertices;
    std::set<TransformedEdge> rEdges;

    const double tol = m_tolerance();

    // Choose a cell size on the order of epsilon, but prevent cell sizes so
    // small as to cause index overflows for objects of size up to 100x100
    // centered at the origin: max int ~10^9 ==> cellSize > 10^-7
    // CollisionGrid<double, Point> cgrid(std::max(tol, 1e-7));
    ::micro::DynamicKdTree3d cTree;

    for (size_t ei = 0; ei < m_baseEdges.size(); ++ei) {
        const auto &e = m_baseEdges[ei];
        const size_t u = e.first, v = e.second;

        auto pu = m_baseVertices.at(u),
             pv = m_baseVertices.at(v);

        Eigen::Matrix3Xd uPosMap = m_baseVertexPositioners.at(u).getPositionMap(numParams(), m_baseVertexVarOffsets.at(u).position);
        Eigen::Matrix3Xd vPosMap = m_baseVertexPositioners.at(v).getPositionMap(numParams(), m_baseVertexVarOffsets.at(v).position);

        for (const auto &isometry : isometries) {
            auto mappedPu = isometry.apply(pu);
            auto mappedPv = isometry.apply(pv);

            // Create a new replicated vertex at "p" unless p coincides
            // with an existing one. Return resulting unique vertex's index
            auto createUniqueVtx = [&](const size_t origVertex, const Point &p, const Eigen::Matrix3Xd &posMap) {
                auto result = cTree.getClosestPoint(p, tol);
                if (result.first < 0) {
                    // Create new replicated vertex at "p" if none exists
                    size_t newVtxIdx = rVertices.size();
                    rVertices.emplace_back(p, origVertex, isometry, posMap);
                    cTree.addPoint(p, newVtxIdx);
                    return newVtxIdx;
                }

                size_t match = result.first;
                // Verify position maps agree for coinciding reflected vertices
                // that originate from a different independent vtx
                if ((rVertices[match].posMap - posMap).squaredNorm() >= 1e-8)
                {
                    std::cerr << m_baseVertices.at(rVertices[match].origVertex).transpose() << std::endl;
                    std::cerr << m_baseVertices.at(origVertex).transpose() << std::endl;
                    std::cerr << rVertices[match].iso << std::endl;
                    std::cerr << rVertices[match].posMap << std::endl;
                    std::cerr << isometry << std::endl;
                    std::cerr << posMap << std::endl;
                    assert(false && "Conflicting maps");
                }
                // Verify repeated vertices originate from same independent vtx
                assert(m_indepVtxForBaseVtx.at(rVertices[match].origVertex) ==
                        m_indepVtxForBaseVtx.at(origVertex));
                return match;
            };

            size_t mappedUIdx = createUniqueVtx(u, mappedPu, isometry.xformMap(uPosMap)),
                   mappedVIdx = createUniqueVtx(v, mappedPv, isometry.xformMap(vPosMap));

            TransformedEdge xfEdge(mappedUIdx, mappedVIdx, ei);
            auto ret = rEdges.insert(xfEdge);
            if (ret.second == false) {
                // reflected already exists; verify it's from the same base edge
                // (unless this is a triply periodic pattern where base edges
                // on the period cell boundary are redundant)
                assert(((ret.first->origEdge == ei) ||
                        std::is_same<Symmetry::TriplyPeriodic<typename PatternSymmetry::Tolerance>, PatternSymmetry>::value ||
                        std::is_same<Symmetry::DoublyPeriodic<typename PatternSymmetry::Tolerance>, PatternSymmetry>::value));
            }
        }
    }

    outVertices = rVertices;
    outEdges.clear();
    outEdges.insert(outEdges.end(), rEdges.begin(), rEdges.end());
}

// Set from embedded graph stored in obj/msh format.
inline void WireMeshBase::
load(const std::string &wirePath) {
    std::vector<MeshIO::IOVertex> inVertices;
    std::vector<MeshIO::IOElement> inElements;
    MeshIO::load(wirePath, inVertices, inElements);
    set(inVertices, inElements);
}

template<class Point>
void _OutputGraph(const std::string &path, const std::vector<Point> &points,
                  const std::vector<std::pair<size_t, size_t>> &edges) {
    std::vector<MeshIO::IOVertex>  outVertices;
    std::vector<MeshIO::IOElement> outElements;
    outVertices.reserve(points.size()), outElements.reserve(edges.size());
    for (const Point &p : points) outVertices.emplace_back(p);
    for (const auto  &e :  edges) outElements.emplace_back(e.first, e.second);
    MeshIO::save(path, outVertices, outElements);
}

// Save the full embedded graph in obj/msh format.
inline void WireMeshBase::
save(const std::string &path) const {
    _OutputGraph(path, m_fullVertices, m_fullEdges);
}

// Save the embedded base unit graph in obj/msh format.
inline void WireMeshBase::
saveBaseUnit(const std::string &path) const {
    _OutputGraph(path, m_baseVertices, m_baseEdges);
}

// For symmetry debugging:
// save the tiled base unit graph in obj/msh format
inline void WireMeshBase::
saveReplicatedBaseUnit(const std::string &path) const {
    std::vector<MeshIO::IOVertex> outVertices;
    std::vector<MeshIO::IOElement> outEdges;
    for (const auto &iso : symmetryGroup()) {
        size_t offset = outVertices.size();
        for (const Point &p : m_baseVertices) { outVertices.emplace_back(iso.apply(p)); }
        for (const Edge  &e :    m_baseEdges) { outEdges.emplace_back(e.first + offset, e.second + offset); }
    }
    MeshIO::save(path, outVertices, outEdges);
}

inline void WireMeshBase::
savePeriodCellGraph(const std::string &path) const {
    std::vector<Point> outVertices;
    std::vector<Edge>  outEdges;
    periodCellGraph(outVertices, outEdges);
    _OutputGraph(path, outVertices, outEdges);
}

// For symmetry debugging:
// save the inflation graph in obj/msh format
inline void WireMeshBase::
saveInflationGraph(const std::string &path, std::vector<double> params) const {
    if (params.size() == 0)
        params = defaultParameters();
    std::vector<Edge> igraphEdges;
    std::vector<Point> igraphVertices;
    std::vector<double> thicknesses, blendingParams;
    std::vector<std::vector<double>> blendingPolyCoeffs;
    inflationGraph(params, igraphVertices, igraphEdges, thicknesses, blendingParams, blendingPolyCoeffs);

    _OutputGraph(path, igraphVertices, igraphEdges);
}

// Find the base vertex within symmetry tolerance of p
// (or throw an exception if none exists).
// The number of base vertices to check is generally small, so fancy
// data structures shouldn't be needed.
template<class Sym>
size_t WireMesh<Sym>::
m_findBaseVertex(const Point &p) const {
    for (size_t i = 0; i < m_baseVertices.size(); ++i) {
        if ((p - m_baseVertices[i]).squaredNorm() < (PatternSymmetry::tolerance * PatternSymmetry::tolerance))
            return i;
    }
    std::cout << "Failed to find " << p.transpose() << std::endl;
    for (const Point &v : m_baseVertices)
        std::cout << "    candidate " << v.transpose() << std::endl;
    throw std::runtime_error("Couldn't find corresponding base vertex.");
}
