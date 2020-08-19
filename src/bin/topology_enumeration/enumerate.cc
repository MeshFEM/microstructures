// Topology enumeration:
//
// Enumerate using barycentric coordinates. 2D enumeration is equivalent to a 3D
// enumeration that considers only patterns lying perfectly in the midplane (and
// ignores isometries mapping out of the midplane)
//
// Brute-force point merging.
#include "intersection_check.hh"
#include <isosurface_inflator/Symmetry.hh>
#include <MeshFEM/MeshIO.hh>
#include <vector>
#include <array>
#include <limits>
#include <utility>
#include <set>
#include <queue>

#include <sstream>
#include <iomanip>

// 3D enumeration by default.
#ifndef DIM
#define DIM 3
#endif

#define MAX_VALENCE 7

using namespace std;
using BaryPoint = array<Real, 4>;
using Sym = Symmetry::Cubic<>;

constexpr size_t NONE = numeric_limits<size_t>::max();

vector<Isometry> symmetryGroup = Sym::symmetryGroup();
Real epsilon = double(Sym::Tolerance::num) / double(Sym::Tolerance::den);
Real epsilonSq = epsilon * epsilon;

vector<BaryPoint> nodesBary;
vector<Point3d> nodes;
vector<pair<size_t, size_t>> edges;

size_t numChecked = 0;
size_t nonSelfIntersecting = 0;
size_t connectedCount = 0;
size_t noCoincidingEdges = 0;
size_t noDanglingEdges = 0;
size_t valenceTresholded = 0;
size_t validCount = 0;

// taken: which edges are chosen
void processGraph(const vector<size_t> &taken) {
    // Extract subgraph
    vector<pair<size_t, size_t>> subEdges;
    vector<Point3d> subNodes;
    vector<size_t> subNodeForNode(nodes.size(), NONE);

    bool verbose = false;
    // if ((taken[0] == 31) && (taken[1] == 58) && (taken[2] == 82)) {
    //     verbose = true;
    //     cout << "verbose mode" << endl;
    // }

    for (size_t ei : taken) {
        size_t u, v;
        tie(u, v) = edges.at(ei);
#if DIM == 2
        // Only consider subgraphs in the 2D base triangle.
        if ((fabs(nodes.at(u)[2]) > 1e-8) || (fabs(nodes.at(v)[2]) > 1e-8)) {
            if (verbose) std::cout << "Not in base triangle" << endl;
            return; // BAIL early (before incrementing numChecked)
        }
#endif
        if (subNodeForNode.at(u) == NONE) {
            subNodeForNode.at(u) = subNodes.size();
            subNodes.push_back(nodes.at(u));
        }

        if (subNodeForNode.at(v) == NONE) {
            subNodeForNode.at(v) = subNodes.size();
            subNodes.push_back(nodes.at(v));
        }

        subEdges.push_back({subNodeForNode.at(u),
                            subNodeForNode.at(v)});
    }

    ++numChecked;
    if (numChecked % 1000 == 0) cout << numChecked << endl;

    // Check for self-intersections. We define these as edges that intersect
    // while not being neighbors.
    if (hasSelfIntersection(subNodes, subEdges)) {
        if (verbose) std::cout << "Has self-intersection" << endl;
        return; // BAIL
    }
    ++nonSelfIntersecting;

    // Replicate subgraph. Note: this will also replicate neighboring
    // subgraphs. This is useful for checking dangling vertices/valences since
    // we don't need special cases for the border vertices.
    // However, we must make sure only to check the graph inside the base
    // cell...
    set<pair<size_t, size_t>> replicatedEdges;
    vector<Point3d>           replicatedNodes;
    replicatedNodes.reserve(symmetryGroup.size() * subNodes.size());
    // replicatedEdges.reserve(symmetryGroup.size() * subEdges.size());

    vector<size_t> nodeRemapper;
    // Brute-force merging of replicated graph
    for (const auto &iso : symmetryGroup) {
        nodeRemapper.assign(subNodes.size(), NONE);
        for (size_t i = 0; i < subNodes.size(); ++i) {
            Point3d pt = iso.apply(subNodes[i]);
#if DIM == 2
            // Ignore vertices created outside the midplane
            if (fabs(pt[2]) > 1e-8) continue;
#endif

            size_t vtxIdx = NONE;
            for (size_t j = 0; j < replicatedNodes.size(); ++j) {
                if ((pt - replicatedNodes[j]).squaredNorm() < epsilonSq) {
                    vtxIdx = j;
                    break;
                }
            }
            if (vtxIdx == NONE) {
                vtxIdx = replicatedNodes.size();
                replicatedNodes.push_back(pt);
            }
            nodeRemapper[i] = vtxIdx;
        }

        // using a set ignores replicated edges
        for (const auto &e : subEdges) {
            size_t u = nodeRemapper.at(e.first),
                   v = nodeRemapper.at(e.second);
            if ((u == NONE) || (v == NONE)) {
                assert(DIM == 2); // only 2D should drop vertices
                continue;
            }
            replicatedEdges.insert({u, v});
        }
    }

    // Build traversible graph representation
    vector<vector<size_t>> adj(replicatedNodes.size());
    for (const auto &e : replicatedEdges) {
        adj.at(e.first).push_back(e.second);
        adj.at(e.second).push_back(e.first);
    }

    // Check if connected
    {
        vector<bool> seen(replicatedNodes.size(), false);
        queue<size_t> bfs;
        bfs.push(0);
        seen[0] = true;
        while (!bfs.empty()) {
            size_t u = bfs.front();
            bfs.pop();
            for (size_t v : adj[u]) {
                if (!seen[v]) {
                    seen[v] = true;
                    bfs.push(v);
                }
            }
        }

        bool disconnected = false;
        for (bool b : seen) disconnected |= !b;
        if (disconnected) {
            if (verbose) cout << "disconnected" << endl;
            return; // BAIL
        }
        ++connectedCount;
    }

    // Check for coinciding edges
    // (This really should be done first, but Nico/Luigi did it after
    // connectivity check, so I need to as well to verify number of patterns.)
    {
        // Check one tet edge at a time
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = i + 1; j < 4; ++j) {
                bool hasSubEdge = false;
                bool hasFullEdge = false;
                // See if any of the pattern edges overlap on this tet edge
                for (size_t ei : taken) {
                    const auto &uBary = nodesBary[edges[ei].first];
                    const auto &vBary = nodesBary[edges[ei].second];
                    // Edges were already generated so that "u" < "v"
                    hasFullEdge |= ((uBary[i] == 1.0) && (vBary[j] == 1.0));

                    // Tet edge nodes are placed after tet vertex nodes, so
                    // again we know the ordering: vertex node is first
                    hasSubEdge |= ((uBary[i] == 1.0) || (uBary[j] == 1.0)) && // One of the edge endpoints
                                  ((vBary[i] == 0.5) && (vBary[j] == 0.5));   // the edge midpoint

                    bool hasOverlap = hasSubEdge && hasFullEdge;
                    if (hasOverlap) {
                        if (verbose) cout << "has overlap" << endl;
                        return; // BAIL
                    }
                }
            }
        }

        // We have coinciding edges exactly when a full tet edge is taken and
        // one of the sub tet edges is taken.
        ++noCoincidingEdges;
    }

    // Check for dangling edges (valence 1 vertices in the period cell)
    // and vertices exceeding max valence.
    {
        bool dangling = false, exceedsValence = false;
        for (size_t i = 0; i < replicatedNodes.size(); ++i) {
            const auto &pt = replicatedNodes[i];
            if (!Symmetry::TriplyPeriodic<>::inBaseUnit(pt)) continue;
            const size_t valence = adj[i].size();
            assert(valence > 0);

            if (valence == 1)  { dangling = true; }
            if (valence >  MAX_VALENCE) { exceedsValence = true; }
        }

        if (dangling) {
            if (verbose) cout << "has dangling" << endl;
            return; // BAIL
        }
        ++noDanglingEdges;
        if (exceedsValence) {
            if (verbose) cout << "Exceeds max valence" << endl;
            return; // BAIL
        }
        ++valenceTresholded;
    }

    ++validCount;

    // Write the period cell (only)
    {
        std::vector<MeshIO::IOVertex> outVertices;
        std::vector<MeshIO::IOElement> outElements;

        std::vector<size_t> outVertexForRepNode(replicatedNodes.size(), NONE);
        // Get nodes inside period cell
        for (size_t i = 0; i < replicatedNodes.size(); ++i) {
            const auto &pt = replicatedNodes[i];
            if (!Symmetry::TriplyPeriodic<>::inBaseUnit(pt)) continue;

            outVertexForRepNode[i] = outVertices.size();
            outVertices.emplace_back(pt);
        }

        // Construct induced graph
        for (const auto &e : replicatedEdges) {
            size_t u = outVertexForRepNode.at(e.first),
                   v = outVertexForRepNode.at(e.second);
            if ((u == NONE) || (v == NONE)) continue;
            outElements.emplace_back(u, v);
        }

        stringstream ss;
        ss << setfill('0') << setw(4) << validCount;
        string outName = ss.str() + ".obj";
        cout << outName << endl;
        MeshIO::save(outName, outVertices, outElements);

        // Output subgraph too
        outVertices.clear(); outElements.clear();
        for (auto &p : subNodes) outVertices.emplace_back(p);
        for (auto &e : subEdges) outElements.emplace_back(e.first, e.second);
        MeshIO::save("subgraph_" + outName, outVertices, outElements);
    }
}

size_t enumerateGraphs(size_t levels, size_t offset = 0, vector<size_t> taken = vector<size_t>()) {
    if (levels == 0) {
        processGraph(taken);
        return 1;
    }

    size_t count = 0;
    for (size_t i = offset; i < edges.size(); ++i) {
        taken.push_back(i);
        count += enumerateGraphs(levels - 1, i + 1, taken);
        taken.pop_back();
    }

    return count;
}

int main(int , char *[])
{
    nodesBary.reserve(15);

    BaryPoint p;
    p.fill(0);

    // Create the tet vertex nodes
    for (size_t i = 0; i < 4; ++i) {
        p[i] = 1;
        nodesBary.push_back(p);
        p[i] = 0;
    }

    // Create the tet edge nodes
    for (size_t i = 0; i < 4; ++i) {
        p[i] = 0.5;
        for (size_t j = i + 1; j < 4; ++j) {
            p[j] = 0.5;
            nodesBary.push_back(p);
            p[j] = 0;
        }
        p[i] = 0;
    }

    // Create the tet face nodes
    p.fill(1.0 / 3.0);
    for (size_t i = 0; i < 4; ++i) {
        p[i] = 0;
        nodesBary.push_back(p);
        p[i] = 1.0 / 3.0;
    }

    // Create the tet center node;
    p.fill(0.25);
    nodesBary.push_back(p);

    // Convert barycentric coordinates to spatial coords
    // (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1)
    nodes.reserve(nodesBary.size());
    for (const BaryPoint &p : nodesBary) {
        nodes.emplace_back(p[1] + p[2] + p[3],
                                  p[2] + p[3],
                                         p[3]);
    }

    // Create all graph edges.
    for (size_t i = 0; i < nodes.size(); ++i) {
        for (size_t j = i + 1; j < nodes.size(); ++j) {
            edges.push_back({i, j});
        }
    }

    vector<MeshIO::IOElement> tetEdges;
    vector<MeshIO::IOVertex> tetNodes;
    for (const auto &n : nodes)
        tetNodes.emplace_back(n);
    for (const auto &e : edges)
        tetEdges.emplace_back(e.first, e.second);
    MeshIO::save("tet.msh", tetNodes, tetEdges);

    enumerateGraphs(1);
    enumerateGraphs(2);
    enumerateGraphs(3);

    cout << "num nodes:\t" << nodes.size() << endl;
    cout << "num edges:\t" << edges.size() << endl;
    cout << "num patterns:\t" << numChecked << endl;

    cout << "num non-self intersecting base tet:\t" << nonSelfIntersecting << endl;
    cout << "connected:\t" << connectedCount << endl;
    cout << "no coinciding edges:\t" << noCoincidingEdges << endl;
    cout << "no dangling edges:\t" << noDanglingEdges << endl;
    cout << "valence less than 8:\t" << valenceTresholded << endl;

    cout << endl;
    cout << "valid patterns: " << validCount << endl;

    return 0;
}
