////////////////////////////////////////////////////////////////////////////////
// match.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Compare enumerated patterns to an existing set. Compare individual
//      patterns by lexicographically sorting.
//
//      Patterns in the first directory are searched for in the second
//      directory.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/12/2016 11:25:24
////////////////////////////////////////////////////////////////////////////////
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/Geometry.hh>
#include <MeshFEM/utils.hh>

#include <stdexcept>
#include <iostream>
#include <fstream>

using namespace std;

using Edge = std::pair<size_t, size_t>;

constexpr size_t NONE = numeric_limits<size_t>::max();

struct Vertex : public MeshIO::IOVertex {
    using Base = MeshIO::IOVertex;
    using Base::Base;

    bool operator<(const Vertex &b) const {
        for (size_t i = 0; i < 3; ++i) {
            if ((*this)[i] < b[i]) return true;
            if ((*this)[i] > b[i]) return false;
        }

        // Equal
        return false;
    }
};

int main(int argc, char *argv[])
{
    if (argc != 3) {
        cerr << "usage: ./match directory1.txt directory2.txt" << endl;
        exit(-1);
    }

    array<vector<vector<Vertex>>, 2> patternVertices;
    array<vector<vector<Edge>>, 2> patternEdges;
    array<vector<string>, 2> names;

    for (size_t d = 1; d <= 2; ++d) {
        ifstream dir(argv[d]);
        if (!dir) throw std::runtime_error("Couldn't open directory " + to_string(d));

        string path;
        while (getline(dir, path)) {
            vector<MeshIO::IOVertex> inVertices;
            vector<MeshIO::IOElement> elements;
            MeshIO::load(path, inVertices, elements);

            BBox<Point3D> bb(inVertices);
            auto badBBox = std::runtime_error("Unexpected bounding box");
            auto dim = bb.dimensions();
            if ((dim - Point3D(1, 1, 1)).norm() < 1e-8) {
                if ((bb.minCorner - Point3D(-0.5, -0.5, -0.5)).norm() > 1e-8) throw badBBox;
                if ((bb.maxCorner - Point3D( 0.5,  0.5,  0.5)).norm() > 1e-8) throw badBBox;
                // Scale to [-1, 1]
                for (auto &v : inVertices) v.point *= 2.0;
                bb = BBox<Point3D>(inVertices);
            }
            if ((bb.minCorner - Point3D(-1.0, -1.0, -1.0)).norm() > 1e-8) throw badBBox;
            if ((bb.maxCorner - Point3D( 1.0,  1.0,  1.0)).norm() > 1e-8) throw badBBox;

            vector<Vertex> vertices;
            vector<Edge> edges;

            vertices.reserve(inVertices.size());
            edges.reserve(elements.size());
            for (auto &v : inVertices) vertices.emplace_back(v.point);
            for (auto &e : elements) {
                if (e.size() != 2) throw runtime_error("Non-edge element encountered.");
                edges.push_back({e[0], e[1]});
            }

            // Determine lexicographic order
            // TODO: first round?
            vector<size_t> perm = sortPermutation(vertices);

            // Sort vertices lexicographically
            //      vertex that was labeled perm[i] is now labeled i
            vector<Vertex> sortedVertices;
            sortedVertices.reserve(vertices.size());
            for (size_t idx : perm)
                sortedVertices.push_back(vertices[idx]);

            // Invert permutation
            vector<size_t> iperm(perm.size());
            for (size_t i = 0; i < perm.size(); ++i)
                iperm[perm[i]] = i;

            // Remap edge indices. vertex that was labeled i is now labeled iperm[i]
            // Also orient edges in sorted order
            for (auto &e : edges) {
                size_t u = iperm[e.first ],
                       v = iperm[e.second];
                e.first  = min(u, v);
                e.second = max(u, v);
            }

            // Sort edges lexicographically
            sort(edges.begin(), edges.end());

            // // Output sorted graph (debug)
            // std::vector<MeshIO::IOVertex> outVertices;
            // std::vector<MeshIO::IOElement> outElements;
            // outVertices.clear(); outElements.clear();
            // for (auto &p : sortedVertices) outVertices.emplace_back(p);
            // for (auto &e : edges)          outElements.emplace_back(e.first, e.second);
            // MeshIO::save(path + ".sorted.obj", outVertices, outElements);

            patternVertices[d - 1].push_back(sortedVertices);
            patternEdges   [d - 1].push_back(edges);
            names          [d - 1].push_back(path);
        }
    }

    size_t dir1Size = patternVertices[0].size();
    size_t dir2Size = patternVertices[1].size();
    vector<bool> taken(patternVertices[1].size(), false);

    for (size_t i = 0; i < dir1Size; ++i) {
        const auto &v1 = patternVertices[0][i];
        const auto &e1 =    patternEdges[0][i];
        size_t matchIdx = NONE;
        for (size_t j = 0; j < dir2Size; ++j) {
            if (taken[j]) continue; // TODO: possibly check for non-injective map instead
            const auto &v2 = patternVertices[1][j];
            const auto &e2 =    patternEdges[1][j];
            if (v1.size() != v2.size()) continue;
            if (e1.size() != e2.size()) continue;
            bool isMatch = true;
            // Compare edges first (cheaper and exact)
            for (size_t ei = 0; ei < e1.size(); ++ei)
                if (e1[ei] != e2[ei]) { isMatch = false; break; }
            if (!isMatch) continue;
            // Then compare vertices
            for (size_t vi = 0; vi < v1.size(); ++vi)
                if ((v1[vi].point - v2[vi].point).squaredNorm() > 1e-8) { isMatch = false; break; }
            if (!isMatch) continue;

            matchIdx = j;
            break;
        }

        cout << names[0][i] << "\t";
        if (matchIdx == NONE) cout << "NONE";
        else cout << names[1][matchIdx];
        cout << endl;

        if (matchIdx != NONE) taken[matchIdx] = true;
    }

    return 0;
}
