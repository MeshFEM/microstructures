#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>

#include "intersection_check.hh"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
using namespace std;

// Check all non-neighboring edges for intersections
bool hasSelfIntersection(const vector<Point3d> &nodes, const vector<pair<size_t, size_t>> &edges) {
    for (size_t i = 0; i < edges.size(); ++i) {
        const auto &ei = edges[i];
        const auto &vi0 = nodes[ei.first];
        const auto &vi1 = nodes[ei.second];

        K::Segment_3 segi(K::Point_3(vi0[0], vi0[1], vi0[2]),
                          K::Point_3(vi1[0], vi1[1], vi1[2]));
        for (size_t j = i + 1; j < edges.size(); ++j) {
            const auto &ej = edges[j];
            const auto &vj0 = nodes[ej.first];
            const auto &vj1 = nodes[ej.second];
            K::Segment_3 segj(K::Point_3(vj0[0], vj0[1], vj0[2]),
                              K::Point_3(vj1[0], vj1[1], vj1[2]));
            // Ignore neighbors
            if ((ei.first  == ej.first) || (ei.first  == ej.second) ||
                (ei.second == ej.first) || (ei.second == ej.second)) continue;
            // cout << "Testing " << vi0 << endl << vi1 << endl << " against " << vj0 << endl << vj1 << endl;

            if (squared_distance(segi, segj) < 1e-8) return true;
        }
    }

    return false;
}
