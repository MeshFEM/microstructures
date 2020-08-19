////////////////////////////////////////////////////////////////////////////////
// SphereConvexHull.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes the signed distance to a convex hull of spheres.
//      The boundary of such a hull consists of triangles (tangent to three
//      spheres), portions of cones ("inflated edges"), and portions of spheres.
//
//      First, we determine this collection of boundary objects by computing the
//      convex hull of a triangulation of the spheres. Assuming the meshing is
//      fine enough to resolve the hull topology, we can determine these
//      boundary objects precisely.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  11/29/2016 17:02:41
////////////////////////////////////////////////////////////////////////////////
#ifndef SPHERECONVEXHULL_HH
#define SPHERECONVEXHULL_HH

#include "InflatorTypes.hh"
#include "SpherePoints.hh"
#include "TriangleClosestPoint.hh"
#include "ConvexHullTriangulation.hh"
#include "SignedDistance.hh"
#include "AutomaticDifferentiation.hh"

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/Geometry.hh>

#ifdef GROUND_TRUTH_DEBUGGING
#include "DisableWarnings.hh"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include "EnableWarnings.hh"
#include <MeshFEM/filters/gen_grid.hh>
#endif

#include <array>
#include <vector>
#include <map>
#include <set>
#include <limits>

namespace SignedDistance { namespace Primitives {

template<typename Real>
class SphereConvexHull {
    // "Oriented sphere triangle:" an ordered triplet of three sphere indices.
    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
public:
    SphereConvexHull(const std::vector<Point3<Real>> &centers,
                     const std::vector<Real> &radii)
        : m_sphereCenters(centers), m_sphereRadii(radii)
    {
#ifdef SPHEREHULLL_VERBOSE
        {
            std::cout << "SphereConvexHull on centers, radii:" << std::endl;
            for (const auto &pt : centers) {
                std::cout << "{"
                    << pt[0] << ", "
                    << pt[1] << ", "
                    << pt[2] << "}, ";
            }
            std::cout << std::endl;
            std::cout << "{";
            for (const Real &r : radii)
                std::cout << r << ", ";
            std::cout << "}" << std::endl;
        }
#endif

        // Determine the objects making up the convex hull:
        //      First generate points on each sphere and compute their convex hull
        std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>> spherePts;
        size_t pointsPerSphere = 1000;
        spherePts.reserve(centers.size() * pointsPerSphere);
        // Note: the sphere corresponding to point index i is
        // floor(i / pointsPerSphere)
        for (size_t i = 0; i < m_sphereCenters.size(); ++i) {
            generateSpherePointsWithExtremeVertices(pointsPerSphere, spherePts,
                    stripAutoDiff(m_sphereRadii[i]), stripAutoDiff(m_sphereCenters[i]));
            assert(spherePts.size() == pointsPerSphere * (i + 1));
        }

#ifdef WRITE_DEBUG_HULL_MESH
        std::vector<MeshIO::IOVertex > allVertices;
        std::vector<MeshIO::IOElement> allTriangles;

        std::vector<size_t> sphereIdx;
        for (size_t i = 0; i < centers.size(); ++i) {
            std::vector<Eigen::Matrix<double, 3, 1, 0, 3, 1>> ptSubset;
            const size_t offset = allVertices.size();

            for (size_t j = offset; j < std::min(offset + pointsPerSphere, spherePts.size()); ++j)
                ptSubset.push_back(spherePts[j]);

            std::vector<MeshIO::IOVertex > sphereVertices;
            std::vector<MeshIO::IOElement> sphereTriangles;
            convexHullFromTriangulation(ptSubset, sphereVertices, sphereTriangles);

            for (const auto &v : sphereVertices) {
                allVertices.push_back(v);
                sphereIdx.push_back(i);
            }
            for (auto e : sphereTriangles) {
                for (size_t &c : e) c += offset;
                allTriangles.emplace_back(e);
            }
        }

        ScalarField<Real> sphereIdxField(allVertices.size());
        for (size_t i = 0; i < allVertices.size(); ++i)
            sphereIdxField[i] = sphereIdx[i];

        MSHFieldWriter writer("sphere_debug.msh", allVertices, allTriangles);
        writer.addField("sphereIdx", sphereIdxField, DomainType::PER_NODE);
#endif

        std::vector<MeshIO::IOVertex > hullVertices;
        std::vector<MeshIO::IOElement> hullTriangles;
        std::vector<size_t> originatingVtx;
        convexHullFromTriangulation(spherePts, hullVertices, hullTriangles, originatingVtx);

#ifdef WRITE_DEBUG_HULL_MESH
        MeshIO::save("hull_debug.msh", hullVertices, hullTriangles);
#endif

        auto sphereForVtx = [&](size_t i) { return originatingVtx[i] / pointsPerSphere; };

        // Determine sphere connectivity from convex hull triangles.
        // Note: can change to std::set here for hulls with large numbers of
        // "supporting sphere triangles." std::vector should be faster for the
        // typical case of only a few triangles.
        std::vector<OrientedTriplet> sphereTris;
        for (const auto &t : hullTriangles) {
            OrientedTriplet st(sphereForVtx(t[0]),
                               sphereForVtx(t[1]),
                               sphereForVtx(t[2]));
            auto it = std::find(sphereTris.begin(), sphereTris.end(), st);
            if (it == sphereTris.end())
                sphereTris.push_back(st);
        }

        // Determine the non-degenerate triangles and the edge- and
        // point-degenerated triangles.
        // Since degenerate sphere triangles whose centers are not degenerate
        // are most easily detected by attempting to construct a tangent plane,
        // we compute the tangent planes here.
        std::set<size_t> degeneratedPt;
        for (const auto &st : sphereTris) {
            // "Combinatorial dimension"
            size_t dim = st.dimension();
            if (dim == 2) {
                // The true geometric/numeric dimension can still be 1 or 0
                // depending on the sphere centers. This happens, e.g., for
                // three collinear spheres of the same radius, when sample points
                // for the middle sphere happen to lie on the convex hull.
                const auto &p1 = m_sphereCenters[st.corners[0]],
                           &p2 = m_sphereCenters[st.corners[1]],
                           &p3 = m_sphereCenters[st.corners[2]];
                //   3
                // 1   2
                Vector3<Real> e[3] = { p3 - p2,
                                       p1 - p3,
                                       p2 - p1 };
                Real lSq[3] = { e[0].squaredNorm(),
                                e[1].squaredNorm(),
                                e[2].squaredNorm() };

                Real sinSq = e[0].cross(e[1]).squaredNorm() / (lSq[0] * lSq[1]);
                if (sinSq < 1e-10) {
                    // Triangle degenerates to longest edge (opposite vtx i) or
                    // to a point.
                    int i = std::distance(lSq, std::max_element(lSq, lSq + 3));
                    // Note: this can create numerically degenerate edges.
                    if (lSq[i] > 1e-10)
                        m_sphereChainEdges.emplace(st.corners[(i + 1) % 3],
                                                   st.corners[(i + 2) % 3]);
                    continue;
                }

                // Now check for the degenerate cases where there is no plane
                // tangent to the three spheres. (Ideally this triangle wouldn't
                // be returned by CGAL, but it happens in practice).
                // We do this by attempting to compute the spheres' tangent
                // plane.
                const Real &r1 = m_sphereRadii[st.corners[0]],
                           &r2 = m_sphereRadii[st.corners[1]],
                           &r3 = m_sphereRadii[st.corners[2]];
                // Note: different edge numbering here from above;
                // TODO: relabel consistently
                //   3
                // 1   2
                Vector3<Real> e1 =  e[2], // p2 - p1 =  e[2]
                              e2 = -e[1]; // p3 - p1 = -e[1]
                Real l1 = sqrt(lSq[2]),
                     l2 = sqrt(lSq[1]);
                e1 /= l1;
                e2 /= l2;

                Vector3<Real> midplaneNormal = e1.cross(e2);
                midplaneNormal /= midplaneNormal.norm();
                Vector3<Real> e1perp = midplaneNormal.cross(e1);
                // This degeneracy should already have been caught.
                if (std::abs(stripAutoDiff(e1perp.norm() - 1.0)) > 1e-10)
                    report_error("Degenerate convex hull triangle uncaught (e1perp)", __LINE__);

                // Determine the (positive) tangent plane for the three spheres.
                // We compute the tangent plane's normal in terms of its components
                // along e1, e1^perp, and midplane normal n_mp:
                //  n = alpha e1 + beta e1perp + gamma n_mp
                // Here e1, e2, e1perp, and n_mp are unit vectors
                Real alpha = (r1 - r2) / l1;
                Real beta = ((r1 - r3) / l2 - alpha * e2.dot(e1)) / e2.dot(e1perp);
                Real aSqbSq = alpha * alpha + beta * beta;
                if (aSqbSq < 1) {
                    Real gamma = sqrt(1 - aSqbSq);
                    Vector3<Real> n = alpha * e1 + beta * e1perp + gamma * midplaneNormal;
                    m_supportingTriangles.emplace_back((p1 + r1 * n).eval(),
                                                       (p2 + r2 * n).eval(),
                                                       (p3 + r3 * n).eval());
                    m_sphereTriangles.push_back(st);
                }
                else {
                    std::cerr << "aSqbSq: " << aSqbSq << std::endl;
                    std::cerr << "Degenerate sphere triangle configuration; triangle ignored." << std::endl;
                }
            }
            // Note: even though (oriented) triangles are unique, a triangle
            // that has degenerated into an edge can be repeated. E.g.:
            // {{1, 1, 2}, {1, 2, 2}, {1, 2, 1}}
            // Thus, we use a set to remove duplicates.
            // Also note that this can produce numerically degenerate edges.
            if (dim == 1) m_sphereChainEdges.insert(st.degenerateAsEdge());
            if (dim == 0) { degeneratedPt.insert(st.corners[0]); }
        }

        // The edge chains consist only of the edges that do not belong to
        // non-degenerate triangles--remove the triangle edges.
        // Also, detect and remove numerically degenerate edges.
        // TODO: Detect and remove double-triangle edge created at the top in the
        // following configuration (this erroneously becomes a chain edge):
        //           +--+--+
        //            \ | /
        //             \|/
        //              *
        for (auto it = m_sphereChainEdges.begin(); it != m_sphereChainEdges.end(); /* advanced inside */) {
            bool isTriEdge = false;
            for (const auto &st : m_sphereTriangles)
                if (st.containsEdge(*it)) { isTriEdge = true; break; }

            bool isDegenerate = (m_sphereCenters[(*it)[0]] - m_sphereCenters[(*it)[1]]).squaredNorm() < 1e-16;
            if (isDegenerate) {
                assert(std::abs(stripAutoDiff(m_sphereRadii[(*it)[0]] - m_sphereRadii[(*it)[1]])) < 1e-8);
                degeneratedPt.insert((*it)[0]);
            }
            bool remove = isTriEdge || isDegenerate;
            if (remove) it = m_sphereChainEdges.erase(it);
            else      ++it;
        }

#ifdef SPHEREHULLL_VERBOSE
        writeDebugInfo();
#endif

        // Check if the convex hull consists of a single sphere.
        if ((m_sphereTriangles.size() == 0) && (m_sphereChainEdges.size() == 0)) {
            // We could have the case of combinatorially non-degenerate spheres but
            // numerically degenerate.
            // TODO: validate this is the case.
            assert(degeneratedPt.size() > 0);
            m_degenerateHullSphereIdx = *degeneratedPt.begin();
        }

        // Create inflated geometry for every edge on the hull.
        for (const auto &st : m_sphereTriangles) {
            for (size_t i = 0; i < 3; ++i) {
                UnorderedPair e(st.corners[i], st.corners[(i + 1) % 3]);
                if (m_inflatedEdges.count(e)) continue;
                m_inflatedEdges.emplace(e,
                    InflatedEdge<Real>(m_sphereCenters.at(e[0]), m_sphereCenters.at(e[1]),
                                       m_sphereRadii.at(e[0]), m_sphereRadii.at(e[1])));
            }
        }
        for (const auto &e : m_sphereChainEdges) {
            // Sphere chain edges are unique and should be distinct from those
            // inflated previously.
            assert(m_inflatedEdges.count(e) == 0);
            auto result = m_inflatedEdges.emplace(e,
                InflatedEdge<Real>(m_sphereCenters.at(e[0]), m_sphereCenters.at(e[1]),
                                   m_sphereRadii.at(e[0]), m_sphereRadii.at(e[1])));
            assert(result.second); // edges are unique; should always be inserted
            m_sphereChainInflatedEdges.push_back(&(result.first->second));
        }

        // Create pointers to the inflated edge opposite each triangle corner.
        m_edgeOppositeTriangleCorner.reserve(3 * numTriangles());
        for (const auto &st : m_sphereTriangles) {
            for (size_t i = 0; i < 3; ++i) {
                size_t nextVertex = st.corners[(i + 1) % 3];
                size_t prevVertex = st.corners[(i + 2) % 3];
                UnorderedPair e_opp(nextVertex, prevVertex);
                m_edgeOppositeTriangleCorner.push_back(&m_inflatedEdges.at(e_opp));
            }
        }
    }

    void report_error(const std::string &msg, size_t line) const {
        std::cerr.precision(19);
        std::cerr << "Error in SphereConvexHull for centers, radii:" << std::endl;
        for (const auto &pt : m_sphereCenters) {
            std::cerr << "{"
                << pt[0] << ", "
                << pt[1] << ", "
                << pt[2] << "}, ";
        }
        std::cerr << std::endl;
        std::cerr << "{";
        for (const Real &r : m_sphereRadii)
            std::cerr << r << ", ";
        std::cerr << "}" << std::endl;

        writeDebugInfo();

        throw std::runtime_error(msg + " at SphereConvexHull.hh:" + std::to_string(line));
    }

    bool hullIsSingleSphere() const { return m_degenerateHullSphereIdx == NONE; }

    void writeDebugInfo() const {
        std::cerr << numTriangles() << " triangles:" << std::endl;
        for (const auto &st : m_sphereTriangles) {
            std::cerr << "{" << st.corners[0]
                << ", " << st.corners[1]
                << ", " << st.corners[2]
                << "}" << std::endl;
        }
        std::cerr << numChainEdges() << " edges:" << std::endl;
        for (const auto &se : m_sphereChainEdges) {
            std::cerr << "{" << se[0]
                << ", " << se[1]
                << "}" << std::endl;
        }
    }

    template<typename Real2, bool DebugOutput = false>
    Real2 signedDistance(const Point3<Real2> &p) const {
        if (DebugOutput)  {
            std::cerr.precision(19);
            std::cerr << "querying distance at point " << stripAutoDiff(p).transpose() << std::endl;
            {
                std::cout << "to SphereConvexHull of centers, radii:" << std::endl;
                for (const auto &pt : m_sphereCenters) {
                    std::cout << "{"
                        << stripAutoDiff(pt[0]) << ", "
                        << stripAutoDiff(pt[1]) << ", "
                        << stripAutoDiff(pt[2]) << "}, ";
                }
                std::cout << std::endl;
                std::cout << "{";
                for (const Real &r : m_sphereRadii)
                    std::cout << stripAutoDiff(r) << ", ";
                std::cout << "}" << std::endl;
            }
        }

        // Fully degenerate case (hull is a single sphere)
        if (m_degenerateHullSphereIdx != NONE) {
            // Auto-diff compatible version:
            return (p - m_sphereCenters[m_degenerateHullSphereIdx].template cast<Real2>()).norm() -
                    m_sphereRadii[m_degenerateHullSphereIdx];
        }

        // Find distance the non-degenerate (volume) regions (if any exist)
        Real2 sd = safe_numeric_limits<Real2>::max();
        if (m_supportingTriangles.size() > 0) {
            Real2 closestTriDistSq = safe_numeric_limits<Real2>::max();
            Vector3<Real2> closestLambda, closestInternalLambda;
            closestLambda.setZero(); closestInternalLambda.setZero(); // silence GCC's warnings
            size_t closestTri = 0;
            for (size_t i = 0; i < m_supportingTriangles.size(); ++i) {
                Vector3<Real2> lambda = m_supportingTriangles[i].baryCoords(p);
                Vector3<Real2> internalLambda = m_supportingTriangles[i].closestInternalBarycoordsToBaryCoords(lambda);
                Real2 distSq = (m_supportingTriangles[i].pointAtBarycoords(internalLambda) - p).squaredNorm();
                if (distSq < closestTriDistSq) {
                    closestTriDistSq = distSq;
                    closestLambda = lambda;
                    closestInternalLambda = internalLambda;
                    closestTri = i;
                }
            }

            // If the projection of the point onto the closest triangle's plane
            // is extremely close to the triangle, treat the point as internal
            // to the triangle. This fixes the corner case where > 3 coplanar
            // spheres lead to a triangulated tangent plane with > 1 triangle,
            // which in turn suffered from incorrect distance on the slice above
            // the shared edge (the code assumed the closest points in this
            // slice lie on a conical frustum patch between the two triangles,
            // but none exists).
            for (size_t i = 0; i < 3; ++i) {
                if ((closestLambda[i] <= 0) && (closestLambda[i] > -1e-12))
                    closestInternalLambda[i] = 1e-16; // perturb inside
            }

            const auto &cst = m_supportingTriangles[closestTri];
            int numZero = 0;
            int zeroCoords[3];
            if (closestInternalLambda[0] == 0.0) { zeroCoords[numZero++] = 0; }
            if (closestInternalLambda[1] == 0.0) { zeroCoords[numZero++] = 1; }
            if (closestInternalLambda[2] == 0.0) { zeroCoords[numZero++] = 2; }
            if (numZero == 3) {
                // this has actually triggered before!!!
                std::cerr << std::setprecision(19) << std::endl;
                std::cerr << "closestTri = " << closestTri << std::endl;
                std::cerr << "closestInternalLambda = " << stripAutoDiff(closestInternalLambda).transpose() << std::endl;
                std::cerr << "p = " << stripAutoDiff(p).transpose() << std::endl;
                std::cerr << "closestTriDist = " << sqrt(closestTriDistSq) << std::endl;
                std::cerr << "closest tri vertices:" << std::endl;
                std::cerr << stripAutoDiff(m_supportingTriangles[closestTri].p0()).transpose() << std::endl;
                std::cerr << stripAutoDiff(m_supportingTriangles[closestTri].p1()).transpose() << std::endl;
                std::cerr << stripAutoDiff(m_supportingTriangles[closestTri].p2()).transpose() << std::endl;

                std::cerr << "Sphere centers: " << std::endl;
                for (auto &c : m_sphereCenters)
                    std::cerr << "{" << c[0] << ", " << c[1] << ", " << c[2] << "}" << std::endl;
                std::cerr << "Sphere radii: {" << std::endl;
                for (Real r : m_sphereRadii)
                    std::cerr << r << ", ";
                std::cerr << "}" << std::endl;
                report_error("numZero == 3", __LINE__);
            }

            // If closest point is in the triangle's interior, it's the closest
            // hull point.
            if (numZero == 0) {
                // Compute signed distance in a way that's differentiable on the
                // triangle. (sqrt(closestTriDistSq) is not!)
                auto v = (p - cst.pointAtBarycoords(closestInternalLambda)).eval();
                const auto &n = cst.normal();
                // dot() isn't working on autodiff types.
                sd = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
                // assert(std::abs(std::abs(stripAutoDiff(sd)) - sqrt(stripAutoDiff(closestTriDistSq))) < 1e-9);

                // if (normalDist >= 0) sd =  sqrt(closestTriDistSq);
                // else                 sd = -sqrt(closestTriDistSq);
                if (DebugOutput) {
                    std::cerr << "Closest hull point in triangle interior" << std::endl;
                    std::cerr << " closestTriDistSq ( " << closestTriDistSq << "):"; reportDerivatives(std::cerr, closestTriDistSq); std::cerr << std::endl;
                }
            }

            if (numZero == 1) {
                if (DebugOutput) { std::cerr << "Closest hull point on triangle edge" << std::endl; }
                // If exactly one barycentric coordinate is zero (i.e. we're on
                // the supporting triangle's edge), the closest surface point to
                // p lies on the convex hull of the two endpoint spheres.
                sd = m_edgeOppositeTriangleCorner[3 * closestTri + zeroCoords[0]]->signedDistance(p);
            }
            if (numZero == 2) {
                if (DebugOutput) { std::cerr << "Closest hull point on triangle edge pair" << std::endl; }
                // If two barycentric coordinates are zero (we're on a vertex),
                // the closest surface point p lies on the union of the two
                // incident inflated edges. These edges are the ones across from
                // the vertices with barycentric coordinate zero.
                // It turns out that std::min computes the correct signed
                // distance in all cases:

                // std::min computes the correct signed distance iff the true
                // closest surface point is the closest point on at least one of
                // the input primitive surfaces. This is always the case for p
                // outside. For p inside, the computed signed distance can be
                // higher (closer to zero) than the true distance--this happens
                // when the closest primitive points are false boundary points
                // (each lies inside another object). In the case of our edge
                // joint geometry, this only happens where a sharp concave
                // corner is formed by the union.
                // Thus std::min will give the correct result for "convex edge
                // pairs." Edge pairs forming the border of "sphere triangles"
                // in the convex hull are **not** convex, but the concave parts
                // should always be covered up by the triangles (i.e. query
                // points closest to these concave parts actually fall into the
                // "numZero == 0" case: closest hull point is on a triangle).
                sd = std::min(m_edgeOppositeTriangleCorner[3 * closestTri + zeroCoords[0]]->signedDistance(p),
                              m_edgeOppositeTriangleCorner[3 * closestTri + zeroCoords[1]]->signedDistance(p));
            }
        }

        if (DebugOutput) { std::cerr << "Pre-chain sd " << sd << " derivatives:"; reportDerivatives(std::cerr, sd); std::cerr << std::endl; }

        // Union in the edge chain geometry.
        // Inflated edge chains must always union into convex geometry and
        // thus the signed distance can be computed accurately with
        // std::min (see numZero == 2 discussion). Further, the edge chains attach to the non-degenerate
        // regions in a convex way, so std::min can be used to union
        // everything together.
        for (const auto edge_ptr : m_sphereChainInflatedEdges) {
            Real2 edge_sd = edge_ptr->template signedDistance<Real2, DebugOutput>(p);
            if (DebugOutput) { std::cerr << "edge_sd " << edge_sd << " derivatives:"; reportDerivatives(std::cerr, edge_sd); std::cerr << std::endl; }
            sd = std::min<Real2>(sd, edge_sd);
        }

        if (DebugOutput) { std::cerr << "Post-chain sd " << sd << " derivatives:"; reportDerivatives(std::cerr, sd); std::cerr << std::endl; }

        return sd;
    }

#ifdef GROUND_TRUTH_DEBUGGING
    // For debugging: output the "ground truth" (approximate) signed distance
    // field computed by meshing the convex hull and using CGAL's AABB tree.
    void writeGroundTruth(size_t sphereResolution, size_t gridResolution, const std::string &path) const {
        // Determine the objects making up the convex hull:
        //      First generate points on each sphere and compute their convex hull
        std::vector<Point3<double>> spherePts;
        for (size_t i = 0; i < m_sphereCenters.size(); ++i) {
            Point3<double> center(m_sphereCenters[i].template cast<double>());
            generateSpherePoints(sphereResolution, spherePts,
                    stripAutoDiff(m_sphereRadii[i]), center);
        }

        std::vector<MeshIO::IOVertex > hullVertices;
        std::vector<MeshIO::IOElement> hullTriangles;
        std::vector<size_t> originatingVtx;
        convexHullFromTriangulation(spherePts, hullVertices, hullTriangles, originatingVtx);

        using K                    = CGAL::Simple_cartesian<double>;
        using Point                = K::Point_3;
        using Triangle             = K::Triangle_3;
        using TriIterator          = std::vector<Triangle>::iterator;
        using AABB_triangle_traits = CGAL::AABB_traits<K, CGAL::AABB_triangle_primitive<K, TriIterator>>;
        using Tree                 = CGAL::AABB_tree<AABB_triangle_traits>;

        std::vector<Triangle> cgal_triangles;
        cgal_triangles.reserve(hullTriangles.size());
        for (const auto &ht : hullTriangles) {
            auto p1 = hullVertices[ht[0]];
            auto p2 = hullVertices[ht[1]];
            auto p3 = hullVertices[ht[2]];
            cgal_triangles.emplace_back(Point(p1[0], p1[1], p1[2]),
                                        Point(p2[0], p2[1], p2[2]),
                                        Point(p3[0], p3[1], p3[2]));
        }
        // construct AABB tree
        Tree tree(cgal_triangles.begin(), cgal_triangles.end());
        tree.accelerate_distance_queries();
        // Compute and output distances using the AABB tree
        {
            std::vector<MeshIO::IOVertex > gv;
            std::vector<MeshIO::IOElement> ge;
            std::vector<size_t> grid_sizes({gridResolution, gridResolution, gridResolution});
            gen_grid(grid_sizes, gv, ge);
            ScalarField<Real> sd(gv.size());
            for (size_t i = 0; i < gv.size(); ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    gv[i].point[j] *= 2.0 / gridResolution;
                    gv[i].point[j] -= 1.0;
                }
                Point samplePt(gv[i].point[0], gv[i].point[1], gv[i].point[2]);
                Point cgal_closestPt;
                TriIterator closest_tri;
                std::tie(cgal_closestPt, closest_tri) = tree.closest_point_and_primitive(samplePt);

                Point3<Real> closestPt(cgal_closestPt.x(), cgal_closestPt.y(), cgal_closestPt.z());
                Real dist = (gv[i].point - closestPt).norm();

                if (closest_tri->supporting_plane().oriented_side(samplePt) == CGAL::ON_NEGATIVE_SIDE)
                    dist *= -1;
                sd[i] = dist;
            }
            MSHFieldWriter writer(path, gv, ge);
            writer.addField("sd", sd, DomainType::PER_NODE);
        }
    }
#endif

    size_t numTriangles()  const { return m_sphereTriangles.size(); }
    size_t numChainEdges() const { return m_sphereChainEdges.size(); }

    const std::vector<Point3<Real>> &sphereCenters() const { return m_sphereCenters; }
    const std::vector<Real>         &sphereRadii()   const { return   m_sphereRadii; }

private:
    // position and radii of all spheres; the triangles/edges/points on the hull
    // will index into these.
    std::vector<Point3<Real>> m_sphereCenters;
    std::vector<Real>         m_sphereRadii;
    // Triplets of spheres whose tangent triangles are the hull's supporting triangles.
    std::vector<OrientedTriplet> m_sphereTriangles;
    // Edges making up the edge chain(s) (edges that do not form triangles)
    std::set<UnorderedPair> m_sphereChainEdges;

    // Signed distance function for all conical frustum (edge) geometry
    // appearing on the hull--including edges of the triangles.
    std::map<UnorderedPair, InflatedEdge<Real>> m_inflatedEdges;

    // Pointer to the inflated edge opposite each of the 3 * numTriangles() triangle corners
    std::vector<const InflatedEdge<Real> *> m_edgeOppositeTriangleCorner;
    // Pointer to the inflated edge for each chain edge
    std::vector<const InflatedEdge<Real> *> m_sphereChainInflatedEdges;

    std::vector<TriangleClosestPoint<Real>> m_supportingTriangles;

    size_t m_degenerateHullSphereIdx = NONE;
};

}} // close namespace SD::Primitives

#endif /* end of include guard: SPHERECONVEXHULL_HH */
