////////////////////////////////////////////////////////////////////////////////
// SnapAndReflect.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Snap representative base cell geometry to exactly fill [0, 1]^3, then
//      reflect into [-1, 1]^3 via the orthotropic symmetry planes.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/30/2015 14:37:00
////////////////////////////////////////////////////////////////////////////////
#ifndef SNAPANDREFLECT_HH
#define SNAPANDREFLECT_HH

#include "Isometries.hh"
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/PeriodicBoundaryMatcher.hh>
#include <MeshFEM/Concepts.hh>
#include <vector>
#include <queue>
#include <ratio>
#include <stdexcept>

#define SNAP_DEBUG 0

inline bool isEq(Real a, Real b, Real tol = 0) {
    return std::abs(a - b) < tol;
}

// TODO: choose tolerance that works with both 2D and 3D inflator?
template<typename Vertex, typename TOL = std::ratio<2, long(1e4)>>
void snapVerticesToUnitCell(std::vector<Vertex> &vertices,
                            std::vector<std::vector<bool>> &onMinFace,
                            std::vector<std::vector<bool>> &onMaxFace) {
    static constexpr double tolerance = double(TOL::num) / double(TOL::den);
    // Snap vertices to [0, 1], determining min/max face membership
    onMinFace.assign(3, std::vector<bool>(vertices.size(), false));
    onMaxFace.assign(3, std::vector<bool>(vertices.size(), false));

    for (size_t vi = 0; vi < vertices.size(); ++vi) {
        auto &v = vertices[vi];
        for (size_t c = 0; c < 3; ++c) {
            if (isEq(v[c], 0, tolerance)) { v[c] = 0; onMinFace[c][vi] = true; }
            if (isEq(v[c], 1, tolerance)) { v[c] = 1; onMaxFace[c][vi] = true; }
        }
    }
}

// The 3D inflator is not guaranteed to place vertices precisely on the base
// cell faces. However, with the CGALClippedVolumeMesher, vertices at the
// intersection of the surface with the base cell will be placed exactly
// (they are extracted with marching squares on the cell faces).
// We can take advantage of this to segment the surface into components
// separated by base cell intersection vertices. Then we can detect and snap
// all vertices in the components lying approximately on the cell faces at
// once, which should be more robust than applying a threshold to each vertex
// independently.
// @param[inout] vertices   mesh vertices to snap
// @param[in]    m          tet mesh over "vertices"
// @param[in]    epsilon    threshold for the distance of a component's
//                          vertices to the cell boundary; if all distances of
//                          component vertices are <= this threshold, the
//                          entire component is snapped to the boundary.
template<class TMesh>
enable_if_not_models_concept_t<Concepts::TetMesh, TMesh, void>
smartSnap3D(std::vector<MeshIO::IOVertex> &/* vertices */, const TMesh &/* mesh */, const BBox<Point3D> &/* cell */,
                 const Real /* epsilon */ = 1e-4) {
    throw std::runtime_error("smartSnap3D must be called on a tet mesh!");
}

template<class TMesh>
enable_if_models_concept_t<Concepts::TetMesh, TMesh, void>
smartSnap3D(std::vector<MeshIO::IOVertex> &vertices, const TMesh &mesh, const BBox<Point3D> &cell,
                 const Real epsilon = 1e-3) {
    // std::cout << "Snapping to meshing cell " << cell << std::endl;
    // Partition the boundary faces into components separated by vertices on the
    // base cell.
    const size_t nbf = mesh.numBoundaryFaces();

    // First mark vertices making up the border between components.
    using FM = PeriodicBoundaryMatcher::FaceMembership<3>;
    std::vector<FM> faceMembership;
    for (auto bv : mesh.boundaryVertices())
        faceMembership.emplace_back(vertices.at(bv.volumeVertex().index()), cell, 0); // zero tolerance!

    auto isSegmentBorderEdge = [&](typename TMesh::template BHEHandle<const TMesh> be) {
        auto  tipFM = faceMembership.at(be. tip().index()),
             tailFM = faceMembership.at(be.tail().index());
        // Treat edges lying perfectly on base cell faces as segment borders.
        // Note: this can cause more segmentation than needed on the boundary
        // faces, but this is fine; we need only the interior to be segmented
        // properly to prevent excessive snapping.
        return (tipFM & tailFM).onAnyFace();
    };

    // Run BFS to partition the boundary faces into components separated by
    // segment border edges. This traversal works on non-manifold meshes (each
    // incident surface patch on a non-manifold boundary edge/vertex will be
    // traversed separately).
    std::vector<size_t> component(nbf, 0);
    std::queue<size_t> bfsQueue;
    size_t componentIdx = 0;
    for (size_t bfi = 0; bfi < nbf; ++bfi) {
        if (component[bfi] != 0) continue;
        component[bfi] = ++componentIdx;
        bfsQueue.push(bfi);
        while (!bfsQueue.empty()) {
            size_t u = bfsQueue.front();
            bfsQueue.pop();
            for (auto bhe_u : mesh.boundaryFace(u).halfEdges()) {
                if (isSegmentBorderEdge(bhe_u)) continue;
                size_t v = bhe_u.opposite().face().index();
                if (component.at(v) == 0) {
                    component[v] = componentIdx;
                    bfsQueue.push(v);
                }
                assert(component[v] == componentIdx);
            }
        }
    }

    // 0th component is unassigned...
    const size_t nComponents = componentIdx + 1;

#if SNAP_DEBUG
    {
        const size_t nbv = mesh.numBoundaryVertices();
        // Output surface mesh marked with connected components for debugging
        ScalarField<Real> segmentBorder(nbv);
        segmentBorder.clear();
        for (auto be : mesh.boundaryHalfEdges()) {
            if (isSegmentBorderEdge(be)) {
                segmentBorder(be.tip().index())  += 1.0;
                segmentBorder(be.tail().index()) += 1.0;
            }
        }

        ScalarField<Real> componentIndicator(nbf);
        for (size_t bfi = 0; bfi < nbf; ++bfi)
            componentIndicator[bfi] = component[bfi];

        std::vector<MeshIO::IOVertex > bVerts;
        std::vector<MeshIO::IOElement> bElems;
        for (auto bv : mesh.boundaryVertices())
            bVerts.push_back(vertices.at(bv.volumeVertex().index()));
        for (auto bs : mesh.boundarySimplices())
            bElems.emplace_back(bs.vertex(0).index(), bs.vertex(1).index(), bs.vertex(2).index());

        MSHFieldWriter writer("cellfaceComponents.msh", bVerts, bElems);
        writer.addField("component indicator", componentIndicator, DomainType::PER_ELEMENT);
        writer.addField("segment border", segmentBorder, DomainType::PER_NODE);

        for (size_t i = 0; i < 6; ++i) {
            ScalarField<Real> bvFM(nbv);
            for (auto bv : mesh.boundaryVertices())
                bvFM(bv.index()) = faceMembership.at(bv.index()).membership[i];
            writer.addField("bvFM " + std::to_string(i), bvFM, DomainType::PER_NODE);
        }

    }
#endif

    // Maximum over component vertices of distance to the 6 cell faces:
    //    min x, min y, min z, max x, max y, max z
    std::vector<std::array<Real, 6>> componentMaxDistances(nComponents, { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0} });
    for (auto bf : mesh.boundaryFaces()) {
        std::array<Real, 6> &cmd = componentMaxDistances.at(component[bf.index()]);
        for (auto bv : bf.vertices()) {
            const auto &v = vertices.at(bv.volumeVertex().index());
            // Compute dist from v to the min/max x/y/z cell faces
            for (size_t c = 0; c < 3; ++c) {
                cmd[    c] = std::max(cmd[    c], std::abs(v[c] - cell.minCorner[c]));
                cmd[3 + c] = std::max(cmd[3 + c], std::abs(v[c] - cell.maxCorner[c]));
            }
        }
    }

    const size_t NONE = std::numeric_limits<size_t>::max();
    std::vector<size_t> cellFaceForComponent(nComponents, NONE);
    for (size_t ci = 1; ci < nComponents; ++ci) {
        size_t count = 0;
        for (size_t f = 0; f < 6; ++f) {
            if (componentMaxDistances[ci][f] <= epsilon) {
                cellFaceForComponent[ci] = f;
                ++count;
            }
        }
        if (count > 1)
            std::cerr << "WARNING: ambiguous component" << std::endl;
    }

    // Snap component faces to their respective base cell faces.
    for (auto bf : mesh.boundaryFaces()) {
        size_t cf = cellFaceForComponent.at(component[bf.index()]);
        if (cf == NONE) continue;
        assert(cf < 6);
        bool minFace = cf < 3;
        if (!minFace) cf -= 3;
        Real snapVal = minFace ? cell.minCorner[cf] : cell.maxCorner[cf];

        for (auto bv : bf.vertices()) {
            auto &vtx = vertices.at(bv.volumeVertex().index());
            assert(std::abs(vtx.point[cf] - snapVal) < epsilon);
            vtx[cf] = snapVal;
        }
    }

    BBox<Point3D> snappedBB(vertices);
    if (snappedBB != cell) throw std::runtime_error("Snapped mesh bounding box does not equal meshing cell.");
}

// For when smartSnap3D fails...
// Tolerance "epsilon" should be based on CGAL meshing parametrs
inline void dumbSnap3D(std::vector<MeshIO::IOVertex> &vertices,
                const BBox<Point3D> &cell, const Real epsilon) {
    for (auto &v : vertices) {
        for (size_t c = 0; c < 3; ++c) {
            if (isEq(v[c], cell.minCorner[c], epsilon)) v[c] = cell.minCorner[c];
            if (isEq(v[c], cell.maxCorner[c], epsilon)) v[c] = cell.maxCorner[c];
        }
    }
}

// Generate full reflected cell in three steps: reflect along x axis,
// reflect along y axis, then reflect along z
//     +------+------+
//    /  reflect z  /|
//   +------+------+ |
//  /      /      /| |
// +------+------+ | +
// |      |      | |/|
// | refx | base | + |
// |      | (0)  |/| |
// +------+------+ | |
// |      |      | |/
// |  reflect y  | |
// |      |      |/
// +------+------+
// Vertices on the reflection planes must not be duplicated, and the
// onReflectionPlane arrays allow us to enforce this.
// onReflectionPlane[d] stores whether vertices lie on reflection plane d.
// We update these arrays with each reflection to keep them valid.
// If Dim = 2, only X and Y reflection are performed.
template<typename Vertex, typename Element>
void reflectXYZ(size_t Dim, // Dimensions to reflect in (length of [x, y, z] prefix)
                const std::vector<Vertex> &vertices,
                const std::vector<Element> &elements,
                std::vector<std::vector<bool>> onReflectionPlane, // copy; changed inside
                std::vector<Vertex>  &reflectedVertices,
                std::vector<Element> &reflectedElements,
                std::vector<size_t>   &vertexOrigin,
                std::vector<Isometry> &vertexIsometry)
{
    assert(onReflectionPlane.size() == 3);
    assert(onReflectionPlane[0].size() == vertices.size() &&
           onReflectionPlane[1].size() == vertices.size() &&
           onReflectionPlane[2].size() == vertices.size());

    reflectedVertices = vertices;
    reflectedElements = elements;

    // We start with the original vertices, so origins/isometries are all the
    // identity.
    vertexOrigin.assign(vertices.size(), 0);
    for (size_t i = 0; i < vertexOrigin.size(); ++i) vertexOrigin[i] = i;
    vertexIsometry.assign(vertices.size(), Isometry());

    for (size_t d = 0; d < Dim; ++d) {
        auto refl = Isometry::reflection(static_cast<Symmetry::Axis>(d));
        // We need a mapping from vertex indices of the new reflected geometry
        // we're about to create to global vertex indices.
        // All vertices except those on the reflection pane are copied.
        std::vector<size_t> globalVertexIndex(reflectedVertices.size());
        size_t numVertices = reflectedVertices.size();
        for (size_t vi = 0; vi < numVertices; ++vi) {
            if (onReflectionPlane[d].at(vi))
                globalVertexIndex[vi] = vi;
            else {
                auto v = reflectedVertices[vi];
                globalVertexIndex[vi] = reflectedVertices.size();
                v[d] *= -1;
                reflectedVertices.push_back(v);

                // Link reflected vertex back to its original in the base cell.
                vertexOrigin.push_back(vertexOrigin.at(vi));
                vertexIsometry.push_back(vertexIsometry.at(vi).compose(refl.get()));

                // Update the onReflectionPlane info for future reflections
                for (size_t d2 = d + 1; d2 < 3; ++d2) {
                    onReflectionPlane[d2].push_back(onReflectionPlane[d2].at(vi));
                }
            }
        }
        size_t numElements = reflectedElements.size();
        for (size_t ei = 0; ei < numElements; ++ei) {
            auto re = reflectedElements[ei];
            // Reindex corner indices.
            // Note: reflection inverts the elements, so we must also permute
            // the corner indices to get positive orientation.
            // This actually matters! The inverted reflected elements cause a
            // cancellation during stiffness matrix assembly resulting in a
            // singular system.
            size_t tmp = re[0];
            re[0] = globalVertexIndex.at(re[1]);
            re[1] = globalVertexIndex.at(tmp);
            for (size_t d = 2; d < re.size(); ++d) re[d] = globalVertexIndex.at(re[d]);
            reflectedElements.push_back(re);
        }
    }
}

// Same as above, but automatically determine vertices on the reflection plane
template<typename Vertex, typename Element>
void reflectXYZ(size_t Dim, // Dimensions to reflect in (length of [x, y, z] prefix)
                const std::vector<Vertex> &vertices,
                const std::vector<Element> &elements,
                std::vector<Vertex>  &reflectedVertices,
                std::vector<Element> &reflectedElements,
                std::vector<size_t>   &vertexOrigin,
                std::vector<Isometry> &vertexIsometry)
{
    std::vector<std::vector<bool>> onReflectionPlane(3, std::vector<bool>(vertices.size(), false));
    for (size_t i = 0; i < vertices.size(); ++i) {
        onReflectionPlane[0][i] = (vertices[i][0] == 0);
        onReflectionPlane[1][i] = (vertices[i][1] == 0);
        onReflectionPlane[2][i] = (vertices[i][2] == 0);
    }
    reflectXYZ(Dim, vertices, elements, onReflectionPlane, reflectedVertices, reflectedElements,
               vertexOrigin, vertexIsometry);
}

// Same as above, but don't output vertexOrigin/vertexIsometry
template<typename Vertex, typename Element>
void reflectXYZ(size_t Dim, // Dimensions to reflect in (length of [x, y, z] prefix)
                const std::vector<Vertex> &vertices,
                const std::vector<Element> &elements,
                std::vector<Vertex>  &reflectedVertices,
                std::vector<Element> &reflectedElements)
{
    std::vector<size_t>   vertexOrigin;
    std::vector<Isometry> vertexIsometry;
    reflectXYZ(Dim, vertices, elements, reflectedVertices, reflectedElements,
               vertexOrigin, vertexIsometry);
}

#endif /* end of include guard: SNAPANDREFLECT_HH */
