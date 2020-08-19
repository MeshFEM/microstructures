#include "PostProcess.hh"
#include "IsosurfaceInflatorImpl.hh"
#include "SnapAndReflect.hh"
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/algorithms/get_element_components.hh>
#include <MeshFEM/filters/remove_dangling_vertices.hh>
#include <MeshFEM/filters/remove_small_components.hh>
#include <MeshFEM/PeriodicBoundaryMatcher.hh>

using namespace std;

// Postprocess:
//    Snap to base cell and then reflect if necessary
//    Compute vertex normals and normal shape velocities
template<size_t N, class Point>
void postProcess(vector<MeshIO::IOVertex>  &vertices,
                 vector<MeshIO::IOElement> &elements,
                 vector<vector<Real>>      &normalShapeVelocities,
                 vector<Point>             &vertexNormals,
                 const IsosurfaceInflator::Impl &inflator,
                 bool                      meshedFullPeriodCell, // whether (vertices, elements) is already a mesh of the full period cell
                 bool                      requestFullPeriodCell,
                 const BBox<Point>         &meshCell,
                 const MeshingOptions      &opts,
                 bool cheapPostProcessing,
                 bool nonPeriodicity) {
    const bool needsReflecting = requestFullPeriodCell && !meshedFullPeriodCell;

    BENCHMARK_START_TIMER_SECTION("postProcess");
    {
        const size_t origSize = vertices.size();
        remove_dangling_vertices(vertices, elements);
        if (vertices.size() != origSize)
            std::cerr << "WARNING: removed " << origSize - vertices.size() << " dangling vertices." << std::endl;
    }

    BENCHMARK_START_TIMER("SnapAndDetermineEvaluationPts");
    vector<vector<bool>> onMinFace, onMaxFace;

    using SMesh = SimplicialMesh<N>;
    std::unique_ptr<SMesh> bcm;
    bool done = false;

    // Remove small components
    while (!done) {
        try {
            bcm = Future::make_unique<SMesh>(elements, vertices.size());
            done = true;
        }
        catch (...) {
            std::cerr << "Exception while building mesh" << std::endl;
            std::cerr << "Dumping debug.msh" << std::endl;
            MeshIO::save("debug.msh", vertices, elements);
            throw;
        }

        // Remove spurious small components. If this is already a full periodic
        // mesh, we need to detect components on the periodically-stitched mesh to
        // avoid discarding parts connected via the adjacent period cells.
        std::vector<size_t> componentIndex;
        std::vector<size_t> componentSize;
        if (!nonPeriodicity && meshedFullPeriodCell) {
            std::vector<MeshIO::IOElement> stitched_elements = elements;
            using Pt = PointND<N>;
            BBox<Point> meshBB(vertices);
            BBox<Pt> cellND(truncateFrom3D<Pt>(meshBB.minCorner), truncateFrom3D<Pt>(meshBB.maxCorner));

            auto pts = MeshFEM::apply(vertices, [](const MeshIO::IOVertex &v) { return truncateFrom3D<Pt>(v.point); });
            std::vector<PeriodicBoundaryMatcher::FaceMembership<N>> fm;
            PeriodicBoundaryMatcher::determineCellBoundaryFaceMembership<N>(pts, cellND, fm, 0 /* epsilon */);
            std::vector<std::vector<size_t>> nodeSets;
            std::vector<size_t> nodeSetForNode;
            // Could be sped up by only matching boundary vertices
            try {
                PeriodicBoundaryMatcher::match(pts, cellND, fm, nodeSets, nodeSetForNode);
            }
            catch (...) {
                std::cerr << "Exception while post processing mesh" << std::endl;
                std::cerr << "Dumping debug.msh" << std::endl;
                MeshIO::save("debug.msh", vertices, elements);
                throw;
            }
            // Relink periodic vertices to the first element of their corresponding node set.
            for (auto &e : stitched_elements) {
                for (size_t &v : e) {
                    size_t ns = nodeSetForNode.at(v);
                    if (ns != PeriodicBoundaryMatcher::NONE) v = nodeSets.at(ns).front();
                }
            }
            auto dummy_vertices = vertices;
            remove_dangling_vertices(dummy_vertices, stitched_elements);
            get_element_components(stitched_elements, componentIndex, componentSize);
        } else {
            get_element_components(elements, componentIndex, componentSize);
        }

        if (remove_small_components(componentIndex, componentSize, vertices, elements)) {
            std::cout << "Removed " << bcm->numVertices() - vertices.size() << " vertices"
                      << " and " << bcm->numSimplices() - elements.size() << " elements"
                      << " (small components)" << std::endl;
            done = false; // Rebuild mesh and try again.
        }
    }

    SMesh &symBaseCellMesh = *bcm;

    if (nonPeriodicity) {
        // We should not compute boundary velocities for elements touching the meshCell boundary.
        using FM = PeriodicBoundaryMatcher::FaceMembership<3>;
        onMinFace.assign(3, std::vector<bool>(vertices.size(), false));
        onMaxFace.assign(3, std::vector<bool>(vertices.size(), false));
        for (size_t i = 0; i < vertices.size(); ++i) {
            FM fm(vertices[i], meshCell, 0); // zero tolerance!
            onMinFace[0][i] = fm.onMinFace(0), onMaxFace[0][i] = fm.onMaxFace(0);
            onMinFace[1][i] = fm.onMinFace(1), onMaxFace[1][i] = fm.onMaxFace(1);
            onMinFace[2][i] = fm.onMinFace(2), onMaxFace[2][i] = fm.onMaxFace(2);
        }
    }
    else {
        try {
            if (N == 3) smartSnap3D(vertices, symBaseCellMesh, meshCell);
        }
        catch (std::exception &e) {
            std::cerr
                    << "WARNING: smartSnap3D failed--probably nonsmooth geometry at cell interface. Resorting to dumbSnap3D."
                    << std::endl;
            std::cerr << "(" << e.what() << ")" << std::endl;
            dumbSnap3D(vertices, meshCell, opts.facetDistance);
        }
        BBox<Point3D> snappedBB(vertices);
        if (N == 2) {
            // We don't care about the z-depth of the bounding box
            snappedBB.minCorner[2] = meshCell.minCorner[2];
            snappedBB.maxCorner[2] = meshCell.maxCorner[2];
        }
        if (snappedBB != meshCell) {
            std::cerr << "snappedBB: " << snappedBB << std::endl;
            std::cerr << "meshCell: " << meshCell << std::endl;
            std::cerr << "Failed to snap mesh. Dumping debug.msh" << std::endl;
            MeshIO::save("debug.msh", vertices, elements);
            throw std::runtime_error("Snapping failed.");
        }

        using FM = PeriodicBoundaryMatcher::FaceMembership<3>;
        onMinFace.assign(3, std::vector<bool>(vertices.size(), false));
        onMaxFace.assign(3, std::vector<bool>(vertices.size(), false));
        std::vector<FM> faceMembership;
        for (size_t i = 0; i < vertices.size(); ++i) {
            FM fm(vertices[i], meshCell, 0); // zero tolerance!
            onMinFace[0][i] = fm.onMinFace(0), onMaxFace[0][i] = fm.onMaxFace(0);
            onMinFace[1][i] = fm.onMinFace(1), onMaxFace[1][i] = fm.onMaxFace(1);
            onMinFace[2][i] = fm.onMinFace(2), onMaxFace[2][i] = fm.onMaxFace(2);
        }
    }
    // MeshIO::save("post_snap.msh", vertices, elements);


    // Mark internal cell-face vertices: vertices on the meshing cell
    // boundary that actually lie inside the object (i.e. they are only mesh
    // boundary vertices because of the intersection of the periodic pattern with
    // the meshing box).
    // This this is not the case if any non-cell-face triangle is incident
    vector<bool> internalCellFaceVertex(symBaseCellMesh.numBoundaryVertices(), true);
    for (auto bs : symBaseCellMesh.boundarySimplices()) {
        bool isCellFace = false;
        for (size_t d = 0; d < N; ++d) {
            bool onMin = true, onMax = true;
            for (auto bv : bs.vertices()) {
                onMin &= onMinFace[d].at(bv.volumeVertex().index());
                onMax &= onMaxFace[d].at(bv.volumeVertex().index());
            }
            isCellFace |= (onMin | onMax);
        }
        if (isCellFace) continue;

        for (auto bv : bs.vertices())
            internalCellFaceVertex.at(bv.index()) = false;
    }

    // Compute parameter shape velocities on all (true) boundary vertices
    vector<Point> evaluationPoints;
    // Tie the boundary vertices to the associated data evaluation point
    vector<size_t> evalPointIndex(symBaseCellMesh.numBoundaryVertices(),
                                  numeric_limits<size_t>::max());
    for (auto bv : symBaseCellMesh.boundaryVertices()) {
        if (internalCellFaceVertex.at(bv.index())) continue;
        evalPointIndex[bv.index()] = evaluationPoints.size();
        evaluationPoints.push_back(vertices.at(bv.volumeVertex().index()));
    }

#if DEBUG_EVALPTS
    {
        MSHFieldWriter debug("debug_evalpts.msh", vertices, elements);
        for (size_t d = 0; d < N; ++d) {
            ScalarField<Real> minFaceIndicator(vertices.size()), maxFaceIndicator(vertices.size());
            for (size_t i = 0; i < vertices.size(); ++i) {
                minFaceIndicator[i] = onMinFace[d].at(i) ? 1.0 : 0.0;
                maxFaceIndicator[i] = onMaxFace[d].at(i) ? 1.0 : 0.0;
            }
            debug.addField("onMinFace[" + std::to_string(d) + "]", minFaceIndicator);
            debug.addField("onMaxFace[" + std::to_string(d) + "]", maxFaceIndicator);
        }
        ScalarField<Real> evalPtIndicator(vertices.size());
        ScalarField<Real> icfv(vertices.size());
        evalPtIndicator.clear();
        icfv.clear();
        for (auto bv : symBaseCellMesh.boundaryVertices()) {
            if (evalPointIndex[bv.index()] < evaluationPoints.size())
                evalPtIndicator[bv.volumeVertex().index()] = 1.0;
            icfv[bv.volumeVertex().index()] = internalCellFaceVertex.at(bv.index()) ? 1.0 : 0.0;
        }
        debug.addField("isEvalPoint", evalPtIndicator);
        debug.addField("internalCellFaceVertex", icfv);
    }
#endif

    BENCHMARK_STOP_TIMER("SnapAndDetermineEvaluationPts");

    vector<vector<Real>> vnp;
    vector<Point> sdGradX;
    vector<Real> sdGradNorms(evaluationPoints.size());

    if (!cheapPostProcessing) {
        BENCHMARK_START_TIMER("SignedDistanceGradientsAndPartials");
        // sd(x, p) = 0
        // grad_x(sd) . dx/dp + d(sd)/dp = 0,   grad_x(sd) = n |grad_x(sd)|
        // ==>  v . n = -[ d(sd)/dp ] / |grad_x(sd)|

        // try {
        vnp = inflator.signedDistanceParamPartials(evaluationPoints);
        sdGradX = inflator.signedDistanceGradient(evaluationPoints);
        // }
        // catch(...) {
        //     MSHFieldWriter debug("debug.msh", vertices, elements);
        //     BENCHMARK_STOP_TIMER_SECTION("SignedDistanceGradientsAndPartials");
        //     BENCHMARK_STOP_TIMER_SECTION("postProcess");
        //     throw;
        // }

#if 0
        {
            // Debug the smoothing modulation field (is it causing bad signed
            // distance partials?)
            vector<Real> smoothness;
            vector<size_t> closestVtx;
            std::tie(smoothness, closestVtx) = inflator.smoothnessAtPoints(evaluationPoints);
            assert(smoothness.size() == evaluationPoints.size());
            {
                Real avg = 0.0;
                for (Real s : smoothness) avg += s;
                avg /= smoothness.size();
                std::cout << "Average smoothness: " << avg << std::endl;
            }

            ScalarField<Real> smoothnessField(symBaseCellMesh.numBoundaryVertices()),
                              closestVtxField(symBaseCellMesh.numBoundaryVertices());
            smoothnessField.clear(), closestVtxField.clear();
            for (auto bv : symBaseCellMesh.boundaryVertices()) {
                size_t idx = evalPointIndex.at(bv.index());
                if (idx > smoothness.size()) continue;
                smoothnessField[bv.index()] = smoothness[idx];
                closestVtxField[bv.index()] = closestVtx[idx];
            }

            // Extract boundary vertices/elements for output
            std::vector<MeshIO::IOVertex > outVertices;
            std::vector<MeshIO::IOElement> outElements;
            for (auto bv : symBaseCellMesh.boundaryVertices())
                outVertices.push_back(vertices.at(bv.volumeVertex().index()));
            for (auto bs : symBaseCellMesh.boundarySimplices()) {
                outElements.emplace_back(bs.vertex(0).index(),
                                         bs.vertex(1).index(),
                                         bs.vertex(2).index());
            }
            MSHFieldWriter writer("smoothness_debug.msh", outVertices, outElements);
            writer.addField("smoothness", smoothnessField, DomainType::PER_NODE);
            writer.addField("closest_vtx", closestVtxField, DomainType::PER_NODE);
        }
#endif
        for (size_t i = 0; i < evaluationPoints.size(); ++i) {
            sdGradNorms[i] = sdGradX[i].norm();
            // We evaluate on the boundary--there should be a well-defined normal
            if (std::isnan(sdGradNorms[i]) || (std::abs(sdGradNorms[i]) < 1e-8)) {
                MSHFieldWriter debug("debug.msh", vertices, elements);
                BENCHMARK_STOP_TIMER("SignedDistanceGradientsAndPartials");
                BENCHMARK_STOP_TIMER_SECTION("postProcess");

                std::cerr << "Undefined normal at evaluation point "
                          << evaluationPoints[i] << std::endl;

                BENCHMARK_STOP_TIMER("SignedDistanceGradientsAndPartials");
                BENCHMARK_STOP_TIMER_SECTION("postProcess");
                throw std::runtime_error("Normal undefined.");
            }
        }

        for (auto &vn : vnp) {
            for (size_t i = 0; i < vn.size(); ++i) {
                vn[i] *= -1.0 / sdGradNorms[i];
                if (std::isnan(vn[i])) {
                    ScalarField<Real> fail(vertices.size());
                    fail.clear();
                    for (const auto bv : symBaseCellMesh.boundaryVertices()) {
                        size_t e = evalPointIndex.at(bv.index());
                        if (e < evaluationPoints.size()) {
                            if (std::isnan(vn.at(e))) fail[bv.volumeVertex().index()] = 1.0;
                        }
                    }

                    MSHFieldWriter debug("debug.msh", vertices, elements);
                    debug.addField("fail", fail);
                    BENCHMARK_STOP_TIMER("SignedDistanceGradientsAndPartials");
                    BENCHMARK_STOP_TIMER_SECTION("postProcess");
                    throw std::runtime_error("nan vn");
                    // assert(false);
                }
            }
        }
        BENCHMARK_STOP_TIMER("SignedDistanceGradientsAndPartials");

        BENCHMARK_START_TIMER("Reflecting");
        // If the mesher only generates the orthotropic base cell, the mesh must
        // be reflected to get the full period cell (if requested).
        if (needsReflecting) {
            vector<MeshIO::IOVertex> reflectedVertices;
            vector<MeshIO::IOElement> reflectedElements;
            vector<size_t> vertexOrigin;
            vector<Isometry> vertexIsometry;
            reflectXYZ(N, vertices, elements, onMinFace,
                       reflectedVertices, reflectedElements,
                       vertexOrigin, vertexIsometry);

            // Copy over normal velocity and vertex normal data.
            // So that we don't rely on consistent boundary vertex numbering
            // (which is nevertheless guaranteed by our mesh datastructure),
            // we store this data per-vertex (setting the data on internal
            // vertices to 0).
            normalShapeVelocities.assign(vnp.size(), vector<Real>(reflectedVertices.size()));
            vertexNormals.assign(reflectedVertices.size(), Point::Zero());

            for (size_t i = 0; i < reflectedVertices.size(); ++i) {
                auto v = symBaseCellMesh.vertex(vertexOrigin[i]);
                assert(v);
                auto bv = v.boundaryVertex();
                if (!bv || internalCellFaceVertex[bv.index()]) continue;
                size_t evalIdx = evalPointIndex.at(bv.index());
                assert(evalIdx < evaluationPoints.size());
                // Transform normals by the reflection isometry and normalize
                vertexNormals[i] = vertexIsometry[i].apply(sdGradX[evalIdx]) / sdGradNorms[evalIdx];
                for (size_t p = 0; p < vnp.size(); ++p)
                    normalShapeVelocities[p][i] = vnp[p][evalIdx];
            }

            reflectedVertices.swap(vertices);
            reflectedElements.swap(elements);
        } else {
            // Otherwise, we still need to copy over the normal shape velocities
            normalShapeVelocities.assign(vnp.size(), vector<Real>(vertices.size()));
            vertexNormals.assign(vertices.size(), Point::Zero());
            for (auto bv : symBaseCellMesh.boundaryVertices()) {
                size_t evalIdx = evalPointIndex[bv.index()];
                if (evalIdx >= evalPointIndex.size()) continue;
                size_t vi = bv.volumeVertex().index();
                vertexNormals[vi] = sdGradX[evalIdx] / sdGradNorms[evalIdx];
                for (size_t p = 0; p < vnp.size(); ++p)
                    normalShapeVelocities[p][vi] = vnp[p][evalIdx];
            }
        }
        BENCHMARK_STOP_TIMER("Reflecting");
    }

    // static int _run_num = 0;
    // MeshIO::save("debug_" + std::to_string(_run_num++) + ".msh", vertices, elements);
    BENCHMARK_STOP_TIMER_SECTION("postProcess");
}
////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations
////////////////////////////////////////////////////////////////////////////////
template void postProcess<2, Eigen::Matrix<Real, 3, 1>>(
                 std::vector<MeshIO::IOVertex>          &vertices,
                 std::vector<MeshIO::IOElement>         &elements,
                 std::vector<std::vector<Real>>         &normalShapeVelocities,
                 std::vector<Eigen::Matrix<Real, 3, 1>> &vertexNormals,
                 const IsosurfaceInflator::Impl         &inflator,
                 bool                                    meshedFullPeriodCell,
                 bool                                    requestFullPeriodCell,
                 const BBox<Eigen::Matrix<Real, 3, 1>>  &meshCell,
                 const MeshingOptions                   &opts,
                 bool cheapPostProcessing,
                 bool nonPeriodicity);

template void postProcess<3, Eigen::Matrix<Real, 3, 1>>(
                 std::vector<MeshIO::IOVertex>          &vertices,
                 std::vector<MeshIO::IOElement>         &elements,
                 std::vector<std::vector<Real>>         &normalShapeVelocities,
                 std::vector<Eigen::Matrix<Real, 3, 1>> &vertexNormals,
                 const IsosurfaceInflator::Impl         &inflator,
                 bool                                    meshedFullPeriodCell,
                 bool                                    requestFullPeriodCell,
                 const BBox<Eigen::Matrix<Real, 3, 1>>  &meshCell,
                 const MeshingOptions                   &opts,
                 bool cheapPostProcessing,
                 bool nonPeriodicity);
