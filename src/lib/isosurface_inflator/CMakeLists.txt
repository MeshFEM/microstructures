################################################################################
# Source files
################################################################################

set(SOURCES
        AABBTree.cc
        AABBTree.hh
        AutomaticDifferentiation.hh
        BilinearMap.cc
        BilinearMap.hh
        BoxIntersection1DFeatures.hh
        BoxIntersectionMesher.cc
        BoxIntersectionMesher.hh
        CGALClippedVolumeMesher.hh
        ConvexHullTriangulation.cc
        ConvexHullTriangulation.hh
        DynamicKdTree.cc
        DynamicKdTree.hh
        FuzzySign.hh
        IGLSurfaceMesherMC.cc
        IGLSurfaceMesherMC.hh
        InflatorTypes.hh
        Isometries.hh
        IsosurfaceInflator.cc
        IsosurfaceInflator.hh
        IsosurfaceInflatorConfig.hh
        IsosurfaceInflatorImpl.hh
        Joint.hh
        MesherBase.hh
        MeshingOptions.cc
        MeshingOptions.hh
        MidplaneMesher.cc
        MidplaneMesher.hh
        NodePositioners.hh
        PaperVisualizationSDFunc.hh
        PatternSignedDistance.hh
        PostProcess.cc
        PostProcess.hh
        rasterize.cc
        rasterize.hh
        ShapeVelocityInterpolator.hh
        SignedDistance.hh
        SignedDistanceRegion.hh
        SnapAndReflect.hh
        SphereConvexHull.hh
        SpherePoints.hh
        StitchedWireMesh.hh
        Symmetry.hh
        TesselateSpheres.hh
        TriangleClosestPoint.hh
        TriplyPeriodicMinimalShell.hh
        VCGSurfaceMesher.cc
        VCGSurfaceMesher.hh
        WireQuadMesh.cc
        WireQuadMesh.hh
        WireMesh.hh
        WireMesh.inl
        quadfoam/jacobians.cpp
        quadfoam/jacobians.h
        quadfoam/navigation.cpp
        quadfoam/navigation.h
        quadfoam/instantiate.cpp
        quadfoam/instantiate.h
)

micro_add_library(isosurface_inflator ${SOURCES})
micro_copy_headers(isosurface_inflator ${SOURCES})

if(MICRO_WITH_CGAL)
    micro_add_library(isosurface_inflator_cgal
            CGALClippedVolumeMesher.cc
    )
endif()

################################################################################
# Dependencies
################################################################################

if(MICRO_WITH_CGAL)
    target_link_libraries(micro_isosurface_inflator_cgal
        PUBLIC
            MeshFEM
            micro::tbb
        PRIVATE
            CGAL::CGAL
    )

    target_compile_definitions(micro_isosurface_inflator PUBLIC MICRO_WITH_CGAL)
endif()

target_link_libraries(micro_isosurface_inflator
    PUBLIC
        MeshFEM
        micro::libigl
        micro::tbb
        nanoflann::nanoflann
    PRIVATE
        $<TARGET_NAME_IF_EXISTS:micro::isosurface_inflator_cgal>
        micro::vcglib
        quickhull::quickhull
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(micro_isosurface_inflator_cgal)
    add_sanitizers(micro_isosurface_inflator)
endif()
