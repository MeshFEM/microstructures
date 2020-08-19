################################################################################
# Source files
################################################################################

set(SOURCES
        MakeInflator.cc
        Inflator.hh
        MakeInflator.hh
        wrappers/BoundaryPerturbationInflator.cc
        wrappers/BoundaryPerturbationInflator.hh
        wrappers/EqualityConstrainedInflator.cc
        wrappers/EqualityConstrainedInflator.hh
        wrappers/IsoinflatorWrapper.cc
        wrappers/IsoinflatorWrapper.hh
        wrappers/JamesInflatorWrapper.cc
        wrappers/JamesInflatorWrapper.hh
        wrappers/LuigiInflatorWrapper.cc
        wrappers/LuigiInflatorWrapper.hh
)

micro_add_library(inflators ${SOURCES})
micro_copy_headers(inflators ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_inflators
    PUBLIC
        MeshFEM
        meshfem::boost
        micro::isosurface_inflator
        micro::pymesh
        micro::vcglib
        micro::clipper
    PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(micro_inflators)
endif()
