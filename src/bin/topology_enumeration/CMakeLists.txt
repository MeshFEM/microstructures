################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

if(MICRO_WITH_CGAL)
    micro_add_library(intersection_check intersection_check.cc)
    target_link_libraries(micro_intersection_check micro::isosurface_inflator MeshFEM CGAL::CGAL warnings::all)
endif()

################################################################################
# Independent binaries
################################################################################

micro_single_app(match MeshFEM)
micro_single_app(sort MeshFEM)
if(MICRO_WITH_CGAL)
    micro_single_app(enumerate MeshFEM micro::intersection_check)
    micro_single_app(enumerate_2D MeshFEM micro::intersection_check)
    micro_single_app(pattern_info MeshFEM micro::isosurface_inflator)
endif()
