################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(isosurface_cli micro::isosurface_inflator meshfem::boost)
micro_single_app(replicate micro::isosurface_inflator meshfem::boost)
micro_single_app(tile micro::isosurface_inflator meshfem::boost)
micro_single_app(stitch_cells_cli micro::isosurface_inflator CLI11::CLI11 json::json)
set_target_properties(stitch_cells_cli PROPERTIES CXX_STANDARD 14)
micro_single_app(VisualizeQuantities micro::isosurface_inflator)
micro_single_app(VisualizeQuantities2D micro::isosurface_inflator)
micro_single_app(PrintabilityCheck micro::isosurface_inflator micro::pattern_optimization)
micro_single_app(PrintabilityConstraints micro::isosurface_inflator meshfem::boost)
micro_single_app(SignedDistance_cli micro::isosurface_inflator meshfem::boost)
micro_single_app(TestWireMesh micro::isosurface_inflator)
micro_single_app(TestStitchedWireMesh micro::isosurface_inflator)
micro_single_app(TriplyPeriodicMinimalShell micro::isosurface_inflator)
if(TARGET igl::core)
	micro_single_app(lookup_and_stitch micro::isosurface_inflator igl::core meshfem::boost)
endif()
