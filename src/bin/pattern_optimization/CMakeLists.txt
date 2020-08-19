################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(PatternOptimization_cli micro::pattern_optimization micro::optimizers micro::inflators meshfem::boost)
micro_single_app(GenIsosurfaceJob micro::pattern_optimization meshfem::boost)
set_target_properties(GenIsosurfaceJob PROPERTIES CXX_STANDARD 14)
micro_single_app(OffsetBoundsToTranslationBounds micro::isosurface_inflator micro::pattern_optimization)
