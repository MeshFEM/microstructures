################################################################################

get_filename_component(CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CURRENT_DIR})

################################################################################

micro_single_app(DiscreteShapeDerivativeValidation_cli MeshFEM micro::pattern_optimization micro::inflators micro::tbb meshfem::boost)
micro_single_app(MacroMicroStress_cli micro::pattern_optimization MeshFEM micro::tbb meshfem::boost)
micro_single_app(WCS_cli micro::pattern_optimization MeshFEM micro::tbb meshfem::boost)
micro_single_app(WCSOptimization_cli micro::pattern_optimization micro::optimizers micro::inflators json::json micro::tbb meshfem::boost)
# micro_single_app(WCSRemeshingOptimization_cli micro::pattern_optimization)
micro_single_app(GenBoundaryPerturbationInflatorJob micro::pattern_optimization micro::inflators json::json meshfem::boost)