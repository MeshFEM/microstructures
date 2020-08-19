################################################################################
# Source files
################################################################################

set(SOURCES
		BaseCellOperations.hh
		IterateFactory.hh
		IterateManager.hh
		LookupTable.hh
		ObjectiveTerm.hh
		ObjectiveTermNormalizations.hh
		PatternOptimizationJob.cc
		PatternOptimizationJob.hh
		PatternOptimizationIterate.hh
		PatternOptimizationConfig.hh
		SDConversions.hh
		SDConversionsOrthoCell.hh
		ShapeVelocityInterpolatorOrthoCell.hh
		objective_terms/IsotropicFit.hh
		objective_terms/IsotropicFitRel.hh
		objective_terms/ProximityRegularization.hh
		objective_terms/TargetVolume.hh
		objective_terms/TensorFit.hh
		objective_terms/SmoothingRegularization.hh
		objective_terms/PeriodicSmoothingRegularization.hh
        objective_terms/ScaleInvariantSmoothingRegularization.hh
		constraints/Printability.hh
		constraints/TensorFit.hh
)

micro_add_library(pattern_optimization ${SOURCES})
micro_copy_headers(pattern_optimization ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_pattern_optimization
	PUBLIC
		micro::inflators
		micro::optimizers
	PRIVATE
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
	add_sanitizers(micro_pattern_optimization)
endif()
