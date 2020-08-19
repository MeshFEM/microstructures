################################################################################
# Source files
################################################################################

set(SOURCES
        BoundConstraints.hh
        Constraint.hh
        EvaluatedConstraint.hh
        EvaluatedObjectiveTerm.hh
        IterateBase.hh
        IterateManagerBase.hh
        OptimizerConfig.hh
        wrappers/ceres.cc
        wrappers/ceres.hh
        wrappers/dlib.cc
        wrappers/dlib.hh
        wrappers/gradient_descent.cc
        wrappers/gradient_descent.hh
        wrappers/knitro.cc
        wrappers/knitro.hh
        wrappers/nlopt.cc
        wrappers/nlopt.hh
)

micro_add_library(optimizers ${SOURCES})
micro_copy_headers(optimizers ${SOURCES})

################################################################################
# Dependencies
################################################################################

target_link_libraries(micro_optimizers
    PUBLIC
        MeshFEM
        micro::inflators
        micro::accelerate
    PRIVATE
        micro::dlib
        micro::ceres
        micro::nlopt
        micro::knitro
        warnings::all
)

if(MICRO_WITH_SANITIZERS)
    add_sanitizers(micro_optimizers)
endif()
