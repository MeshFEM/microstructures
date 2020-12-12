#ifndef DLIB_HH
#define DLIB_HH

#include "../IterateManagerBase.hh"
#include "../BoundConstraints.hh"
#include "../OptimizerConfig.hh"

#include <MeshFEM/EdgeFields.hh>
#include <string>

// bfgs when oconfig.lbfgs_memory == 0, lbfgs otherwise.
void optimize_dlib_bfgs(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);

// bfgs with a custom line search
void optimize_dlib_custom_bfgs(ScalarField<Real> &params,
                        const PatternOptimization::BoundConstraints &bds,
                        PatternOptimization::IterateManagerBase &im,
                        const PatternOptimization::OptimizerConfig &oconfig,
                        const std::string &outPath);

#endif /* end of include guard: DLIB_HH */
