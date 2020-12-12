#ifndef GRADIENT_DESCENT_HH
#define GRADIENT_DESCENT_HH

#include "../IterateManagerBase.hh"
#include "../BoundConstraints.hh"
#include "../OptimizerConfig.hh"

#include <MeshFEM/EdgeFields.hh>
#include <string>

void optimize_gd(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);

// gradient descent with custom line search
void optimize_gd_smartstep(ScalarField<Real> &params,
                           const PatternOptimization::BoundConstraints &bds,
                           PatternOptimization::IterateManagerBase &im,
                           const PatternOptimization::OptimizerConfig &oconfig,
                           const std::string &outPath);

#endif /* end of include guard: GRADIENT_DESCENT_HH */
