#ifndef MY_NLOPT_H
#define MY_NLOPT_H

#include "../IterateManagerBase.hh"
#include "../BoundConstraints.hh"
#include "../OptimizerConfig.hh"

#include <MeshFEM/EdgeFields.hh>
#include <string>

void optimize_nlopt_slsqp(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);

void optimize_nlopt_lbfgs(ScalarField<Real> &params,
                          const PatternOptimization::BoundConstraints &bds,
                          PatternOptimization::IterateManagerBase &im,
                          const PatternOptimization::OptimizerConfig &oconfig,
                          const std::string &outPath);

#endif /* end of include guard: MY_NLOPT_H */
