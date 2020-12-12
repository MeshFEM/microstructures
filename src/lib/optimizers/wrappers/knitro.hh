#ifndef MY_KNITRO_HH
#define MY_KNITRO_HH

#include "../IterateManagerBase.hh"
#include "../BoundConstraints.hh"
#include "../OptimizerConfig.hh"

#include <MeshFEM/EdgeFields.hh>
#include <string>

void optimize_knitro_active_set(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);

#endif /* end of include guard: MY_KNITRO_HH */

