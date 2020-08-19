#ifndef ITERATEMANAGERBASE_HH
#define ITERATEMANAGERBASE_HH

#include "IterateBase.hh"

namespace PatternOptimization {

// Enables creation and sharing of a iterates amongst several objects (e.g.
// CeresCostWrapper(s) and IterationCallback)
struct IterateManagerBase {
    virtual IterateBase &get(size_t nParams, const double * const params) = 0;

    virtual       IterateBase    &get()       = 0;
    virtual const IterateBase    &get() const = 0;
    virtual       IterateBase *getPtr()       = 0;
    virtual const IterateBase *getPtr() const = 0;

    virtual size_t numParameters() const = 0;
    virtual bool    isParametric() const = 0;

    virtual void update() = 0;

    // Functions to manage solution
    virtual void setOutPath(std::string outPath) = 0;
    virtual bool areNewParameters(const std::vector<Real> &params) = 0;
    virtual bool isFeasible(const std::vector<Real> &params) = 0;
    virtual bool isImprovement(Real val, bool respectConstraints) = 0;
    virtual void updateAndReport(const std::vector<Real> &x) = 0;
    virtual size_t numberIterations() = 0;
    virtual size_t bestIteration() = 0;

    virtual ~IterateManagerBase() { }
};

}

#endif /* end of include guard: ITERATEMANAGERBASE_HH */
