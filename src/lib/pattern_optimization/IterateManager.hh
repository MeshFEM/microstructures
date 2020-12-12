////////////////////////////////////////////////////////////////////////////////
// IterateManager.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Wraps the creation and re-use of iterates for particular parameter
//      values. (Repeated calls to get() with the same parameters only create
//      one iterate).
//
//      IterateManagerBase provides a non-templated interface to the iterate
//      manager.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  04/28/2016 17:30:15
////////////////////////////////////////////////////////////////////////////////
#ifndef ITERATEMANAGER_HH
#define ITERATEMANAGER_HH

#include "PatternOptimizationIterate.hh"
#include <optimizers/IterateManagerBase.hh>
#include <memory>

namespace PatternOptimization {

// Implementation of IterateManagerBase for a particular IterateFactory.
template<class _ItFactory>
struct IterateManager : public IterateManagerBase {
    using Iterate  = typename _ItFactory::Iterate;

    IterateManager(std::unique_ptr<_ItFactory> itFactory)
        : m_iterFactory(std::move(itFactory)), m_best(std::numeric_limits<Real>::max()) { }

    Iterate &get(size_t nParams, const double * const params) override {
        m_currIterate = m_iterFactory->getIterate(std::move(m_currIterate), nParams, &params[0]);
        return *m_currIterate;
    }

    virtual       Iterate &get()          override { assert(m_currIterate); return *m_currIterate; }
    virtual const Iterate &get()    const override { assert(m_currIterate); return *m_currIterate; }

    virtual       Iterate *getPtr()       override { return m_currIterate.get(); }
    virtual const Iterate *getPtr() const override { return m_currIterate.get(); }

    virtual size_t numParameters()  const override { return m_iterFactory->numParameters(); }
    virtual bool    isParametric()  const override { return m_iterFactory->isParametric(); }

    // Use information of current solution to update inflator
    virtual void update() override {
        m_iterFactory->update();
    }

    virtual void setOutPath(std::string outPath) override {
        m_outPath = outPath;
    }

    virtual bool areNewParameters(const std::vector<Real> &params) override {
        bool differ = true;
        if (m_prevParams.size() == params.size()) {
            differ = false;
            for (size_t i = 0; i < params.size(); ++i) {
                if (std::abs(m_prevParams[i] - params[i]) > 1e-9) {
                    differ = true;
                    break;
                }
            }
        }
        m_prevParams = params;
        if (!differ) return false;

        ++m_niters;
        return true;
    }

    virtual bool isFeasible(const std::vector<Real> &params) override {
        auto &it = this->get(params.size(), params.data());

        return it.isFeasible(EQ_TOL, INEQ_TOL);
    }

    virtual bool isImprovement(Real val, bool respectConstraints) override {
        if (val < m_best && respectConstraints) {
            m_best = val;
            m_bestIter = m_niters;
            return true;
        }
        return false;
    }

    virtual void updateAndReport(const std::vector<Real> &x) override {
        auto &it = this->get();
        double val = it.evaluate();

        if (areNewParameters(x) && it.shouldReport()) {
            std::cout << "Iteration: " << m_niters << std::endl;
            it.writeDescription(std::cout);

            bool successful = isImprovement(val, it.isFeasible(EQ_TOL, INEQ_TOL));
            if (successful) {
                std::cout << "\\o/ GREAT! Found a good solution that is also the best so far..." << std::endl;

                // Uses the opportunity to update the inflators, if needed
                update();
            }
            else {
                std::cout << "=( Not the best solution so far..." << std::endl;
            }

            std::cout << "Write mesh" << std::endl;
            if (m_outPath != "") {
                it.writeMeshAndFields(m_outPath + "_" + std::to_string(m_niters));

                // If succeeded, then also print two extra files with information about current best solution
                if (successful) {
                    std::string infoPath = m_outPath + "_sol.txt";
                    std::ofstream out(infoPath);
                    out << "Iteration: " << m_niters << std::endl;
                    it.writeDescription(out);
                    it.writeMeshAndFields(m_outPath + "_best");
                    out.close();
                }
            }

            std::cout << std::endl;
        }
    }

    virtual size_t numberIterations() override {
        return m_niters;
    }

    virtual size_t bestIteration() override {
        return m_bestIter;
    }

    virtual ~IterateManager() { }
private:
    std::unique_ptr<_ItFactory> m_iterFactory;
    std::unique_ptr<Iterate>    m_currIterate;

    std::vector<Real> m_prevParams;
    size_t m_niters = 0;
    size_t m_bestIter = 0;
    Real m_best;
    std::string m_outPath;

    static constexpr Real EQ_TOL = 1e-2;
    static constexpr Real INEQ_TOL = 1e-14;
};

template<class _IF>
std::shared_ptr<IterateManager<_IF>> make_iterate_manager(std::unique_ptr<_IF> itFactory) {
    return std::make_shared<IterateManager<_IF>>(std::move(itFactory));
}

}

#endif /* end of include guard: ITERATEMANAGER_HH */
