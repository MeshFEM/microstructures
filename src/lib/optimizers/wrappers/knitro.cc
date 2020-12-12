////////////////////////////////////////////////////////////////////////////////
// knitro.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Interface with knitro's active set (SQP) implementation.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  02/07/2018 15:37:24
////////////////////////////////////////////////////////////////////////////////
#include "knitro.hh"

#if HAS_KNITRO

#include <KTRSolver.h>
#include <KTRProblem.h>
#include <limits>

#include <vector>
#include <cassert>

#include "../Constraint.hh"

using namespace PatternOptimization;
using std::vector;

struct MicrostructureDesignProblem : public knitro::KTRProblem {
    // constructor: pass number of variables and constraints to base class
    MicrostructureDesignProblem(IterateManagerBase &im, const size_t nParams, const size_t nConstraints, const std::string &outPath)
            : KTRProblem(nParams, nConstraints), m_im(im), m_nParams(nParams), m_nConstraints(nConstraints), m_outPath(outPath)
    {
        setObjType(KTR_OBJTYPE_GENERAL);
        setObjGoal(KTR_OBJGOAL_MINIMIZE);
    }

    double evaluateFC(const vector<double> &x,
                            vector<double>& cval,
                            vector<double>& objGrad,
                            vector<double>& jac) {
        assert(x.size() == m_nParams);
        assert(cval.size() == m_nConstraints);
        assert(objGrad.size() == m_nParams);
        assert(jac.size() == m_nConstraints * m_nParams);
        const auto &it = m_im.get(m_nParams, &x[0]);

        // Evaluate constraints and constraint jacobian
        for (size_t cnum = 0, offset = 0; cnum < it.numConstraints(); ++cnum) {
            const auto &c = it.evaluatedConstraint(cnum);
            for (size_t comp = 0; comp < c.dimension(); ++comp)
                cval.at(offset + comp) = c.values[comp];
            offset += c.dimension();
        }

        // Evaluate objective
        double val = it.evaluate();

        // Evaluate objective gradient
        it.gradp().getFlattened(objGrad);

        // Evaluate constraint jacobian
        for (size_t cnum = 0, offset = 0; cnum < it.numConstraints(); ++cnum) {
            const auto &c = it.evaluatedConstraint(cnum);
            for (size_t comp = 0; comp < c.dimension(); ++comp) {
                const size_t i = offset + comp;
                for (size_t j = 0; j < m_nParams; ++j)
                     jac.at(i * m_nParams + j) = c.jacobian(comp, j);
            }
            offset += c.dimension();
        }

        if (isImprovement(val, x) && it.shouldReport()) {
            it.writeDescription(std::cout);
            std::cout << std::endl;
            if (m_outPath != "")
                it.writeMeshAndFields(m_outPath + "_" + std::to_string(m_niters));
        }

        return val;
    }

    // Gradient is evaluated in evaluateFC
    int evaluateGA(const vector<double> &x, vector<double> &objGrad, vector<double>&jac) {
        // According to Knitro's documentation, it suffices to simply implement
        // evaluateFC, but this isn't working for some reason...
        vector<double> cval(m_nConstraints);
        evaluateFC(x, cval, objGrad, jac);
        return 0;
    }

    void setVariableBounds(const BoundConstraints &bds) {
        for (size_t p = 0; p < m_nParams; ++p) {
            if (!(bds.hasLowerBound.at(p) && bds.hasUpperBound.at(p)))
                throw std::runtime_error("Currently all parameters must be bounded.");
        }
        assert(bds.lowerBound.size() == m_nParams);
        assert(bds.upperBound.size() == m_nParams);
        setVarLoBnds(bds.lowerBound);
        setVarUpBnds(bds.upperBound);
    }

    void configureConstraints(const vector<double> &lb, const vector<double> &ub) {
        assert(lb.size() == m_nConstraints);
        assert(ub.size() == m_nConstraints);
        setConTypes(knitro::KTREnums::ConstraintType::ConGeneral);
        setConLoBnds(lb);
        setConUpBnds(ub);
    }

private:
    // TODO: When constraints are involved, the iterate can be an improvement if
    // either infeasibility decreases or objective decreases.
    // (Currently just always reports, unless the iterate is a repeat)
    bool isImprovement(Real /* val */, const vector<Real> &x) {
        bool differ = true;
        if (m_prevParams.size() == x.size()) {
            differ = false;
            for (size_t i = 0; i < x.size(); ++i) {
                if (std::abs(m_prevParams[i] - x[i]) > 1e-9) {
                    differ = true;
                    break;
                }
            }
        }
        m_prevParams = x;
        if (!differ) return false;
#if 0
        if (val < best) {
            best = val;
            ++m_niters;
            return true;
        }
        return false;
#endif
		++m_niters;
        return true;
    }

    IterateManagerBase &m_im;
    const size_t m_nParams, m_nConstraints;
    size_t m_niters = 0;
    std::string m_outPath;
    vector<double> m_prevParams;
};

void optimize_knitro_active_set(ScalarField<Real> &params,
        const BoundConstraints &bds,
        IterateManagerBase &im,
        const OptimizerConfig &oconfig,
        const std::string &outPath)
{
    const size_t nParams = params.domainSize();
    const auto &it = im.get(nParams, &params[0]);
    const size_t nVecConstraints = it.numConstraints();

    size_t nConstraints = 0;
    vector<double> constraintLB, constraintUB, feasTol;
    // Count constraint components and determine constraint properties
    // We use a lower tolerance for equality (tensor fitting) constraints than
    // inequality (printability).
    for (size_t i = 0; i < nVecConstraints; ++i) {
        const auto &c = it.evaluatedConstraint(i);
        nConstraints += c.dimension();
        for (size_t j = 0; j < c.dimension(); ++j) {
            // Inequality constraints are of the form c_i <= 0
            // Eqality    constraints are of the form c_i  = 0
            constraintUB.push_back(0.0);
            if      (c.type == ConstraintType::  EQUALITY) { constraintLB.push_back(0.0);           feasTol.push_back( 1e-3); }
            else if (c.type == ConstraintType::INEQUALITY) { constraintLB.push_back(-KTR_INFBOUND); feasTol.push_back(1e-14); }
            else throw std::runtime_error("Unexpected constraint type");
        }
    }

    MicrostructureDesignProblem problem(im, nParams, nConstraints, outPath);
    problem.setVariableBounds(bds);
    problem.configureConstraints(constraintLB, constraintUB);
    for (size_t p = 0; p < nParams; ++p)
        problem.setXInitial(p, params[p]);

    // Create a solver - optional arguments:
    // exact first derivatives
    // BFGS approximation for second derivatives
    knitro::KTRSolver solver(&problem, KTR_GRADOPT_EXACT, KTR_HESSOPT_BFGS);
    solver.setParam(KTR_PARAM_HONORBNDS, KTR_HONORBNDS_ALWAYS); // always respect bounds during optimization
    solver.setParam(KTR_PARAM_MAXIT, int(oconfig.niters));
    solver.setParam(KTR_PARAM_PRESOLVE, KTR_PRESOLVE_NONE);
    solver.setParam(KTR_PARAM_ALGORITHM, KTR_ALG_ACT_SQP);   // active set with BFGS Hessian approximation
    // solver.setParam(KTR_PARAM_ALGORITHM, KTR_ALG_IPDIRECT);
    // solver.setParam(KTR_PARAM_BAR_FEASIBLE, KTR_BAR_FEASIBLE_NO);

    solver.setFeasTols(feasTol, vector<double>(), vector<double>());

    try {
        int solveStatus = solver.solve();

        if (solveStatus != 0) {
            std::cout << std::endl;
            std::cout << "KNITRO failed to solve the problem, final status = ";
            std::cout << solveStatus << std::endl;
        }
        else {
            std::cout << std::endl << "KNITRO successful" << std::endl;
        }
    params = solver.getXValues();

    }
    catch (knitro::KTRException &e) {
        e.printMessage();
        throw e;
    }
}

#else // !HAS_KNITRO

void optimize_knitro_active_set(ScalarField<Real> &/* params */,
        const PatternOptimization::BoundConstraints &/* bds */,
        PatternOptimization::IterateManagerBase &/* im */,
        const PatternOptimization::OptimizerConfig &/* oconfig */,
        const std::string &/* outPath */) {
    throw std::runtime_error("Built without Knitro");
}

#endif // HAS_KNITRO
