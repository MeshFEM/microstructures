////////////////////////////////////////////////////////////////////////////////
// IterateBase.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Non-templated, dynamic interface to the pattern optimization iterate.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  05/02/2016 19:08:21
////////////////////////////////////////////////////////////////////////////////
#ifndef ITERATEBASE_HH
#define ITERATEBASE_HH

#include "EvaluatedObjectiveTerm.hh"
#include "EvaluatedConstraint.hh"
#include <MeshFEM/EdgeFields.hh>
#include <vector>
#include <memory>
#include <stdexcept>

namespace PatternOptimization {

struct IterateBase {
    using EOTPtr = std::unique_ptr<EvaluatedObjectiveTerm>;
    using  ECPtr = std::unique_ptr<EvaluatedConstraint>;

    using SField = ScalarField<Real>;
    // parametricOptimization: whether we're running a lower-dimensional
    // parametric optimization where gradients must be computed by taking inner
    // products with the shape velocity fields.
    //
    // For free boundary shape optimization problems, this inner product
    // approach to gradients is intractable and unnecessary.
    IterateBase(bool parametricOptimization) : m_parametricOptimization(parametricOptimization) { }

    const EvaluatedObjectiveTerm &evaluatedObjectiveTerm(const std::string &name) const {
        for (const auto &term : m_evaluatedObjectiveTerms)
            if (term->name == name) return *term;
        throw std::runtime_error("Objective term not found: " + name);
    }

    const EvaluatedConstraint &evaluatedConstraint(const std::string &name) const {
        for (const auto &c : m_evaluatedConstraints)
            if (c->name == name) return *c;
        throw std::runtime_error("Constraint not found: " + name);
    }

    const EvaluatedObjectiveTerm &evaluatedObjectiveTerm(size_t i) const {
        return *m_evaluatedObjectiveTerms.at(i);
    }

    const EvaluatedConstraint &evaluatedConstraint(size_t i) const {
        return *m_evaluatedConstraints.at(i);
    }

    const std::vector<EOTPtr> &evaluatedObjectiveTerms() const {
        return m_evaluatedObjectiveTerms;
    }

    const std::vector<ECPtr> &evaluatedConstraints() const {
        return m_evaluatedConstraints;
    }

    virtual size_t numObjectiveTerms() const { return m_evaluatedObjectiveTerms.size(); }
    virtual size_t numConstraints()    const { return m_evaluatedConstraints.size(); }

    ////////////////////////////////////////////////////////////////////////////
    // Full objective parametric gradient and steepest descent: linear
    // combination of evaluated objective terms' gradients.
    ////////////////////////////////////////////////////////////////////////////
    SField gradp() const {
        SField full;
        for (const auto &term : m_evaluatedObjectiveTerms) {
            if (full.domainSize() == 0) full  = term->gradContribution();
            else                        full += term->gradContribution();
        }
        return full;
    }

    // Steepest descent direction with respect to the object's L^2 boundary
    // metric.
    // Works for both parametric and non-parametric optimization.
    SField steepestDescent() const {
        SField full;
        for (const auto &term : m_evaluatedObjectiveTerms) {
            if (full.domainSize() == 0) full  = term->descentContribution();
            else                        full += term->descentContribution();
        }
        return full;
    }

    bool paramsDiffer(size_t nParams, const Real *params) const {
        assert(nParams == m_params.size());
        for (size_t i = 0; i < nParams; ++i)
            if (m_params[i] != params[i])
                return true;
        return false;
    }

    const std::vector<Real> &params() const { return m_params; }
    ////////////////////////////////////////////////////////////////////////////
    // Implemented by subclass
    ////////////////////////////////////////////////////////////////////////////
    virtual Real evaluate()                                  const = 0;
    virtual Real evaluateNormalized(const std::string &name) const = 0;
    virtual void writeMeshAndFields(const std::string &path) const = 0;
    virtual void writeDescription(std::ostream &os)          const = 0;
    virtual std::vector<std::string> objectiveTermNames()    const = 0;

    virtual bool hasConstraint(const std::string &name) const = 0;

    // Verifies if current solution is feasible: this means deciding if it
    // respects or not the imposed constraints
    virtual bool isFeasible(Real eq_tol = 1e-2, Real ineq_tol = 1e-14) const {

        for (auto &c : m_evaluatedConstraints) {
            SField result = c->values;

            if (c->type == ConstraintType::EQUALITY) {
                if (std::abs(result.maxMag()) > eq_tol)
                    return false;
            }
            else {
                if (result.min() < -ineq_tol)
                    return false;
            }

        }

        return true;
    }

    bool isParametric() const { return m_parametricOptimization; }

    virtual bool shouldReport() const = 0;

    virtual ~IterateBase() { }

protected:
    // Current params
    std::vector<Real> m_params;

    // Whether we're running a lower-dimensional parametric optimization where
    // gradients must be computed by taking inner products with the shape
    // velocity fields. This also determines whether gradients are printed by
    // the writeDescription method.
    bool m_parametricOptimization;

    // Filled out by subclass
    std::vector<EOTPtr> m_evaluatedObjectiveTerms;
    std::vector< ECPtr> m_evaluatedConstraints;
};

}

#endif /* end of include guard: ITERATEBASE_HH */
