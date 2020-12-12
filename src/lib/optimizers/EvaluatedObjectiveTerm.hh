////////////////////////////////////////////////////////////////////////////////
// EvaluatedObjectiveTerm.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Non-templated data structure to store the quantities needed by the
//      various optimization algorithms (objective and gradient info).
//
//      This is the primary interface between (parametric) optimization
//      algorithms and the optimization iterates.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  04/27/2016 11:08:21
////////////////////////////////////////////////////////////////////////////////
#ifndef EVALUATEDOBJECTIVETERM_HH
#define EVALUATEDOBJECTIVETERM_HH

#include <MeshFEM/Types.hh>
#include <MeshFEM/EdgeFields.hh>
#include <vector>
#include <memory>
#include <string>

namespace PatternOptimization {

struct EvaluatedObjectiveTerm {
    std::string name;

    Real normalization = 1.0;
    Real weight = 1.0;

    void setValue(Real val)       { m_value = val; }

    Real normalizedWeight() const { return normalization * weight; }
    Real contribution()     const { return value() * normalizedWeight(); }

    ScalarField<Real> gradContribution() const {
        ScalarField<Real> result(gradp);
        result *= normalizedWeight();
        return result;
    }

    ScalarField<Real> descentContribution() const {
        ScalarField<Real> result(descentp);
        result *= normalizedWeight();
        return result;
    }

    virtual Real value() const {
        Real result = m_value;
        if (isEstimating()) result += gradp.values().dot(m_estimationDeltaP.values());
        return result;
    }

    bool isEstimating() const { return m_estimationDeltaP.domainSize() != 0; }

    void setEstimateWithDeltaParams(const ScalarField<Real> &delta) {
        if (delta.domainSize() != gradp.domainSize())
            throw std::runtime_error("Illegal delta param size");
        m_estimationDeltaP = delta;
    }

    void disableEstimation() { m_estimationDeltaP.resizeDomain(0); }

    // Note: also outputs estimated term, if we're estimating
    virtual void writeGradientDescription(std::ostream &os, bool isParametric) {
        if (isEstimating())
            os << "estimated " << name << ":\t" << value() << std::endl;
        if (isParametric) {
            os << "grad_p " << name << ":\t";
            gradp.print(os, "", "", "", "\t");
            os << std::endl;
        }
        os << "||grad_p " << name << "||:\t" << gradp.norm() << std::endl;
    }

    // All objective terms provide gradient information.
    ScalarField<Real> gradp, descentp;
    virtual ~EvaluatedObjectiveTerm() { }
protected:
    Real m_value;
    ScalarField<Real> m_estimationDeltaP;
};

// Only Nonlinear Least Squares terms support computation of residuals and
// Jacobians.
struct EvaluatedObjectiveTermNLLS : public EvaluatedObjectiveTerm {
    Real residual(size_t r) const { return residualComponents[r]; }
    Real jacobian(size_t r, size_t p) const { return jacobianComponents(r, p); }

    size_t numResiduals() const { return residualComponents.domainSize(); }

    virtual Real value() const {
        if (isEstimating()) {
            // Note: linearly approximating the least squares objective is
            // different from linearly approximating the residual, and we prefer
            // to do the latter.
            ScalarField<Real> estimatedResidual(residualComponents);
            estimatedResidual += ScalarField<Real>(jacobianComponents * this->m_estimationDeltaP.values().transpose());
            Real result = 0;
            for (size_t ri = 0; ri < residualComponents.domainSize(); ++ri)
                result += residualComponents[ri] * residualComponents[ri];
            return result;
        }
        return m_value;
    }

    ScalarField<Real> residualComponents;
    Eigen::MatrixXd   jacobianComponents;
    virtual ~EvaluatedObjectiveTermNLLS() { }
};

}

#endif /* end of include guard: EVALUATEDOBJECTIVETERM_HH */
