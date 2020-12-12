////////////////////////////////////////////////////////////////////////////////
// ObjectiveTerm.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Organizes the individual terms in a compound pattern optimization objective
//  and provides a unified interface for evaluating various types of
//  derivatives.
//
//  Derivative information is represented by boundary velocity one-forms (a
//  boundary vector field that, when dotted with a boundary velocity, gives the
//  change in objective).
//
//  All objective terms store gradients, but Nonlinear Least Squares (NLLS)
//  objectives can also compute residuals and jacobians.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  04/07/2016 16:14:46
////////////////////////////////////////////////////////////////////////////////
#ifndef OBJECTIVETERM_HH
#define OBJECTIVETERM_HH

#include <MeshFEM/Future.hh>
#include <MeshFEM/Flattening.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <stdexcept>

namespace PatternOptimization {

// Base class for objective terms. At a minimum, subclasses must implement:
//    evaluate(): evaluate sub-objective value
// and initialize member variable:
//    m_differential: boundary vector field whose inner product with a boundary
//                    velocity evaluates the shape derivative one-form on this
//                    velocity field
template<size_t N>
struct ObjectiveTerm {
    ObjectiveTerm() { }

    using  OForm = ScalarOneForm<N>;
    using VField = VectorField<Real, N>;
    using SField = ScalarField<Real>;

    virtual Real evaluate() const = 0;

    Real evaluateNormalized() const {
        Real result = evaluate();
        return result * m_normalization;
    }

    // Differential (one-forms) of the sub-objective
    // This is a BOUNDARY differential!
    const OForm &differential() const { return m_differential; }

    // Partial derivatives with respect to pattern parameters inducing boundary
    // shape velocities bdrySVels
    virtual SField gradp(const std::vector<VField> &bdrySVels) const {
        size_t np = bdrySVels.size();
        SField g(np);
        for (size_t p = 0; p < np; ++p)
            g[p] = m_differential[bdrySVels[p]];
        return g;
    }

    // Access weight for combining into a compound objective--this weight is
    // *NOT* applied in the evaluate() and gradient routines.
    void setWeight(Real weight) { m_weight = weight; }
    Real weight() const { return m_weight; }

    void setNormalization(Real n) { m_normalization = n; }
    Real    normalization() const { return m_normalization; }

    Real normalizedWeight() const { return m_weight * m_normalization; }

    virtual void writeContinuousGradientInfo(std::ostream &/*os*/, const std::string &/*name*/) const {
        // TODO: non-parametric gradient norm info:
        // M_norm(steepestDescent), since steepest descent is the Riesz representative of
        // the differential, and we want it's norm. This ends up being
        // sqrt(g^T M^-1 g) where g is the differential
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const {
        os << name << ":\t" << evaluate() << std::endl;
        os << "normalized " << name << ":\t" << evaluateNormalized() << std::endl;
        writeContinuousGradientInfo(os, name);
    }

    virtual void writeFields(MSHFieldWriter &/*writer*/) const { }

    virtual ~ObjectiveTerm() { }
protected:
    Real m_weight = 1.0;
    Real m_normalization = 1.0;
    OForm m_differential;
};

// Abstract base class for Nonlinear Least Squares objective terms.
// On top of ObjectiveTerm's requirements, subclasses must implement:
//      - jacobian
//      - residual
template<size_t N>
struct NLLSObjectiveTerm : public ObjectiveTerm<N> {
    using SField = typename ObjectiveTerm<N>::SField;
    using VField = typename ObjectiveTerm<N>::VField;

    NLLSObjectiveTerm() { }

    // Precompute jacobians with respect to the pattern parameters inducing
    // boundary shape velocities bdrySVels
    virtual Eigen::MatrixXd jacobian(const std::vector<VField> &bdrySVels) const = 0;
    virtual SField residual() const = 0;

    virtual ~NLLSObjectiveTerm() { }
};

}

#endif /* end of include guard: OBJECTIVETERM_HH */
