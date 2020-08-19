////////////////////////////////////////////////////////////////////////////////
// Constraint.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Organizes nonlinear inequality and equality constraints and their
//      derivatives.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/07/2016 17:05:05
////////////////////////////////////////////////////////////////////////////////
#ifndef CONSTRAINT_HH
#define CONSTRAINT_HH

#include <MeshFEM/Future.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/OneForm.hh>
#include <MeshFEM/Flattening.hh>
#include <stdexcept>

namespace PatternOptimization {

enum class ConstraintType {
    INEQUALITY, EQUALITY
};

// Vector-valued equality/inequality constraint
template<size_t N>
struct Constraint {
    using Jacobian = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    using  OForm = ScalarOneForm<N>;

    // If type == INEQUALITY, the constraint is "evaluate() >= 0"
    // If type ==   EQUALITY, the constraint is "evaluate() == 0"
    ConstraintType type;

    Constraint(ConstraintType t) : type(t) { }

    virtual ScalarField<Real> evaluate() const = 0;

    virtual size_t dimension() const { return m_componentDifferentials.size(); }

    // Partial derivatives with respect to pattern parameters inducing boundary
    // shape velocities bdrySVels
    virtual Jacobian jacobian(const std::vector<VectorField<Real, N>> &bdrySVels) const {
        const size_t  m = dimension();
        const size_t np = bdrySVels.size();
        Jacobian result(m, np);
        for (size_t c = 0; c < m; ++c) {
            for (size_t p = 0; p < np; ++p)
                result(c, p) = m_componentDifferentials[c][bdrySVels[p]];
        }
        return result;
    }

    virtual void writeContinuousGradientInfo(std::ostream &/*os*/, const std::string &/*name*/) const {
        // TODO: non-parametric gradient norm info:
        // M_norm(steepestDescent), since steepest descent is the Riesz representative of
        // the differential, and we want it's norm. This ends up being
        // sqrt(g^T M^-1 g) where g is the differential
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const {
        os << name << " violation:\t";
        evaluate().print(os, "", "", "", "\t");
        os << std::endl;
        writeContinuousGradientInfo(os, name);
    }

    virtual void writeFields(MSHFieldWriter &/*writer*/) const { }

    virtual ~Constraint() { }
protected:
    std::vector<OForm> m_componentDifferentials;
};

} // end namespace PatternOptimization


#endif /* end of include guard: CONSTRAINT_HH */
