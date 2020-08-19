////////////////////////////////////////////////////////////////////////////////
// EvaluatedConstraint.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Non-templated data structure to store the quantities needed by the
//      various optimization algorithms (objective and gradient info).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/07/2016 21:22:43
////////////////////////////////////////////////////////////////////////////////
#ifndef EVALUATEDCONSTRAINT_HH
#define EVALUATEDCONSTRAINT_HH

#include "Constraint.hh"
#include <MeshFEM/Types.hh>
#include <vector>
#include <memory>
#include <string>

namespace PatternOptimization {

struct EvaluatedConstraint {
    using Jacobian = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    ConstraintType type;
    std::string name;
    ScalarField<Real> values;
    Jacobian jacobian;

    EvaluatedConstraint(ConstraintType t, const std::string &n,
                        const ScalarField<Real> &v, const Jacobian &j)
        : type(t), name(n), values(v), jacobian(j) { }

    size_t dimension() const { return jacobian.rows(); }

    // Note: also outputs estimated term, if we're estimating
    void writeGradientDescription(std::ostream &os, bool isParametric) const {
        assert(isParametric); // constraints currently only implemented for parametric optimization
        if (isParametric) {
            os << "jacobian " << name << " row norms:";
            for (int i = 0; i < jacobian.rows(); ++i)
                os << '\t' << jacobian.row(i).norm();
            os << std::endl;
            // os << "jacobian " << name << ":\t" << std::endl;
            // os << jacobian << std::endl;
        }
    }
};

}

#endif /* end of include guard: EVALUATEDCONSTRAINT_HH */
