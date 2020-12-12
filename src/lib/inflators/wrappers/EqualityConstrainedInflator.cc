#include "EqualityConstrainedInflator.hh"

#include "../rref.h"
#include "../ParameterConstraint.hh"

using namespace std;

template<size_t N>
EqualityConstrainedInflator<N>::
EqualityConstrainedInflator(std::unique_ptr<Inflator<N>> infl,
                            const std::vector<std::string> &constraintStrings)
    : m_infl(std::move(infl))
{
    const size_t fullNumParams = m_infl->numParameters();
    std::vector<ParameterConstraint> constraints;
    for (const auto &cs : constraintStrings)
        constraints.emplace_back(fullNumParams, cs);
    m_setConstraints(constraints);
}

template<size_t N>
EqualityConstrainedInflator<N>::
EqualityConstrainedInflator(std::unique_ptr<Inflator<N>> infl,
                            const std::vector<ParameterConstraint> &constraints)
        : m_infl(std::move(infl))
{
    m_setConstraints(constraints);
}

////////////////////////////////////////////////////////////////////////////
// Shape velocity computation
////////////////////////////////////////////////////////////////////////////
// Translate shape velocities from full to reduced parameters (chain rule)
// (Effectively apply the transpose of the change of variables matrix.)
template<size_t N>
std::vector<VectorField<Real, N>> 
EqualityConstrainedInflator<N>::
volumeShapeVelocities() const {
    std::vector<VectorField<Real, N>> fullVelocities = m_infl->volumeShapeVelocities();

    const size_t nip = numParameters();
    std::vector<VectorField<Real, N>> result(nip);
    VectorField<Real, N> tmp;
    for (size_t ip = 0; ip < nip; ++ip) {
        auto &vvel = result[ip];
        vvel.resizeDomain(vertices().size());
        vvel.clear();
        for (const std::pair<size_t, Real> &contrib : m_indepParamInfluence[ip]) {
            tmp = fullVelocities.at(contrib.first);
            tmp *= contrib.second;
            vvel += tmp;
        }
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////
// EqualityConstrainedInflator-specific
////////////////////////////////////////////////////////////////////////////
// Effectively apply the change of variables matrix.
template<size_t N>
std::vector<Real>
EqualityConstrainedInflator<N>::
fullParametersForReduced(const std::vector<Real> &rparams) {
    std::vector<Real> fullParams = m_paramConstOffset;

    const size_t nip = numParameters();
    assert(rparams.size() == nip);

    for (size_t ip = 0; ip < nip; ++ip) {
        for (const std::pair<size_t, Real> &contrib : m_indepParamInfluence[ip])
            fullParams.at(contrib.first) += contrib.second * rparams.at(ip);
    }

    return fullParams;
}

////////////////////////////////////////////////////////////////////////////////
// Private methods
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
void
EqualityConstrainedInflator<N>::
m_setConstraints(const std::vector<ParameterConstraint> &constraints) {
    if (!m_infl->isParametric())
        std::cerr << "WARNING: equality constraints intended only for parametric inflators" << std::endl;
    // Impose all equality constraints exactly:
    // Solve for a set of independent variables and express the
    // dependent variables in terms of these.
    const size_t fullNumParams = m_infl->numParameters();
    const size_t expectedCols = fullNumParams + 1;
    vector<vector<Real>> augmentedConstraintSystem;
    for (const ParameterConstraint &pc : constraints) {
        if (!pc.isEqualityConstraint()) continue;
        augmentedConstraintSystem.emplace_back(pc.augmentedRow());
        if (augmentedConstraintSystem.back().size() != expectedCols)
            throw std::runtime_error("Invalid constraint row");
    }
    if (augmentedConstraintSystem.size() != constraints.size())
        std::cerr << "WARNING: Only equality constraints are currently implemented." << std::endl;

    // depIdx[p]: index of full parameter corresponding to dependent var p
    // depRow[p]: row of augmentedConstraintSystem specifying contributions
    //            to dependent var p (indep var coeffs and constant offset)
    std::vector<size_t> depIdx, depRow;

    // By construction, the augmented column will be made independent.
    // If this doesn't happen, the system is inconsistent.
    rref(augmentedConstraintSystem, depIdx, depRow);

    // Determine the independent variables (complement of dep)
    std::vector<size_t> indepIdx;
    {
        size_t currDepVar = 0;
        for (size_t i = 0; i < expectedCols; ++i) {
            if ((currDepVar < depIdx.size()) && (i == depIdx[currDepVar]))
                ++currDepVar;
            else indepIdx.push_back(i);
        }
        assert(currDepVar == depIdx.size());

        // The augmented column vector must have been made independent.
        if ((indepIdx.size() == 0) || (indepIdx.back() != expectedCols - 1))
            throw std::runtime_error("Couldn't make RHS independent (inconsistent constraint system)");

        // The last column isn't actually a variable
        indepIdx.pop_back();
    }

    // Determine the dependent parameter specified by each row
    std::vector<size_t> depParamForRow(augmentedConstraintSystem.size(),
                                       std::numeric_limits<size_t>::max());
    assert(depRow.size() == depIdx.size());
    for (size_t dpi = 0; dpi < depRow.size(); ++dpi)
        depParamForRow.at(depRow[dpi]) = depIdx[dpi];

    // Determine the effect of each independent variable
    m_indepParamInfluence.clear();
    m_indepParamInfluence.reserve(indepIdx.size());
    for (size_t i = 0; i < indepIdx.size(); ++i) {
        size_t indi = indepIdx[i];
        // Independent parameters control their corresponding full param
        m_indepParamInfluence[i].push_back({indi, 1.0});
        // And full params corresponding to affected dependent params
        for (size_t rowi = 0; rowi < augmentedConstraintSystem.size(); ++rowi) {
            Real coeff = augmentedConstraintSystem[rowi][indi];
            if (std::abs(coeff) > 1e-8) {
                size_t dep = depParamForRow[rowi];
                assert(dep < fullNumParams);
                m_indepParamInfluence[i].push_back({dep, coeff});
            }
        }
    }

    // Determine the constant offset for each full parameter
    // (This is 0 for the independent variables)
    m_paramConstOffset.assign(fullNumParams, 0.0);
    for (size_t rowi = 0; rowi < augmentedConstraintSystem.size(); ++rowi) {
        size_t dp = depParamForRow[rowi];
        if (dp < fullNumParams) m_paramConstOffset.at(dp) = augmentedConstraintSystem[rowi].back();
        else {
            // Rows not specifying dependent parameters better be all-zeros.
            assert(dp == std::numeric_limits<size_t>::max());
            assert(std::abs(augmentedConstraintSystem[rowi].back()) < 1e-8);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: 2D and 3D inflators.
////////////////////////////////////////////////////////////////////////////////
template class EqualityConstrainedInflator<2>;
template class EqualityConstrainedInflator<3>;
