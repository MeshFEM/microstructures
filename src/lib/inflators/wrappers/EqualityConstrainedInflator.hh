////////////////////////////////////////////////////////////////////////////////
// EqualityConstrainedInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Applies linear parameter equality constraints to an inflator by constructing
//  a set of reduced variables.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  05/22/2016 13:53:44
////////////////////////////////////////////////////////////////////////////////
#ifndef EQUALITYCONSTRAINEDINFLATOR_HH
#define EQUALITYCONSTRAINEDINFLATOR_HH

#include "../Inflator.hh"

#include <memory>
#include <utility>
#include <vector>
#include <string>
#include <limits>

class ParameterConstraint;

// Note: we don't use Inflator<N>::{m_elements,m_vertices}, but rather access
// the geometry from m_infl;
template<size_t N>
class EqualityConstrainedInflator : public Inflator<N> {
public:
    EqualityConstrainedInflator(std::unique_ptr<Inflator<N>> infl,
                                const std::vector<std::string> &constraintStrings);

    EqualityConstrainedInflator(std::unique_ptr<Inflator<N>> infl,
                                const std::vector<ParameterConstraint> &constraints);

    ////////////////////////////////////////////////////////////////////////////
    // Geometry access (dimension agnostic)
    ////////////////////////////////////////////////////////////////////////////
    virtual const std::vector<MeshIO::IOElement> &elements() const override { return m_infl->elements(); }
    virtual const std::vector<MeshIO::IOVertex>  &vertices() const override { return m_infl->vertices(); }
    virtual void clear() override { m_infl->clear(); }

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override { m_infl->inflate(fullParametersForReduced(params)); }
public:

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    // Translate shape velocities from full to reduced parameters (chain rule)
    // (Effectively apply the transpose of the change of variables matrix.)
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const override;

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return m_infl->isParametric(); }
    virtual size_t numParameters() const override { return m_indepParamInfluence.size(); }

    // If only one full parameter is controlled, p's type is given by inflator.
    // Otherwise, p is a metaparameter
    virtual ParameterType parameterType(size_t p) const override {
        const auto &influence = m_indepParamInfluence.at(p);
        assert(influence.size() > 0);
        if (influence.size() == 1)
            return m_infl->parameterType(influence[0].first);
        return ParameterType::Meta;
    }

    virtual bool isPrintable(const std::vector<Real> &params) override { return m_infl->isPrintable(fullParametersForReduced(params)); }

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override { m_infl->loadMeshingOptions(moptsPath); }

    virtual void setMaxElementVolume(Real maxElementVol) override { m_infl->setMaxElementVolume(maxElementVol); }
    virtual Real getMaxElementVolume() const             override { return m_infl->getMaxElementVolume(); }
    virtual void setReflectiveInflator(bool use)         override { m_infl->setReflectiveInflator(use); }
    virtual void setDumpSurfaceMesh(bool dump = true)    override { m_infl->setDumpSurfaceMesh(dump); }
    virtual void configureSubdivision(const std::string &algorithm, size_t levels) override { m_infl->configureSubdivision(algorithm, levels); }
    virtual void setOrthoBaseCell(bool ortho)            override { m_infl->setOrthoBaseCell(ortho); }

    ////////////////////////////////////////////////////////////////////////////
    // EqualityConstrainedInflator-specific
    ////////////////////////////////////////////////////////////////////////////
    // Effectively apply the change of variables matrix.
    std::vector<Real> fullParametersForReduced(const std::vector<Real> &rparams);

    ~EqualityConstrainedInflator() { }

private:
    void m_setConstraints(const std::vector<ParameterConstraint> &constraints);

    ////////////////////////////////////////////////////////////////////////////
    // Data members
    ////////////////////////////////////////////////////////////////////////////
    std::unique_ptr<Inflator<N>> m_infl;

    // Sparse representation of the effect (contribution) of each (independent)
    // parameter on the full parameters
    //      m_indepParamInfluence[p] is a list of (full_param_idx, coeff)
    //          where each entry means parameter p contributes
    //          reduced_params[p] * coeff to full_param_idx
    std::vector<std::vector<std::pair<size_t, Real>>> m_indepParamInfluence;

    // The constant offset of each full (output) parameter. For the full
    // parameters corresponding to each reduced (independent) parameter, this is
    // 0. For the full parameters that are dependent, e.g. p1 in:
    //      p1 = a * p2 + b
    // this is b.
    std::vector<Real> m_paramConstOffset;
};

#endif /* end of include guard: EQUALITYCONSTRAINEDINFLATOR_HH */
