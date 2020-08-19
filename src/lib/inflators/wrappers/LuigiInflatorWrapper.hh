#ifndef LUIGIINFLATORWRAPPER_HH
#define LUIGIINFLATORWRAPPER_HH

#include "../Inflator.hh"

#include <vector>
#include <memory>
#include <string>

#include <MeshFEM/EdgeFields.hh>

#define LUIGI_INFLATOR_ENABLED (HAS_VCGLIB && HAS_CLIPPER)

#if LUIGI_INFLATOR_ENABLED

// Forward-declare some of Luigi's types.
class WireInflator2D;
struct TessellationParameters;

class LuigiInflatorWrapper : public Inflator<2> {
public:
    // MHS JUL14, 2015: optional "symmetryMode" is passed to WireMesh2DMorteza
    LuigiInflatorWrapper(const std::string &wireMeshPath,
                         const int symmetryMode = -1);

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override;
public:

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    virtual std::vector<VectorField<Real, 2>> volumeShapeVelocities() const override { return m_vvels; }

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override { return m_paramTypes.size(); }
    virtual ParameterType parameterType(size_t p) const override { return m_paramTypes.at(p); }
    // 2D is always printable.
    virtual bool isPrintable(const std::vector<Real> &/* params */) override { return true; }

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void setMaxElementVolume(Real maxElementVol) override;
    virtual Real getMaxElementVolume() const override;
    virtual void setReflectiveInflator(bool use) override { if (use) throw std::runtime_error("Luigi's inflator is not reflective"); }

    // Out-of-line destructor needed due to incompletely-typed members.
    virtual ~LuigiInflatorWrapper();

private:
    std::shared_ptr<WireInflator2D> m_inflator;
    std::unique_ptr<TessellationParameters> m_tparams;

    std::vector<ParameterType> m_paramTypes;
    std::vector<VectorField<Real, 2>> m_vvels; // velocity field induced by each param
};

#else // !!LUIGI_INFLATOR_ENABLED
class LuigiInflatorWrapper : public Inflator<2> {
public:
    LuigiInflatorWrapper(const std::string &/* wireMeshPath */, const int /* symmetryMode */ = -1) {
        throw std::runtime_error("Luigi's inflator is disabled (check HAS_VCGLIB and HAS_CLIPPER)");
    }

    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override { return 0; }
    virtual ParameterType parameterType(size_t /* p */) const override { return ParameterType::Thickness; }
    // 2D is always printable.
    virtual bool isPrintable(const std::vector<Real> &/* params */) override { return true; }

    virtual ~LuigiInflatorWrapper() { }
private:
    virtual void m_inflate(const std::vector<Real> &/* params */) override { }
};
#endif // LUIGI_INFLATOR_ENABLED

#endif /* end of include guard: LUIGIINFLATORWRAPPER_HH */
