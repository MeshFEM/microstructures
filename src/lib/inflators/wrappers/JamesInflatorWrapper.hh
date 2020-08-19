#ifndef JAMESINFLATORWRAPPER_HH
#define JAMESINFLATORWRAPPER_HH

#include "../Inflator.hh"
#include <memory>
#include <MeshFEM/EdgeFields.hh>

#define JAMES_INFLATOR_ENABLED HAS_PYMESH

#if JAMES_INFLATOR_ENABLED

// Forward declare some of James' types
namespace PyMesh {
class PeriodicExploration;
}

class JamesInflatorWrapper : public Inflator<3> {
public:
    JamesInflatorWrapper(const std::string &wireMeshPath,
             Real cell_size = 5.0, Real default_thickness = 0.5 * sqrt(2),
             bool isotropic_params = false, bool vertex_thickness = false);

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override;
public:
    using Inflator<3>::inflate; // Don't hide primary interface with what follows!
    // Not part of the Inflator interface
    void inflate(const std::string &dofFile);

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    virtual std::vector<VectorField<Real, 3>> volumeShapeVelocities() const override;

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override;

    virtual ParameterType parameterType(size_t p) const override;

    // the printability check actually modifies the inflator's internal state,
    // so we can't mark this const.
    // Also, we must change the current dofs to run the printability test
    virtual bool isPrintable(const std::vector<Real> &params) override;

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void setMaxElementVolume(Real maxElementVol) override { m_maxElementVol = maxElementVol; }
    virtual Real getMaxElementVolume() const override { return m_maxElementVol; }
    virtual void configureSubdivision(const std::string &algorithm, size_t levels) override;

    void setReflectiveInflator(bool use) override { m_useReflectiveInflator = use; }
    void setDumpSurfaceMesh(bool dump = true) override { m_dumpSurfaceMesh = dump; }
    
    ////////////////////////////////////////////////////////////////////////////
    // DoF Stuff (not part of Inflator interface)
    ////////////////////////////////////////////////////////////////////////////
    // Configure automatic logging of every set of inflated DoF parameters.
    // If pathPrefix is nonempty, a dof file will be written at
    // pathPrefix_$inflationNumber.dof
    void setDoFOutputPrefix(const std::string &pathPrefix) { m_dofOutPathPrefix = pathPrefix; }
    // Note, overwrites the dofs in m_inflator
    void loadPatternDoFs(const std::string &path, std::vector<Real> &params);

    // Write parameters in James' DoF format.
    void writePatternDoFs(const std::string &path,
                          const std::vector<Real> &params);

    // Out-of-line destructor needed due to incompletely-typed members.
    virtual ~JamesInflatorWrapper();
private:
    std::unique_ptr<PyMesh::PeriodicExploration> m_inflator;

    std::string m_dofOutPathPrefix;
    size_t m_inflationCount = 0;

    Real m_maxElementVol;
    bool m_useReflectiveInflator = true;

    // Used for debugging when tetgen fails.
    bool m_dumpSurfaceMesh = false;

    // Inflate the DoFs already stored in the inflator.
    // (optionally written to surface_debug.msh)
    void m_inflate_dofs();
};

#else // !JAMES_INFLATOR_ENABLED
class JamesInflatorWrapper : public Inflator<3> {
public:
    JamesInflatorWrapper(const std::string &/* wireMeshPath */,
             Real /* cell_size */ = 5.0, Real /* default_thickness */ = 0.5 * sqrt(2),
             bool /* isotropic_params */ = false, bool /* vertex_thickness */ = false) {
        throw std::runtime_error("James' inflator is disabled (check HAS_PYMESH)");
    }

    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override { return 0; }
    virtual ParameterType parameterType(size_t /* p */) const override { return ParameterType::Thickness; }
    // 2D is always printable.
    virtual bool isPrintable(const std::vector<Real> &/* params */) override { return true; }
    virtual ~JamesInflatorWrapper() { }
private:
    virtual void m_inflate(const std::vector<Real> &/* params */) override { }
};

#endif // JAMES_INFLATOR_ENABLED

#endif /* end of include guard: JAMESINFLATORWRAPPER_HH */
