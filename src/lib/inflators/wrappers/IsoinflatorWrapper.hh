#ifndef ISOINFLATORWRAPPER_HH
#define ISOINFLATORWRAPPER_HH

#include "../Inflator.hh"

#include <memory>
#include <vector>
#include <string>
#include <stdexcept>

#include <MeshFEM/EdgeFields.hh>
#include <MeshFEM/BaseCellType.hh>

// Forward-declare IsosurfaceInflator
class IsosurfaceInflator;

template<size_t N>
class IsoinflatorWrapper : public Inflator<N> {
public:
    IsoinflatorWrapper(const std::string &wireMeshPath,
                       const std::string &symmetryType, bool vertex_thickness,
                       size_t inflationGraphRadius);

    ////////////////////////////////////////////////////////////////////////////
    // Override geometry access
    ////////////////////////////////////////////////////////////////////////////
    virtual const std::vector<MeshIO::IOElement> &elements() const override;
    virtual const std::vector<MeshIO::IOVertex>  &vertices() const override;
    virtual void clear() override;

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override;
public:

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    // Get a per-vertex perturbation vector field induced by changing
    // each parameter. Internal vertices get 0 velocity.
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const override;

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override    { return true; }
    virtual size_t numParameters() const override;
    virtual ParameterType parameterType(size_t p) const override;
    virtual std::vector<Real> defaultParameters() const override;
    virtual bool isPrintable(const std::vector<Real> &params) override;
    virtual Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
        selfSupportingConstraints(const std::vector<double> &params) const override;

    virtual bool hasOrthotropicSymmetry() const override;

    virtual BBox<Vector3D> meshingCell() override;

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &moptsPath) override;
    virtual void setMaxElementVolume(Real maxElementVol) override;
    virtual void setReflectiveInflator(bool use) override { if (!use) throw std::runtime_error("IsosurfaceInflator is always reflective."); }

    virtual MeshingOptions &meshingOptions() override;

    virtual void setOrthoBaseCell(bool ortho) override;
    virtual BaseCellType baseCellType() const override;

    // Avoid computing normal shape velocities related to all parameters.
    // Gives an important performance advantage when using only the inflator binary (without optimization)
    void disableCheapPostprocess();
    void  enableCheapPostprocess();

    // Out-of-line destructor needed due to incompletely-typed members.
    virtual ~IsoinflatorWrapper();
private:
    std::unique_ptr<IsosurfaceInflator> m_inflator;
};

#endif /* end of include guard: ISOINFLATORWRAPPER_HH */
