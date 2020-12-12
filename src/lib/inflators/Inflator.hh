////////////////////////////////////////////////////////////////////////////////
// Inflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Abstract base class for inflators.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/10/2015 22:16:39
////////////////////////////////////////////////////////////////////////////////
#ifndef INFLATOR_HH
#define INFLATOR_HH

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/Functions.hh>
#include <MeshFEM/BaseCellType.hh>
#include <MeshFEM/Fields.hh>
#include <isosurface_inflator/MeshingOptions.hh>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Meta parameters are for EqualityConstrainedInflator (custom types allow you to use the optimization with different inflators)
enum class ParameterType { Thickness, Offset, Blending, Meta, Custom1, Custom2, Custom3, Custom4, Custom5, Custom6, Custom7, Custom8 };

inline std::string parameterTypeString(const ParameterType &type) {
    switch(type) {
        case ParameterType::Thickness: return "Thickness";
        case ParameterType::Offset:    return "Offset";
        case ParameterType::Blending:  return "Blending";
        case ParameterType::Meta:      return "Meta";
        case ParameterType::Custom1:      return "Custom1";
        case ParameterType::Custom2:      return "Custom2";
        case ParameterType::Custom3:      return "Custom3";
        case ParameterType::Custom4:      return "Custom4";
        case ParameterType::Custom5:      return "Custom5";
        case ParameterType::Custom6:      return "Custom6";
        case ParameterType::Custom7:      return "Custom7";
        case ParameterType::Custom8:      return "Custom8";
        default:                       return "Invalid";
    }
}

// The dimension-independent part of the interface
class InflatorBase {
public:
    virtual size_t dimension() const = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
    void inflate(const std::vector<Real> &params) {
        if (m_inflationDumpPath.size()) {
            // also be more verbose when dumping
            std::cerr << "INFLATING PARAMS:";
            std::cerr << std::setprecision(19);
            for (Real p : params) std::cerr << "\t" << p;
            std::cerr << std::endl;
        }

        try { this->m_inflate(params); }
        catch (const std::exception &e) {
            std::cerr << std::setprecision(20);
            std::cerr << "Exception \"" << e.what() << "\" while inflating parameters:" << std::endl;
            for (size_t i = 0; i < params.size(); ++i) std::cerr << params[i] << "\t";
            std::cerr << std::endl;
            throw;
        }
        if (m_inflationDumpPath.size())
            MeshIO::save(m_inflationDumpPath, vertices(), elements());
    }

    // Configure to dump the output geometry immediately after inflation
    void setInflationDumpPath(const std::string &path) { m_inflationDumpPath = path; }
    void disableInflationDump() { m_inflationDumpPath.clear(); }

    ////////////////////////////////////////////////////////////////////////////
    // Geometry access (dimension agnostic)
    ////////////////////////////////////////////////////////////////////////////
    virtual const std::vector<MeshIO::IOElement> &elements() const { return m_elements; }
    virtual const std::vector<MeshIO::IOVertex>  &vertices() const { return m_vertices; }
    virtual void clear() { m_elements.clear(), m_vertices.clear(); }

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const = 0;
    virtual size_t numParameters() const = 0;
    virtual ParameterType parameterType(size_t p) const = 0;
    virtual std::vector<Real> defaultParameters() const { throw std::runtime_error("This inflator doesn't implement default parameter query."); }
    // James' inflator's printability check modifies state :(
    virtual bool isPrintable(const std::vector<Real> &params) = 0;

    // Self-supporting constraints evaluated at the current parameters.
    virtual Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    selfSupportingConstraints(const std::vector<double> &/* params */) const {
        throw std::runtime_error("This inflator doesn't implement self-supporting constraints.");
    }

    virtual bool hasOrthotropicSymmetry() const {
        std::cerr << "Warning: inflator doesn't implement hasOrthotropicSymmetry(); assuming true" << std::endl;
        return true;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void loadMeshingOptions(const std::string &/* moptsPath */) {
        throw std::runtime_error("This inflator doesn't support meshing options files");
    }
    virtual MeshingOptions &meshingOptions() {
        throw std::runtime_error("This inflator doesn't support MeshingOptions");
    }

    virtual void setMaxElementVolume(Real /* maxElementVol */) { throw std::runtime_error("Unimplemented"); }
    virtual Real getMaxElementVolume() const                   { throw std::runtime_error("Unimplemented"); }
    virtual void setReflectiveInflator(bool /* use */)         { throw std::runtime_error("Unimplemented"); }
    virtual void setDumpSurfaceMesh(bool dump = true)          { if (dump) throw std::runtime_error("This inflator does not support surface meshing."); }

    // Triply periodic base cell is used by default, but we can use the
    // orhotropic base cell instead:
    virtual void setOrthoBaseCell(bool ortho)                  { if (ortho) throw std::runtime_error("This inflator does not support orthotropic base cells."); }
    virtual BaseCellType baseCellType() const { return BaseCellType::TriplyPeriodic; }

    virtual void configureSubdivision(const std::string &/* algorithm */, size_t levels) {
        if (levels != 0)
            throw std::runtime_error("This inflator doesn't support subdivision");
    }

    virtual void update() {
        // Most inflators do not update!! However, some of them, like Slimflator, can use information of current
        // solution to update their structures. In case of Slimflator, the inflator needs a reference structure
        // to work (cause it uses SLIM to move internal vertices, following boundary perturbation). The algorithm
        // works better if the current reference shape is updated with the best solution found so far.
    }

    virtual ~InflatorBase() { }

protected:
    std::vector<MeshIO::IOElement> m_elements;
    std::vector<MeshIO::IOVertex>  m_vertices;
    std::string m_inflationDumpPath;
private:
    virtual void m_inflate(const std::vector<Real> &params) = 0;
};

// The dimension-dependent part of the interface
template<size_t _N>
class Inflator : public InflatorBase {
public:
    static constexpr size_t N = _N;

    virtual size_t dimension() const override { return N; }

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation/access
    // All inflators must implement per-volume-vertex parameter shape
    // velocities (volumeShapeVelocities). Then a per-boundary-vertex shape
    // velocity is extracted by ignoring the internal values.
    //
    // This is done because the inflators themselves might not know which are
    // boundary vertices (or might not have the same boundary vertex indexing)
    ////////////////////////////////////////////////////////////////////////////
    virtual std::vector<VectorField<Real, N>> volumeShapeVelocities() const {
        throw std::runtime_error("Shape velocities only supported by parametric inflators.");
    }

    // For non-parametric inflators (where individual vertex positions are
    // variables), it's inefficient to work with parameter velocities. However,
    // we can efficiently extract the parameter vector corresponding to a boundary
    // descent direction/differential form.
    virtual ScalarField<Real> paramsFromBoundaryVField(const VectorField<Real, N> &/*bdryVField*/) const {
        throw std::runtime_error("Boundary => param descent conversion only supported by non-parametric inflators.");
    }

    // Compute the per-boundary-vertex velocity field corresponding to a parameter
    // velocity field.
    virtual VectorField<Real, N> boundaryVFieldFromParams(const ScalarField<Real> &/*pvel*/) const {
        throw std::runtime_error("boundaryVFieldFromParams() unimplemented by this inflator");
    }

    // Extract (true) boundary shape velocities from volume shape velocities.
    // Non-virtual!
    template<class _FEMMesh>
    std::vector<VectorField<Real, N>> shapeVelocities(const _FEMMesh &mesh) const {
        auto vvels = volumeShapeVelocities();
        const size_t nvels = vvels.size();
        std::vector<VectorField<Real, N>> result(nvels, VectorField<Real, N>(mesh.numBoundaryVertices()));

        std::vector<bool> isTrueBoundaryVertex(mesh.numBoundaryVertices(), false);
        for (auto be : mesh.boundaryElements()) {
            if (be->isInternal) continue;
            for (auto bv : be.vertices())
                isTrueBoundaryVertex.at(bv.index()) = true;
        }

        for (size_t i = 0; i < nvels; ++i) {
            result[i].clear();
            for (auto bv : mesh.boundaryVertices()) {
                if (!isTrueBoundaryVertex.at(bv.index())) continue;
                result[i](bv.index()) = vvels[i](bv.volumeVertex().index());
            }
        }
        return result;
    }

    // The meshing cell box.
    virtual BBox<Vector3D> meshingCell() {
        throw std::runtime_error("meshingCell() unimplemented by this inflator");
    }

    virtual ~Inflator() { }
};

#endif /* end of include guard: INFLATOR_HH */
