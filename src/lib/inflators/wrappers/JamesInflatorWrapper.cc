#include "JamesInflatorWrapper.hh"

#if JAMES_INFLATOR_ENABLED

#include <MeshFEM/Future.hh>
#include <Wires/Interfaces/PeriodicExploration.h>
#include <stdexcept>

using namespace PyMesh;

namespace {
    VectorF ToVectorF(const std::vector<Real> &vec) {
        VectorF result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) result(i) = vec[i];
        return result;
    }
    std::vector<Real> FromVectorF(const VectorF &vec) {
        std::vector<Real> result(vec.rows());
        for (int i = 0; i < vec.rows(); ++i) result[i] = vec(i);
        return result;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////////////
JamesInflatorWrapper::JamesInflatorWrapper(const std::string &wireMeshPath,
             Real cell_size, Real default_thickness,
             bool isotropic_params, bool vertex_thickness) {
    m_inflator = Future::make_unique<PeriodicExploration>(wireMeshPath,
                    cell_size, default_thickness);
    ParameterCommon::TargetType thickness_type =
        vertex_thickness ? ParameterCommon::VERTEX : ParameterCommon::EDGE;
    if (isotropic_params) m_inflator->with_all_isotropic_parameters(thickness_type);
    else                  m_inflator->with_all_parameters(thickness_type);
    setMaxElementVolume(0.0);
}

////////////////////////////////////////////////////////////////////////////////
// Inflation
////////////////////////////////////////////////////////////////////////////////
// Not part of the Inflator interface
void JamesInflatorWrapper::inflate(const std::string &dofFile) {
    m_inflator->load_dofs(dofFile);
    inflate(FromVectorF(m_inflator->get_dofs()));
}

void JamesInflatorWrapper::m_inflate(const std::vector<Real> &params) {
    m_inflate_dofs();
    m_inflator->set_dofs(ToVectorF(params));
}

////////////////////////////////////////////////////////////////////////////////
// Shape velocity computation
////////////////////////////////////////////////////////////////////////////////
std::vector<VectorField<Real, 3>> JamesInflatorWrapper::volumeShapeVelocities() const {
    std::vector<MatrixFr> vertexVelocities = m_inflator->get_shape_velocities();
    const size_t nv = vertices().size();
    size_t nParams = numParameters();
    assert(vertexVelocities.size() == nParams);
    assert(vertexVelocities.at(0).rows() == int(nv));
    std::vector<VectorField<Real, 3>> result(nParams, VectorField<Real, 3>(nv));
    for (size_t p = 0; p < nParams; ++p) {
        for (size_t vi = 0; vi < nv; ++vi)
            result[p](vi) = vertexVelocities[p].row(vi);
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Queries
////////////////////////////////////////////////////////////////////////////////
bool JamesInflatorWrapper::isPrintable(const std::vector<Real> &params) {
    m_inflator->set_dofs(ToVectorF(params));
    return m_inflator->is_printable();
}

size_t JamesInflatorWrapper::numParameters() const { return m_inflator->get_num_dofs(); }

////////////////////////////////////////////////////////////////////////////////
// Configuration
////////////////////////////////////////////////////////////////////////////////
void JamesInflatorWrapper::configureSubdivision(const std::string &algorithm, size_t levels) {
    m_inflator->with_refinement(algorithm, levels);
}

////////////////////////////////////////////////////////////////////////////////
// DoF Stuff (not part of Inflator interface)
////////////////////////////////////////////////////////////////////////////////
ParameterType JamesInflatorWrapper::parameterType(size_t p) const {
    return m_inflator->is_thickness_dof(p) ? ParameterType::Thickness
                                           : ParameterType::Offset;
}

void JamesInflatorWrapper::loadPatternDoFs(const std::string &path, std::vector<Real> &params) {
    m_inflator->load_dofs(path);
    params = FromVectorF(m_inflator->get_dofs());
}

// Write parameters in James' DoF format.
void JamesInflatorWrapper::writePatternDoFs(const std::string &path,
                                            const std::vector<Real> &params) {
    m_inflator->set_dofs(ToVectorF(params));
    m_inflator->save_dofs(path);
}

////////////////////////////////////////////////////////////////////////////////
// Private methods
////////////////////////////////////////////////////////////////////////////////
// Inflate the DoFs already stored in the inflator.
// (optionally written to surface_debug.msh)
void JamesInflatorWrapper::m_inflate_dofs() {
    if (m_dofOutPathPrefix != "") {
        m_inflator->save_dofs(m_dofOutPathPrefix + "_" +
                std::to_string(m_inflationCount) + ".dof");
    }

    m_inflator->periodic_inflate(m_useReflectiveInflator);

    if (m_dumpSurfaceMesh) {
        // Debug surface mesh (for when tetgen fails)
        std::vector<MeshIO::IOElement> triangles;
        std::vector<MeshIO::IOVertex>  vertices;
        MatrixFr verts = m_inflator->get_vertices();
        MatrixIr facs = m_inflator->get_faces();
        for (size_t i = 0; i < size_t(facs.rows()); ++i)
            triangles.emplace_back(facs(i, 0), facs(i, 1), facs(i, 2));
        for (size_t i = 0; i < size_t(verts.rows()); ++i)
            vertices.emplace_back(Point3D(verts.row(i)));
        MeshIO::save("surface_debug.msh", vertices, triangles);
    }

    m_inflator->run_tetgen(m_maxElementVol);

    ++m_inflationCount;

    // Convert to MeshIO format.
    this->clear();
    MatrixIr els = m_inflator->get_voxels();
    for (size_t i = 0; i < size_t(els.rows()); ++i)
        m_elements.emplace_back(els(i, 0), els(i, 1), els(i, 2), els(i, 3));
    MatrixFr verts = m_inflator->get_vertices();
    for (size_t i = 0; i < size_t(verts.rows()); ++i)
        m_vertices.emplace_back(Point3D(verts.row(i)));
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////
JamesInflatorWrapper::~JamesInflatorWrapper() { }

#endif // JAMES_INFLATOR_ENABLED
