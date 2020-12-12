#include "LuigiInflatorWrapper.hh"

#if LUIGI_INFLATOR_ENABLED

#include <luigi_wire_inflator/src/WireInflator2D.h>

#include <MeshFEM/Future.hh>

////////////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////////////
LuigiInflatorWrapper::LuigiInflatorWrapper(const std::string &wireMeshPath, const int symmetryMode) {
    // symmetryMode < 0 reverts to Luigi's symmetry parameters (WireMesh2D)
    if (symmetryMode < 0) m_inflator = WireInflator2D::construct<WireMesh2D       >(wireMeshPath);
    else                  m_inflator = WireInflator2D::construct<WireMesh2DMorteza>(wireMeshPath, symmetryMode);

    // Translate parameter operation types into our ParameterType
    for (const ParameterOperation &op : m_inflator->getParameterOperations()) {
        switch (op.type) {
            case ParameterOperation::Radius:
                m_paramTypes.push_back(ParameterType::Thickness);
                break;
            case ParameterOperation::Translation:
                m_paramTypes.push_back(ParameterType::Offset);
                break;
            default: assert(false);
        }
    }

    m_tparams = Future::make_unique<TessellationParameters>();

    setMaxElementVolume(0.0001);
}

////////////////////////////////////////////////////////////////////////////////
// Inflation
////////////////////////////////////////////////////////////////////////////////
void LuigiInflatorWrapper::m_inflate(const std::vector<Real> &params) {
    const size_t np = numParameters();
    assert(params.size() == np);
    CellParameters p_params = m_inflator->createParameters();
    for (size_t i = 0; i < np; ++i)
        p_params.parameter(i) = params[i];

    if (!m_inflator->parametersValid(p_params))
        throw runtime_error("Invalid parameters specified.");
    WireInflator2D::OutMeshType inflatedMesh;
    m_inflator->generatePattern(p_params, *m_tparams, inflatedMesh);

    // Convert to MeshIO format
    this->clear();
    for (const auto &p : inflatedMesh.nodes)
        m_vertices.emplace_back(p[0], p[1], 0);
    for (const auto &e : inflatedMesh.elements)
        m_elements.emplace_back(e[0], e[1], e[2]);

    const size_t nv = m_vertices.size();

    // Convert vertex velocity fields into our format.
    assert(inflatedMesh.vertex_velocities.size() == np);
    m_vvels.assign(np, VectorField<Real, 2>(nv));
    for (size_t p = 0; p < np; ++p) {
        assert(inflatedMesh.vertex_velocities[p].size() == nv);
        for (size_t vi = 0; vi < nv; ++vi) {
            m_vvels[p](vi)[0] = inflatedMesh.vertex_velocities[p][vi][0];
            m_vvels[p](vi)[1] = inflatedMesh.vertex_velocities[p][vi][1];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Configuration
////////////////////////////////////////////////////////////////////////////////
void LuigiInflatorWrapper::setMaxElementVolume(Real maxElementVol) { m_tparams->max_area = maxElementVol; }
Real LuigiInflatorWrapper::getMaxElementVolume() const { return m_tparams->max_area; }

////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////
LuigiInflatorWrapper::~LuigiInflatorWrapper() { }

#endif // LUIGI_INFLATOR_ENABLED
