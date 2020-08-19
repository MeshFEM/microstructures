#include "WireInflator2DFacade.h"

#include <iostream>
#include <cassert>
#include <sstream>

#include <Core/Exception.h>

WireInflatorFacade::~WireInflatorFacade() {
    for (size_t row=0; row<m_rows; row++) {
        for (size_t col=1; col<m_cols; col++) {
            if ((*m_p_params)(col, row) != NULL) {
                delete (*m_p_params)(col, row);
            }
        }
    }
}

void WireInflatorFacade::set_dimension(size_t rows, size_t cols) {
    m_rows = rows;
    m_cols = cols;
    m_p_params = ParameterGridPtr(new ParameterGrid(m_cols, m_rows));
    for (size_t row = 0; row < m_rows; row++) {
        for (size_t col=0; col < m_cols; col++) {
            (*m_p_params)(col, row) = NULL;
        }
    }
}

WireInflatorFacade::ParameterType WireInflatorFacade::get_parameter_type(size_t i) {
    const ParameterOperation& param_op = m_inflator->getParameterOperations()[i];
    switch (param_op.type) {
        case ParameterOperation::Radius:
            return THICKNESS;
        case ParameterOperation::Translation:
            return VERTEX_OFFSET;
        default:
            throw NotImplementedError("Unknown parameter type detected.");
    }
}

VectorI WireInflatorFacade::get_affected_vertex_orbit(size_t i) {
    VectorI vertex_orbit;
    const ParameterOperation& param_op = m_inflator->getParameterOperations()[i];
    size_t count = 0;
    switch (param_op.type) {
        case ParameterOperation::Radius:
            vertex_orbit.resize(param_op.nodes.size());
            std::copy(param_op.nodes.begin(), param_op.nodes.end(), vertex_orbit.data());
            break;
        case ParameterOperation::Translation:
            vertex_orbit.resize(param_op.nodes_displ.size());
            for (auto itr : param_op.nodes_displ) {
                vertex_orbit[count] = itr.first;
                count++;
            }
            break;
        default:
            throw NotImplementedError("Unknown parameter type detected.");
    }
    assert(vertex_orbit.size() > 0);
    return vertex_orbit;
}

MatrixF WireInflatorFacade::get_offset_direction(size_t i) {
    const ParameterOperation& param_op = m_inflator->getParameterOperations()[i];
    if (param_op.type == ParameterOperation::Translation) {
        const size_t num_nodes = param_op.nodes_displ.size();
        MatrixF offset_dir(num_nodes, 2);
        size_t count = 0;
        for (auto itr : param_op.nodes_displ) {
            offset_dir.coeffRef(count, 0) = itr.second[0];
            offset_dir.coeffRef(count, 1) = itr.second[1];
            count++;
        }
        return offset_dir;
    } else {
        std::stringstream err_msg;
        err_msg << "Parameter " << i << " is not offset parameter";
        throw RuntimeError(err_msg.str());
    }
}

void WireInflatorFacade::set_parameter(size_t row, size_t col, const VectorF& param) {
    assert(row < m_rows && col < m_cols);
    assert(param.size() == get_num_parameters());
    if (row >= m_rows) {
        std::stringstream err_msg;
        err_msg << "Out of bound: " << row << " >= " << m_rows;
        throw RuntimeError(err_msg.str());
    }
    if (col >= m_cols) {
        std::stringstream err_msg;
        err_msg << "Out of bound: " << col << " >= " << m_cols;
        throw RuntimeError(err_msg.str());
    }

    const size_t num_param = m_inflator->numberOfParameters();
    CellParameters* p = new CellParameters(m_inflator->numberOfParameters());

    for (size_t i=0; i<num_param; i++) {
        p->parameter(i) = param[i];
    }

    (*m_p_params)(col, row) = p;
}

void WireInflatorFacade::generate_pattern_with_guide_mesh(
        const VectorF& vertices, const VectorI& faces,
        const MatrixFr& raw_parameters) {
    const size_t num_cells = raw_parameters.rows();
    const size_t num_params = get_num_parameters();
    assert(raw_parameters.cols() == num_params);

    ParameterVector parameters(num_cells);
    for (size_t i=0; i<num_cells; i++) {
        const VectorF& param = raw_parameters.row(i);
        parameters[i] = CellParameters(num_params);
        for (size_t j=0; j<num_params; j++) {
            parameters[i].parameter(j) = param[j];
        }
    }

    m_inflator->generateQuadsPattern(vertices, faces, parameters, m_t_params, m_mesh);
}

VectorF WireInflatorFacade::get_vertices() {
    const size_t dim = 2;
    typedef WireInflator2D::OutMeshType MeshType;
    MeshType::NodeVector nodes = m_mesh.nodes;
    VectorF vertices(nodes.size() * dim);
    for (size_t i=0; i<nodes.size(); i++) {
        for (size_t j=0; j<dim; j++) {
            vertices[i*dim + j] = nodes[i][j];
        }
    }
    return vertices;
}

VectorI WireInflatorFacade::get_triangles() {
    const size_t dim = 2;
    const size_t nodes_per_elem = 3;
    typedef WireInflator2D::OutMeshType MeshType;
    MeshType::ElementVector elems = m_mesh.elements;

    VectorI faces(elems.size() * nodes_per_elem);
    for (size_t i=0; i<elems.size(); i++) {
        for (size_t j=0; j<nodes_per_elem; j++) {
            faces[i*nodes_per_elem + j] = elems[i][j];
        }
    }
    return faces;
}

MatrixFr WireInflatorFacade::get_boundary_velocity() {
    const size_t dim = 2;
    typedef WireInflator2D::OutMeshType MeshType;

    MatrixFr bd_velocity = MatrixF::Zero(m_mesh.nodes.size(), get_num_parameters());
    size_t count = 0;
    for (auto itr = m_mesh.edge_fields.begin();
            itr != m_mesh.edge_fields.end(); itr++) {
        MeshType::EdgeType edge = itr->first;
        MeshType::Fields edge_velocity = itr->second;
        std::copy(edge_velocity.begin(), edge_velocity.end(),
                bd_velocity.row(edge.first).data());
        std::copy(edge_velocity.begin(), edge_velocity.end(),
                bd_velocity.row(edge.second).data());
        count++;
    }

    return bd_velocity;
}
