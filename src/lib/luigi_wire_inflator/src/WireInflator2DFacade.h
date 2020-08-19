#pragma once
#include <memory>
#include <string>
#include <vector>

#include "InflatorParameters.h"
#include "WireInflator2D.h"

#include <Core/EigenTypedef.h>

class WireInflatorFacade {
    public:
        enum ParameterType {
            THICKNESS,
            VERTEX_OFFSET
        };

    public:
        WireInflatorFacade(std::string wire_file) :
            m_rows(0), m_cols(0) {
                m_inflator = WireInflator2D::construct(wire_file);
                m_t_params.max_area = 0.001;
        }

        virtual ~WireInflatorFacade();

        void set_dimension(size_t rows, size_t cols);

        size_t get_num_parameters() {
            return m_inflator->numberOfParameters();
        }

        ParameterType get_parameter_type(size_t i);
        VectorI get_affected_vertex_orbit(size_t i);
        MatrixF get_offset_direction(size_t i);
        void set_parameter(size_t row, size_t col, const VectorF& param);

        void set_max_triangle_area(Float max_area) {
            m_t_params.max_area = max_area;
        }

        void generate_periodic_pattern() {
            assert((*m_p_params)(0, 0) != NULL);
            m_inflator->generatePattern(*(*m_p_params)(0, 0), m_t_params, m_mesh);
        }

        void generate_tiled_pattern() {
            m_inflator->generateTiledPattern(*m_p_params, m_t_params, m_mesh);
        }

        void generate_pattern_with_guide_mesh(
                const VectorF& vertices, const VectorI& faces,
                const MatrixFr& raw_parameters);

        VectorF get_vertices();
        VectorI get_triangles();
        MatrixFr get_boundary_velocity();

    private:
        typedef Array2D<CellParameters*> ParameterGrid;
        typedef std::vector<CellParameters> ParameterVector;
        typedef std::shared_ptr<ParameterGrid> ParameterGridPtr;
        size_t m_rows;
        size_t m_cols;
        TessellationParameters m_t_params;
        ParameterGridPtr  m_p_params;
        WireInflator2D::OutMeshType m_mesh;
        WireInflator2D::Ptr m_inflator;
};
