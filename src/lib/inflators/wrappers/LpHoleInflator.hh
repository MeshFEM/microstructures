////////////////////////////////////////////////////////////////////////////////
// LpHoleInflator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Simple 2D inflator for debugging: an p-norm hole inflator.
//      There are two parameters: the radius, r in (0, 1), and the Lp-norm
//      parameter, p in (0, inf). The base cell is the square [-1, 1]^2 with a
//      hole of radius r and centered at (0, 0) subtracted.
//
//      To ensure a periodic boundary with minimal effort, the base cell is then
//      reflected into a 4x4 grid of cells.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/10/2015 22:09:49
////////////////////////////////////////////////////////////////////////////////
#ifndef LPHOLEINFLATOR_HH
#define LPHOLEINFLATOR_HH

#include "../Inflator.hh"
#include <MeshFEM/Triangulate.h>
#include <MeshFEM/filters/reflect.hh>
#include <MeshFEM/EdgeFields.hh>
#include <cmath>

struct LpHoleInflator : Inflator<2> {
    LpHoleInflator() { }

    ////////////////////////////////////////////////////////////////////////////
    // Inflation
    ////////////////////////////////////////////////////////////////////////////
private:
    virtual void m_inflate(const std::vector<Real> &params) override {
        assert(params.size() == numParameters());
        m_radius = params[0];
        m_p = params[1];

        // Create the square
        std::vector<MeshIO::IOVertex> inVertices = { {-1, -1, 0},
                                                     { 1, -1, 0},
                                                     { 1,  1, 0},
                                                     {-1,  1, 0} };
        std::vector<std::pair<size_t, size_t>> inEdges = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };

        size_t firstHoleVertex = inVertices.size();
        inVertices.emplace_back(m_radius, 0);

        double degreesPerSubdiv = 2 * M_PI / m_nsubdiv;
        // Create all hole boundary segments except the last
        for (size_t i = 1; i < m_nsubdiv; ++i) {
            double theta = degreesPerSubdiv * i;
            // https://www.mathworks.com/matlabcentral/newsreader/view_thread/279050
            inVertices.emplace_back(
                m_radius * copysign(pow(fabs(cos(theta)), 2 / m_p), cos(theta)),
                m_radius * copysign(pow(fabs(sin(theta)), 2 / m_p), sin(theta)));
            inEdges.push_back({inVertices.size() - 2,
                               inVertices.size() - 1});
        }

        // Close the hole path.
        inEdges.push_back({inVertices.size() - 1, firstHoleVertex});

        // Pick a point in the hole: the first segment forms a triangle with the
        // origin that lies entirely within the hole. Choose its barycenter.
        std::vector<MeshIO::IOVertex> holes;
        holes.emplace_back(((1 / 3.0) * (inVertices.at(firstHoleVertex    ).point +
                                         inVertices.at(firstHoleVertex + 1).point)).eval());

        triangulatePSLG(inVertices, inEdges, holes, m_vertices, m_elements,
                        m_triangleArea, "Q");

        reflect(2, m_vertices, m_elements, m_vertices, m_elements);
    }
public:

    ////////////////////////////////////////////////////////////////////////////
    // Shape velocity computation
    ////////////////////////////////////////////////////////////////////////////
    // Compute velocity induced on each vertex by changing radius and p.
    // Vertices are assumed to move only in the normal direction (since
    // tangential velocity is left undefined by levelset equations).
    //      phi(x, r, p) = 0 = r^p - |x|^p - |y|^p
    //  ==> grad phi . dx/dp + dphi/dp = 0,           n = grad phi / |grad phi|
    //  ==> |grad phi| n . dx/dp = -dphi/dp
    //  ==> n . dx/dp = -1/|grad phi| dphi/dp
    //  ==> dx/dp = -n/|grad phi| dphi/dp    (zero tangential motion)
    //  ==> dx/dp = - grad phi/|grad phi|^2 dphi/dp
    // WARNING: computes nonsense values on internal vertices (OK since these
    // are ignored by Inflator<2>::shapeVelocities())
    virtual std::vector<VectorField<Real, 2>> volumeShapeVelocities() const override {
        const size_t nv = this->m_vertices.size();
        std::vector<VectorField<Real, 2>> result(2, VectorField<Real, 2>(nv));
        // Velocity wrt radius (param 0) and p (param 1)
        VectorField<Real, 2> &vvel_r = result[0], &vvel_p = result[1];

        for (size_t vi = 0; vi < nv; ++vi) {
            std::array<Real, 2> coords, flipper;

            // Reflect back to [-1, 1]^2 base cell
            for (size_t c = 0; c < 2; ++c) {
                coords[c] = this->m_vertices[vi][c];
                if (coords[c] < -1) { flipper[c] = -1.0; coords[c] += 2; coords[c] *= -1; }
                else                { flipper[c] =  1.0; }
                assert((coords[c] >= -1) && (coords[c] <= 1.0));
            }

            Real absX = fabs(coords[0]),
                 absY = fabs(coords[1]);

            // Outward normal components (unnormalized)
            Real phi_x = -copysign(m_p * pow(absX, m_p - 1.0), coords[0]);
            Real phi_y = -copysign(m_p * pow(absY, m_p - 1.0), coords[1]);

            Real phi_r = m_p * pow(m_radius, m_p - 1.0);
            Real phi_p = pow(m_radius, m_p) * log(m_radius)
                       - ((absX > 1e-9) ? pow(absX, m_p) * log(absX) : 0.0)
                       - ((absY > 1e-9) ? pow(absY, m_p) * log(absY) : 0.0);

            Real inv_grad_norm_sq = 1.0 / (phi_x * phi_x + phi_y * phi_y);
            vvel_r(vi)[0] = -flipper[0] * inv_grad_norm_sq * phi_x * phi_r;
            vvel_r(vi)[1] = -flipper[1] * inv_grad_norm_sq * phi_y * phi_r;

            vvel_p(vi)[0] = -flipper[0] * inv_grad_norm_sq * phi_x * phi_p;
            vvel_p(vi)[1] = -flipper[1] * inv_grad_norm_sq * phi_y * phi_p;
        }
        
        return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Queries
    ////////////////////////////////////////////////////////////////////////////
    virtual bool isParametric() const override { return true; }
    virtual size_t numParameters() const override { return 2; }
    // Parameters: (radius, p)
    ParameterType parameterType(size_t p) const override {
        if (p == 0) return ParameterType::Thickness;
        if (p == 1) return ParameterType::Blending; // Sorta...
        throw std::runtime_error("Invalid parameter: " + std::to_string(p));
    }

    // 2D is always printable.
    virtual bool isPrintable(const std::vector<Real> &/* params */) override { return true; }

    ////////////////////////////////////////////////////////////////////////////
    // Configuration
    ////////////////////////////////////////////////////////////////////////////
    virtual void setMaxElementVolume(Real maxElementVol) override { m_triangleArea = maxElementVol; }
    void setNumSubdiv(size_t ns) { m_nsubdiv = ns; }
    virtual void setReflectiveInflator(bool use) override { if (!use) throw std::runtime_error("LpHoleInflator is always reflective."); }

    ////////////////////////////////////////////////////////////////////////////
    // Analytic normals (not part of Inflator interface)
    ////////////////////////////////////////////////////////////////////////////
    // Isosurface normal:
    // phi(x, p) = 0 = r^p - x^p - y^p
    // grad phi / |grad phi|
    template<class _FEMMesh>
    VectorField<Real, 2> analyticNormals(const _FEMMesh &mesh) const {
        VectorField<Real, 2> normals;
        normals.resizeDomain(mesh.numBoundaryVertices());

        for (auto bv : mesh.boundaryVertices()) {
            auto v = bv.volumeVertex().node()->p;
            Vector2D n;
            for (size_t c = 0; c < 2; ++c) {
                bool flip = false;
                // map back to base cell
                if (v[c] < -1) { v[c] += 2; v[c] *= -1; flip = true; }
                n[c] = -copysign(m_p * pow(fabs(v[c]), m_p - 1.0), v[c]); // grad_c phi
                if (flip) n[c] *= -1;
            }
            normals(bv.index()) = n.normalized();
        }

        // Clear normals on the periodic boundary
        for (auto be : mesh.boundaryElements()) {
            if (be->isInternal) {
                normals(be.vertex(0).index()).setZero();
                normals(be.vertex(1).index()).setZero();
            }
        }

        return normals;
    }

    virtual ~LpHoleInflator() { }

private:
    size_t m_nsubdiv = 64;
    Real   m_triangleArea = 0.001;
    Real   m_radius, m_p; // Set by inflation
};

#endif /* end of include guard: LPHOLEINFLATOR_HH */
