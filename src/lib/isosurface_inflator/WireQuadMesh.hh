#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#pragma once
////////////////////////////////////////////////////////////////////////////////
#include "BilinearMap.hh"
#include "Symmetry.hh"
#include "WireMesh.hh"
#include <json.hpp>
#include <stdexcept>
#include <memory>
#include <set>
////////////////////////////////////////////////////////////////////////////////

//
// By default, if no active quad is chosen (index < 0), then we should mesh
// the entire graph!!
//
class WireQuadMesh {
public:
    using PatternSymmetry = Symmetry::Null<>;
    using WireMeshBasePtr = std::shared_ptr<WireMeshBase>;

    using Point = WireMeshBase::Point; // Point3<double>;
    using Edge  = WireMeshBase::Edge; // std::pair<size_t, size_t>

    struct MapToBaseUnit {
        BilinearMap func_;
        MapToBaseUnit(BilinearMap f = BilinearMap()) : func_(f) { }

        template<typename Real>
        Point3<Real> operator() (Point3<Real> p) const {
            Point3<Real> q = func_.apply(p[0], p[1]);
            q[2] = p[2];
            return q;
        }
    };

public:
    WireQuadMesh(
        const std::vector<MeshIO::IOVertex> &V,
        const std::vector<MeshIO::IOElement> &F,
        const nlohmann::json &params);

    ThicknessType thicknessType() const { return m_thicknessType; }

    // Set currently active quad
    void setActiveQuad(int idx);

    // Return wiremesh associated to the currently active quad
    const WireMeshBase &activeWireMesh() const;

    // Inflation parameters the whole graph (simple concatenation)
    std::vector<double> params() const;

    // Area of each quad of the background mesh
    Eigen::VectorXd areas() const;

    MapToBaseUnit mapFunctor() const { return MapToBaseUnit(m_bilinearMap); }

    BBox<Point3d> boundingBox() const { return m_bbox; }

    BilinearMap getBilinearMap(int i) const;

    double getScalingFactor(int i) const;

    // Build the inflation graph for the whole quad mesh, stitching together adjacent nodes
    // (averaging stitched points' locations, thicknesses, and blending params).
    //
    // @param[in]  allParams               { Inflation parameters for the whole mesh }
    // @param[out] stitchedPoints          { Graph vertices positions }
    // @param[out] stitchedEdges           { Graph edges indices }
    // @param[out] stitchedThicknesses     { Graph edge thicknesses }
    // @param[out] stitchedBlendingParams  { Graph vertices blending parameters }
    //
    void inflationGraph(const std::vector<double> &allParams,
        std::vector<Point>  &stitchedPoints,
        std::vector<Edge>   &stitchedEdges,
        std::vector<double> &stitchedThicknesses,
        std::vector<double> &stitchedBlendingParams,
        std::vector<std::vector<double>> &stitchedBlendingPolyParams) const;

private:
    // Fill m_allJacobians based on the background quad mesh (m_V, m_F)
    void compute_jacobians();

private:
    Eigen::MatrixXd m_V;
    Eigen::MatrixXi m_F;

    int m_activeQuad = -1;
    ThicknessType m_thicknessType = ThicknessType::Vertex;

    BilinearMap m_bilinearMap;
    BBox<Point3d> m_bbox;

    std::vector<WireMeshBasePtr> m_allTopologies;
    std::vector<std::vector<double>> m_allParameters;
    std::vector<Eigen::Matrix3d> m_allJacobians; // ref square [-1,1]² to mapped parallelogram
};

#endif
