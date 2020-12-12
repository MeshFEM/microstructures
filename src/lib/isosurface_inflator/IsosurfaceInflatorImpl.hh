////////////////////////////////////////////////////////////////////////////////
// IsosurfaceInflatorImpl.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Helper classes for IsosurfaceInflator that do all the actual implementation.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  08/23/2017 17:24:39
////////////////////////////////////////////////////////////////////////////////
#ifndef ISOSURFACEINFLATORIMPL_HH
#define ISOSURFACEINFLATORIMPL_HH

#include "AutomaticDifferentiation.hh"
#include "IsosurfaceInflator.hh"
#include "IsosurfaceInflatorConfig.hh"
#include "MesherBase.hh"
#include "PatternSignedDistance.hh"
#include "PostProcess.hh"
#include "rasterize.hh"

#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/StringUtils.hh>

#include <type_traits>


class IsosurfaceInflator::Impl {
public:
    virtual std::vector<Real> defaultParameters(Real thickness = 0.07) const = 0;
    virtual size_t                numParams() const = 0;
    virtual bool   isThicknessParam(size_t p) const = 0;
    virtual bool    isPositionParam(size_t p) const = 0;
    virtual bool    isBlendingParam(size_t p) const = 0;
    virtual int     whichBlendingPolyParam(size_t p) const = 0;

    virtual       MeshingOptions &meshingOptions()       = 0;
    virtual const MeshingOptions &meshingOptions() const = 0;

    // Delegates to derived IsosurfaceInflatorImpl for WMesh-dependent stuff
    // (via the virtual functions below).
    void inflate(const std::vector<Real> &params) {
        inflatedParams = params;

        if (meshingOptions().debugLoadMeshPath.size()) {
            MeshIO::load(meshingOptions().debugLoadMeshPath, vertices, elements);
            // The normal/signed computation expects parameters to have been set
            m_setParameters(params);
        }
        else {
            meshPattern(params);
        }

        // Terminate after initial meshing, for debugging.
        if (m_noPostprocess) return;

        // std::cout << "Meshed params:";
        // for (Real p : params) std::cout << "\t" << p;
        // std::cout << std::endl;

        // Determine if meshed domain is 2D or 3D and postprocess accordingly
        BBox<Point> bbox(vertices);
        if (std::abs(bbox.dimensions()[2]) < 1e-8) postProcess<2>(vertices, elements, normalShapeVelocities, vertexNormals, *this, !_meshingOrthoCell(), generateFullPeriodCell, meshingCell(), meshingOptions(), m_cheapPostProcessing, m_nonPeriodicity);
        else                                       postProcess<3>(vertices, elements, normalShapeVelocities, vertexNormals, *this, !_meshingOrthoCell(), generateFullPeriodCell, meshingCell(), meshingOptions(), m_cheapPostProcessing, m_nonPeriodicity);

        if (meshingOptions().debugSVelPath.size()) {
            MSHFieldWriter writer(meshingOptions().debugSVelPath, vertices, elements);
            VectorField<Real, 3> normals(vertices.size());
            for (size_t i = 0; i < vertices.size(); ++i)
                normals(i) = vertexNormals[i];
            writer.addField("normals", normals, DomainType::PER_NODE);

            ScalarField<Real> vn(vertices.size());
            for (size_t p = 0; p < normalShapeVelocities.size(); ++p) {
                for (size_t i = 0; i < vertices.size(); ++i)
                    vn(i) = normalShapeVelocities[p][i];
                writer.addField("svel " + std::to_string(p), vn, DomainType::PER_NODE);
            }

#if 0
            // Get interpolated velocities
            using Mesh = LinearElasticity::Mesh<2, 1, LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter>;
            LinearElasticity::Simulator<Mesh> sim(elements, vertices);
            const auto &mesh = sim.mesh();
            ShapeVelocityInterpolator interpolator(sim);
            VectorField<Real, 2> bvel(mesh.numBoundaryVertices());
            for (size_t p = 0; p < normalShapeVelocities.size(); ++p) {
                for (auto bv : mesh.boundaryVertices()) {
                    size_t vi = bv.volumeVertex().index();
                    bvel(bv.index()) = truncateFrom3D<Point2D>(vertexNormals.at(vi)) * normalShapeVelocities[p][vi];
                }
                auto vvel = interpolator.interpolate(sim, bvel);
                writer.addField("vvel " + std::to_string(p), vvel, DomainType::PER_NODE);
            }
#endif
        }
    }

    // Mesh the param (fills vertices, elements member arrays)
    virtual void meshPattern(const std::vector<Real> &params) = 0;

    // Rasterize to a density field on a 2D/3D grid
    // (Infer dimension from resolutionString, which specifies rasterization grid size along each dimension)
    virtual void rasterize(const std::vector<Real> &params, const std::string &resolutionString, const std::string &outPath) = 0;

    // Dump inflation graph
    virtual void dumpInflationGraph(const std::string &path, const std::vector<Real> &params) const = 0;

    // Derivative of signed distance function with respect to evaluation point
    virtual std::vector<Point> signedDistanceGradient(const std::vector<Point> &evalPoints) const = 0;

    // Derivative of signed distance function wrt each pattern parameter
    virtual std::vector<std::vector<Real>> signedDistanceParamPartials(const std::vector<Point> &evalPoints) const = 0;

    // Determine the degree of smoothing at each evaluation point (for debugging purposes)
    virtual std::pair<std::vector<Real>,std::vector<size_t>> smoothnessAtPoints(const std::vector<Point> &evalPoints) const = 0;

    // Debug the signed distance gradient, e.g. to diagnose nonzero z components
    // of 2D normals
    virtual Point trackSignedDistanceGradient(const Point &pt) const = 0;

    // Printability checks/constraints
    virtual bool isPrintable(const std::vector<Real> &params) const = 0;
    virtual Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
        selfSupportingConstraints(const std::vector<Real> &params) const = 0;

    // Whether the inflator generates only the orthotropic symmetry base cell
    // (by default--reflectiveInflator will override this)
    virtual bool _meshingOrthoCell() const = 0;

    // Whether the inflator generates the full doubly- or triply-periodic mesh.
    virtual bool _meshingFullPeriodCell() const = 0;

    virtual bool hasOrthotropicSymmetry() const = 0;

    // The meshing cell box.
    virtual BBox<Point> meshingCell() const = 0;

    void clear() {
        vertices.clear(), elements.clear();
        normalShapeVelocities.clear();
        vertexNormals.clear();
        inflatedParams.clear();
    }

    virtual ~Impl() { }

    ////////////////////////////////////////////////////////////////////////////
    // "Private" data memebers
    ////////////////////////////////////////////////////////////////////////////
    // Chooses whether the full period cell is created in cases where the
    // representative cell is smaller (e.g. orthotropic patterns, whose base
    // cell is 1/8th the period cell).
    bool generateFullPeriodCell = true;
    // Chooses whether we mesh only the symmetry base cell and reflect to the
    // full period cell if necessary (default), or we mesh the full period cell
    // when generateFullPeriodCell == true
    bool reflectiveInflator = true;
    bool m_noPostprocess = false; // for debugging
    bool m_cheapPostProcessing = false; // only valid when using inflation alone
    std::vector<MeshIO::IOVertex>  vertices;
    std::vector<MeshIO::IOElement> elements;
    std::vector<std::vector<Real>> normalShapeVelocities;
    std::vector<Point>             vertexNormals;

    // The params that were most recently inflated (to which
    // (vertices, elements) correspond).
    std::vector<Real> inflatedParams;

    bool m_nonPeriodicity = false; // use together with NonPeriodic symmetry
protected:
    // Manually set the parameters in PatternSignedDistance instance
    // (useful if bypassing meshing for debugging).
    virtual void m_setParameters(const std::vector<Real> &params) = 0;
};

// The WMesh-dependent implementation details.
// E.g.: WMesh = WireMesh<Symmetry::Cubic<>>
template<class WMesh>
class IsosurfaceInflatorImpl : public IsosurfaceInflator::Impl {
public:
    using Point = IsosurfaceInflator::Point;
    typedef PatternSignedDistance<Real, WMesh> PSD;
    typedef typename WMesh::PatternSymmetry PatternSymmetry;
    IsosurfaceInflatorImpl(const std::string &wireMeshPath, std::unique_ptr<MesherBase> &&m, size_t inflationNeighborhoodEdgeDist, size_t blendingPolySize)
        : wmesh(wireMeshPath, inflationNeighborhoodEdgeDist, blendingPolySize), pattern(wmesh), mesher(std::move(m)) {

        m_nonPeriodicity = std::is_same<PatternSymmetry, Symmetry::NonPeriodic<typename PatternSymmetry::Tolerance, 2>>::value ||
                           std::is_same<PatternSymmetry, Symmetry::NonPeriodic<typename PatternSymmetry::Tolerance, 3>>::value;
    }

    virtual void meshPattern(const std::vector<Real> &params) override {
#if 0
        std::cerr << "Meshing parameters:";
        for (auto p : params)
            std::cout << "\t" << p;
        std::cerr << std::endl;
#endif

        BENCHMARK_START_TIMER_SECTION("meshPattern");
        // Optional debugging graph output.
        const auto &config = IsosurfaceInflatorConfig::get();
        if (config.dumpInflationGraph())  { wmesh.saveInflationGraph(    config.inflationGraphPath, params); }
        if (config.dumpReplicatedGraph()) { wmesh.saveReplicatedBaseUnit(config.replicatedGraphPath); }
        if (config.dumpBaseUnitGraph())   { wmesh.saveBaseUnit(          config.baseUnitGraphPath); }

        pattern.setParameters(params, meshingOptions().jacobian, meshingOptions().jointBlendingMode, meshingOptions().jointBlendingFunction);

        // Change the pattern's meshing domain if we're forcing meshing of the
        // full TriplyPeriodic base cell.
        if (generateFullPeriodCell && !reflectiveInflator)
            pattern.setBoundingBox(Symmetry::TriplyPeriodic<>::representativeMeshCell<Real>());

        mesher->meshInterfaceConsistently = !_meshingOrthoCell() ||
                                            meshingOptions().forceConsistentInterfaceMesh;

        mesher->mesh(pattern, this->vertices, this->elements);
        BENCHMARK_STOP_TIMER_SECTION("meshPattern");
        // cout << vertices.size() << " vertices, " << elements.size() << " elements" << endl;
    }

    // Rasterize to a indicator scalar field on a 2D/3D grid
    // (Infer dimension from resolutionString, which specifies rasterization grid size along each dimension)
    virtual void rasterize(const std::vector<Real> &params, const std::string &resolutionString, const std::string &outPath) override {
        pattern.setParameters(params, meshingOptions().jacobian, meshingOptions().jointBlendingMode, meshingOptions().jointBlendingFunction);

        std::vector<MeshIO::IOVertex > vertices;
        std::vector<MeshIO::IOElement> elements;
        ScalarField<Real> indicator;


        std::vector<std::string> sizeStrings = MeshFEM::split(resolutionString, ",x");
        std::vector<size_t> sizes;
        for (auto &s : sizeStrings) sizes.push_back(std::stoul(s));
        ::rasterize(pattern, sizes, vertices, elements, indicator);

        auto type = (sizes.size() == 2) ? MeshIO::MESH_QUAD : MeshIO::MESH_HEX;
        MSHFieldWriter writer(outPath, vertices, elements, type);
        writer.addField("indicator", indicator, DomainType::PER_ELEMENT);
    }

    // Dump inflation graph
    virtual void dumpInflationGraph(const std::string &path, const std::vector<Real> &params) const override {
        wmesh.saveInflationGraph(path, params);
    }

    // Note: when reflectiveInflator = false, the mesher generates the full
    // period cell.
    virtual bool _meshingOrthoCell() const override {
        return std::is_base_of<Symmetry::Orthotropic<typename PatternSymmetry::Tolerance>,
                               PatternSymmetry>::value
               && (reflectiveInflator || !generateFullPeriodCell);
    }

    // Note: when reflectiveInflator = false, the mesher generates the full
    // period cell.
    virtual bool _meshingFullPeriodCell() const override {
        return generateFullPeriodCell || !hasOrthotropicSymmetry();
    }

    virtual bool hasOrthotropicSymmetry() const override {
        return !(std::is_same<PatternSymmetry, Symmetry::Diagonal<typename PatternSymmetry::Tolerance>>::value
              || std::is_same<PatternSymmetry, Symmetry::TriplyPeriodic<typename PatternSymmetry::Tolerance>>::value
              || std::is_same<PatternSymmetry, Symmetry::DoublyPeriodic<typename PatternSymmetry::Tolerance>>::value);
    }

    // The meshing cell box.
    virtual BBox<Point> meshingCell() const override { return PatternSymmetry::template representativeMeshCell<Real>(); }

    // Derivative of signed distance function with respect to evaluation point
    // (autodiff-based).
    virtual std::vector<Point> signedDistanceGradient(const std::vector<Point> &evalPoints) const override {
        std::vector<Point> distGradX(evalPoints.size());
#if !SKIP_SVEL
        // Scalar supporting 3D spatial gradient
        using ADScalar = Eigen::AutoDiffScalar<Vector3<Real>>;
        using Pt = Point3<ADScalar>;
        const size_t nEvals = evalPoints.size();

        auto evalAtPtIdx = [&](size_t p) {
            Pt x(ADScalar(evalPoints[p][0], 3, 0),
                 ADScalar(evalPoints[p][1], 3, 1),
                 ADScalar(evalPoints[p][2], 3, 2));
            ADScalar dist = pattern.signedDistance(x);
            distGradX[p]  = dist.derivatives();
        };

#if MICRO_WITH_TBB
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, nEvals),
            [&](const tbb::blocked_range<size_t> &r) {
                for (size_t p = r.begin(); p < r.end(); ++p) evalAtPtIdx(p);
            });
#else
        for (size_t p = 0; p < nEvals; ++p) evalAtPtIdx(p);
#endif
#else   // SKIP_SVEL
        for (auto &p : distGradX)
            p = Point(1, 0, 0);
#endif // !SKIP_SVEL
        return distGradX;
    }

    virtual Point trackSignedDistanceGradient(const Point &pt) const override {
#if !SKIP_SVEL
        // Scalar supporting 3D spatial gradient
        using ADScalar = Eigen::AutoDiffScalar<Vector3<Real>>;
        using Pt = Point3<ADScalar>;
        Pt x(ADScalar(pt[0], 3, 0),
             ADScalar(pt[1], 3, 1),
             ADScalar(pt[2], 3, 2));
        ADScalar dist = pattern.template signedDistance<ADScalar, true>(x);
        return dist.derivatives();
#else
        return Point(1, 0, 0);
#endif
    }

    // Derivative of signed distance function with respect to each pattern
    // parameter (autodiff-based).
    virtual std::vector<std::vector<Real>> signedDistanceParamPartials(const std::vector<Point> &evalPoints) const override {
        const size_t nEvals = evalPoints.size(),
                    nParams = pattern.numParams();
        std::vector<std::vector<Real>> partials(nParams, std::vector<Real>(nEvals));
#if !SKIP_SVEL
        // Scalar supporting derivatives with respect to each pattern parameter
        using PVec = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
        using ADScalar = Eigen::AutoDiffScalar<PVec>;

        assert(inflatedParams.size() == nParams);

        PatternSignedDistance<ADScalar, WMesh> patternAutodiff(wmesh);
        std::vector<ADScalar> params;
        params.reserve(params.size());
        for (size_t p = 0; p < nParams; ++p)
            params.emplace_back(inflatedParams[p], nParams, p);
        patternAutodiff.setParameters(params, meshingOptions().jacobian, meshingOptions().jointBlendingMode, meshingOptions().jointBlendingFunction);

        auto evalAtPtIdx = [&](size_t e) {
            ADScalar sd = patternAutodiff.signedDistance(evalPoints[e].template cast<ADScalar>().eval());
            Real sdOrig = pattern.signedDistance(evalPoints[e]);
            if (std::abs(stripAutoDiff(sd) - sdOrig) > 1.25e-2) {
                throw std::runtime_error("Incorrect signed distance computed by autodiff version: differ by " + std::to_string(std::abs(stripAutoDiff(sd) - sdOrig)));
            }

            for (size_t p = 0; p < nParams; ++p) {
                partials[p][e] = sd.derivatives()[p];
                if (std::isnan(partials[p][e])) {
                    std::cerr << "nan sd partial " << p << " at evaluation point "
                              << evalPoints[e] << std::endl;
                    std::cerr << "sd at pt:\t" << pattern.signedDistance(evalPoints[e]) << std::endl;
                    std::cerr << "autodiff sd at pt:\t" << sd << std::endl;
                    // throw std::runtime_error("nan sd");
                }
            }
        };
#if MICRO_WITH_TBB
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, nEvals),
            [&](const tbb::blocked_range<size_t> &r) {
                for (size_t e = r.begin(); e < r.end(); ++e) evalAtPtIdx(e);
            });
#else
        for (size_t e = 0; e < nEvals; ++e) evalAtPtIdx(e);
#endif
#endif // !SKIP_SVEL
        return partials;
    }

    // Determine the degree of smoothing at each evaluation point (for debugging purposes)
    virtual std::pair<std::vector<Real>,std::vector<size_t>> smoothnessAtPoints(const std::vector<Point> &evalPoints) const override {
        std::vector<Real> smoothness;
        std::vector<size_t> closestVtx;
        smoothness.reserve(evalPoints.size()), closestVtx.reserve(evalPoints.size());
        for (const auto &pt : evalPoints) {
            Real s;
            size_t vtx;
            std::tie(s, vtx) = pattern.smoothnessAndClosestVtx(pt);
            smoothness.push_back(s);
            closestVtx.push_back(vtx);
        }
        return std::make_pair(smoothness, closestVtx);
    }

    virtual ~IsosurfaceInflatorImpl() { }

    virtual std::vector<Real> defaultParameters(Real t) const override { return wmesh.defaultParameters(t); }
    virtual size_t               numParams() const override { return wmesh.numParams(); }
    virtual bool  isThicknessParam(size_t p) const override { return wmesh.isThicknessParam(p); }
    virtual bool   isPositionParam(size_t p) const override { return wmesh.isPositionParam(p); }
    virtual bool   isBlendingParam(size_t p) const override { return wmesh.isBlendingParam(p); }

    // Answers which blending poly coefficient index p corresponds to. If index p does not correspond to blending poly
    // param, then returns -1.
    virtual int    whichBlendingPolyParam(size_t p) const override { return wmesh.whichBlendingPolyParam(p); }

    virtual       MeshingOptions &meshingOptions()       override { return mesher->meshingOptions; }
    virtual const MeshingOptions &meshingOptions() const override { return mesher->meshingOptions; }

    // Printability checks/constraints
    virtual bool isPrintable(const std::vector<Real> &params) const override {
        return wmesh.isPrintable(params);
    }

    virtual Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
        selfSupportingConstraints(const std::vector<Real> &params) const override {
        return wmesh.selfSupportingConstraints(params);
    }

    WMesh wmesh;
    PSD pattern;
    std::unique_ptr<MesherBase> mesher;
protected:
    virtual void m_setParameters(const std::vector<Real> &params) override {
        pattern.setParameters(params, meshingOptions().jacobian, meshingOptions().jointBlendingMode);
    }
};

#endif /* end of include guard: ISOSURFACEINFLATORIMPL_HH */
