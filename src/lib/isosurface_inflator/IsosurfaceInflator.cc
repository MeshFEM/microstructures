#include "IsosurfaceInflator.hh"

#include <memory>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <map>

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/SimplicialMesh.hh>
#include <MeshFEM/Future.hh>
#include <MeshFEM/Utilities/apply.hh>

#include "WireMesh.hh"
#include "PatternSignedDistance.hh"
#include "Isometries.hh"
#include "ShapeVelocityInterpolator.hh"

#include <MeshFEM/LinearElasticity.hh>

#include "IsosurfaceInflatorImpl.hh"

#define DEBUG_EVALPTS 0

// We can get faster builds for debugging/experimenting with the signed
// distance function if we disable the autodiff sections.
#ifndef SKIP_SVEL
#define SKIP_SVEL 0
#if SKIP_SVEL
#warning "Skipping shape velocity/normal computation!"
#endif
#endif

#include "CGALClippedVolumeMesher.hh"
#include "BoxIntersectionMesher.hh"
#include "MidplaneMesher.hh"
#include "IGLSurfaceMesherMC.hh"

using namespace std;

void IsosurfaceInflator::inflate(const vector<Real> &params) { m_imp->inflate(params); }
void IsosurfaceInflator::dumpInflationGraph(const std::string &path, const std::vector<Real> &params) const { m_imp->dumpInflationGraph(path, params); }

void IsosurfaceInflator::rasterize(const vector<Real> &params, const std::string &resolution, const std::string &outPath) { m_imp->rasterize(params, resolution, outPath); }
size_t       IsosurfaceInflator::numParams()                const { return m_imp->numParams(); }
vector<Real> IsosurfaceInflator::defaultParameters(Real t)  const { return m_imp->defaultParameters(t); }
bool         IsosurfaceInflator::isThicknessParam(size_t p) const { return m_imp->isThicknessParam(p); }
bool         IsosurfaceInflator:: isPositionParam(size_t p) const { return m_imp->isPositionParam(p); }
bool         IsosurfaceInflator:: isBlendingParam(size_t p) const { return m_imp->isBlendingParam(p); }
int          IsosurfaceInflator:: whichBlendingPolyParam(size_t p) const { return m_imp->whichBlendingPolyParam(p); }

bool IsosurfaceInflator::hasOrthotropicSymmetry() const { return m_imp->hasOrthotropicSymmetry(); }

void IsosurfaceInflator::setGenerateFullPeriodCell(bool onoff) { m_imp->generateFullPeriodCell = onoff; }
void IsosurfaceInflator::setReflectiveInflator(bool onoff) { m_imp->reflectiveInflator = onoff; }
BaseCellType IsosurfaceInflator::baseCellType() const {
    if      (m_imp->_meshingFullPeriodCell()) return BaseCellType::TriplyPeriodic;
    else if (m_imp->_meshingOrthoCell())      return BaseCellType::Orthotropic;
    else throw std::runtime_error("Unknown/incompatible base cell type.");
}

const vector<MeshIO::IOVertex > &IsosurfaceInflator::vertices()              const { return m_imp->vertices; }
const vector<MeshIO::IOElement> &IsosurfaceInflator::elements()              const { return m_imp->elements; }
const vector<vector<Real>>      &IsosurfaceInflator::normalShapeVelocities() const { return m_imp->normalShapeVelocities; }
const vector<Real>              &IsosurfaceInflator::inflatedParams()        const { return m_imp->inflatedParams; }

void                             IsosurfaceInflator::clear()                       { m_imp->clear(); }

auto IsosurfaceInflator::vertexNormals() const -> const vector<Point> & { return m_imp->vertexNormals; }

IsosurfaceInflator::Point
IsosurfaceInflator::trackSignedDistanceGradient(const IsosurfaceInflator::Point &pt) const {
    return m_imp->trackSignedDistanceGradient(pt);
}

void IsosurfaceInflator::disablePostprocess() { m_imp->m_noPostprocess = true; }
void IsosurfaceInflator:: enablePostprocess() { m_imp->m_noPostprocess = false; }

void IsosurfaceInflator::disableCheapPostprocess() { m_imp->m_cheapPostProcessing = false; }
void IsosurfaceInflator:: enableCheapPostprocess() { m_imp->m_cheapPostProcessing = true; }

// Printability checks/constraints
bool IsosurfaceInflator::isPrintable(const std::vector<Real> &params) const { return m_imp->isPrintable(params); }
Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>
IsosurfaceInflator::selfSupportingConstraints(const std::vector<Real> &params) const {
    return m_imp->selfSupportingConstraints(params);
}

MeshingOptions &IsosurfaceInflator::meshingOptions() { return m_imp->meshingOptions(); }

BBox<Point3<Real>> IsosurfaceInflator::meshingCell() const { return m_imp->meshingCell(); }

IsosurfaceInflator::~IsosurfaceInflator() {
    delete m_imp;
}

////////////////////////////////////////////////////////////////////////////////
// Factory/instantiations
////////////////////////////////////////////////////////////////////////////////
IsosurfaceInflator::IsosurfaceInflator(const string &type, bool vertexThickness, const string &wireMeshPath, size_t inflationNeighborhoodEdgeDist, size_t blendingPolySize) {
    if (!vertexThickness) throw runtime_error("Only per-vertex thickness is currently supported.");
    string name = type;
    transform(name.begin(), name.end(), name.begin(), ::tolower);

    // Decode symmetry type and mesher from the inflator name, which is
    // composed of three parts:
    // 1) a possible "2D_" prefix indicating a 2D inflator; this will set the mesher to MidplaneMesher
    // 2) a symmetry type
    // 3) a possible suffix selecting a custom mesher (only valid for 3D inflators):
    //      "_preview"  (mesh with marching cubes instead of CGAL)
    //      "_features" (mesh only the sharp feature curves that will be passed to CGAL)
    std::unique_ptr<MesherBase> mesher;
    size_t pos;
    bool shouldDisablePostprocess = false;
    if (name.find("2d_") == 0) {
        mesher = Future::make_unique<MidplaneMesher>();
        if (name.find("non_periodic") == std::string::npos) {
            name = name.substr(3, string::npos);
        }
    }
    else if ((pos = name.find("_preview")) != string::npos) {
        mesher = Future::make_unique<IGLSurfaceMesherMC>();
        name = name.substr(0, pos);
        shouldDisablePostprocess = true;
    }
    else if ((pos = name.find("_features")) != string::npos) {
        mesher = Future::make_unique<BoxIntersectionMesher>();
        name = name.substr(0, pos);
        shouldDisablePostprocess = true;
    }

    // The default mesher for 3D is CGAL.
    if (!mesher) mesher = Future::make_unique<CGALClippedVolumeMesher>();

    map<string, function<void()>> makeImplForSymmetry = {
        {"cubic",           [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::Cubic<>         >>(wireMeshPath, std::move(mesher), inflationNeighborhoodEdgeDist, blendingPolySize); }},
        {"orthotropic",     [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::Orthotropic<>   >>(wireMeshPath, std::move(mesher), inflationNeighborhoodEdgeDist, blendingPolySize); }},
        {"diagonal",        [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::Diagonal<>      >>(wireMeshPath, std::move(mesher), inflationNeighborhoodEdgeDist, blendingPolySize); }},
        {"2d_non_periodic", [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::NonPeriodic<DEFAULT_TOL, 2>>>(wireMeshPath, std::move(mesher), 0, blendingPolySize); }},
        {"non_periodic",    [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::NonPeriodic<DEFAULT_TOL, 3>>>(wireMeshPath, std::move(mesher), 0, blendingPolySize); }},
        {"triply_periodic", [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::TriplyPeriodic<>>>(wireMeshPath, std::move(mesher), inflationNeighborhoodEdgeDist, blendingPolySize); }},
        {"doubly_periodic", [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::DoublyPeriodic<>>>(wireMeshPath, std::move(mesher), inflationNeighborhoodEdgeDist, blendingPolySize); }},
        {"square",          [&]() { m_imp = new IsosurfaceInflatorImpl<WireMesh<Symmetry::Square<>        >>(wireMeshPath, std::move(mesher), inflationNeighborhoodEdgeDist, blendingPolySize); }}
    };

    if (makeImplForSymmetry.count(name) == 0) throw std::runtime_error("Invalid inflator name: '" + type + "'");
    makeImplForSymmetry[name]();

    if (shouldDisablePostprocess) disablePostprocess();
}
