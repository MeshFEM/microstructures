#include "IsoinflatorWrapper.hh"

#include <isosurface_inflator/IsosurfaceInflator.hh>
#include <isosurface_inflator/IsosurfaceInflatorConfig.hh>

#include <boost/algorithm/string.hpp>
#include <MeshFEM/Future.hh>
#include <MeshFEM/MSHFieldWriter.hh>

#include <iostream>
#include <iomanip>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
// Override geometry access
////////////////////////////////////////////////////////////////////////////////
template<size_t N> const std::vector<MeshIO::IOElement> &IsoinflatorWrapper<N>::elements() const { return m_inflator->elements(); }
template<size_t N> const std::vector<MeshIO::IOVertex>  &IsoinflatorWrapper<N>::vertices() const { return m_inflator->vertices(); }
template<size_t N> void IsoinflatorWrapper<N>::clear() { m_inflator->clear(); }

// cheap post processing can be used in isosurface_cli, skipping computation of normals and shape velocities
template<size_t N> void IsoinflatorWrapper<N>::disableCheapPostprocess() { m_inflator->disableCheapPostprocess();}
template<size_t N> void IsoinflatorWrapper<N>::enableCheapPostprocess() { m_inflator->enableCheapPostprocess();}


////////////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
IsoinflatorWrapper<N>::IsoinflatorWrapper(const std::string &wireMeshPath,
                                          const std::string &symmetryType, bool vertex_thickness,
                                          size_t inflationGraphRadius) {
    std::string inflatorName(symmetryType);
    boost::algorithm::to_lower(inflatorName);
    inflatorName = ((N == 2) ? "2D_" : "") + inflatorName;
    m_inflator = Future::make_unique<IsosurfaceInflator>(
                inflatorName,
                vertex_thickness, wireMeshPath, inflationGraphRadius);
}

////////////////////////////////////////////////////////////////////////////////
// Inflation
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
void IsoinflatorWrapper<N>::m_inflate(const std::vector<Real> &params)
{
    m_inflator->inflate(params);
}

////////////////////////////////////////////////////////////////////////////////
// Shape velocity computation
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
std::vector<VectorField<Real, N>>
IsoinflatorWrapper<N>::volumeShapeVelocities() const
{
    std::vector<std::vector<Real>> nsv = m_inflator->normalShapeVelocities();
    std::vector<Point3D>             n = m_inflator->vertexNormals();

    const size_t np = numParameters();
    const size_t nv = this->vertices().size();
    assert(nv == n.size());
    assert(nv == nsv.at(0).size());

    std::vector<VectorField<Real, N>> result(np);
    try {
        double maxZMag = 0;
        for (size_t p = 0; p < np; ++p) {
            result[p].resizeDomain(nv);
            for (size_t vi = 0; vi < nv; ++vi) {
                maxZMag = std::max(maxZMag, std::abs(n[vi][2]));
                // With the sphere convex hull blending region, the signed
                // distance function is often spatially nondifferentiable at
                // the midplane in the z direction (the signed distance function
                // is only C0 on the medial axis). This causes us to get slighly
                // incorrect normal components in the z direction.
                // This cannot be fixed easily in 3D, but the error is hopefully
                // small.
                if (N == 2) {
                    n[vi][2] = 0;
                    n[vi] /= n[vi].norm();
                }

                result[p](vi)  = truncateFrom3D<VectorND<N>>(n[vi]);
                result[p](vi) *= nsv[p][vi];
                //assert(!std::isnan(result[p](vi)));
                //assert(!std::isnan(result[p](vi)));
            }
        }

        if ((N == 2) && (maxZMag > 1e-1)) {
            std::cerr << "Large normal z component: "
                      << maxZMag << " (probably ok)"
                      << std::endl;
        }
        if ((N == 2) && (maxZMag > 0.4)) {
            std::cerr << "Extremely large z component in 2D normal: " << maxZMag
                      << "; throwing error" << std::endl;
            throw std::runtime_error("Extremely large z component in 2D normal.");
        }

        //MSHFieldWriter writer("debug_velocities.msh", vertices(), elements());
        //size_t i = 0;
        //for (const auto &vv : result) {
        //    writer.addField("vv " + std::to_string(i++), vv, DomainType::PER_NODE);
        //}
    }
    catch (...) {
        std::cerr.precision(19);
        std::cerr << "Error for pattern parameters: " << std::endl;
        for (Real pval : m_inflator->inflatedParams())
            std::cerr << pval << "\t";
        std::cerr << std::endl;

        MSHFieldWriter writer("debug_velocities.msh", vertices(), elements());
        {
            VectorField<double, 3> n_out(n.size());
            ScalarField<double> normal_z(n.size());
            for (size_t i = 0; i < n.size(); ++i) {
                n_out(i) = n[i];
                normal_z[i] = n[i][2];
            }
            writer.addField("normals", n_out);
            writer.addField("normal_zcomp", normal_z);
        }
        size_t i = 0;
        for (const auto &nsv_p : nsv) {
            ScalarField<double> nsv_out(nsv_p);
            writer.addField("nsv " + std::to_string(i++), nsv_out);
        }
        throw;
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Queries
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
size_t
IsoinflatorWrapper<N>::numParameters() const { return m_inflator->numParams(); }

template<size_t N>
ParameterType
IsoinflatorWrapper<N>::parameterType(size_t p) const
{
    if (m_inflator->isThicknessParam(p)) return ParameterType::Thickness;
    if (m_inflator-> isPositionParam(p)) return ParameterType::Offset;
    if (m_inflator-> isBlendingParam(p)) return ParameterType::Blending;
    throw std::runtime_error("Unknown parameter type for param " + std::to_string(p));
}

template<size_t N>
std::vector<Real>
IsoinflatorWrapper<N>::defaultParameters() const {
    return m_inflator->defaultParameters();
}

template<size_t N>
bool IsoinflatorWrapper<N>::isPrintable(const std::vector<Real> &params) {
    if (N == 2) return true; // 2D patterns are always printable.
    return m_inflator->isPrintable(params); }

template<size_t N>
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
IsoinflatorWrapper<N>::selfSupportingConstraints(const std::vector<double> &params) const {
    return m_inflator->selfSupportingConstraints(params);
}
template<size_t N>
bool
IsoinflatorWrapper<N>::hasOrthotropicSymmetry() const {
    return m_inflator->hasOrthotropicSymmetry();
}

template<size_t N>
BBox<Vector3D>
IsoinflatorWrapper<N>::meshingCell() {
    return m_inflator->meshingCell();
}

////////////////////////////////////////////////////////////////////////////////
// Configuration
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
void IsoinflatorWrapper<N>::loadMeshingOptions(const std::string &moptsPath) { m_inflator->meshingOptions().load(moptsPath); }

template<size_t N>
MeshingOptions &IsoinflatorWrapper<N>::meshingOptions() { return m_inflator->meshingOptions(); }

template<size_t N>
void IsoinflatorWrapper<N>::setMaxElementVolume(Real maxElementVol)
{
    if (N == 2) m_inflator->meshingOptions().maxArea = maxElementVol;
    else        m_inflator->meshingOptions().cellSize = maxElementVol;
}

template<size_t N> void IsoinflatorWrapper<N>::setOrthoBaseCell(bool ortho) { m_inflator->setGenerateFullPeriodCell(!ortho); }
template<size_t N> BaseCellType IsoinflatorWrapper<N>::baseCellType() const { return m_inflator->baseCellType(); }

////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
IsoinflatorWrapper<N>::~IsoinflatorWrapper() { }

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: 2D and 3D inflators.
////////////////////////////////////////////////////////////////////////////////
template class IsoinflatorWrapper<2>;
template class IsoinflatorWrapper<3>;
