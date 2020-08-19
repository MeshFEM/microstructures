#include "IGLSurfaceMesherMC.hh"

#if HAS_LIBIGL
#include <igl/copyleft/marching_cubes.h>

void IGLSurfaceMesherMC::
mesh(const SignedDistanceRegion<3> &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements,
     const double isolevel) const
{
    size_t gs = meshingOptions.marchingCubesGridSize;
    const auto bbox = sdf.boundingBox();

    const size_t nsamples = gs * gs * gs;

    // Evaluation point locations;
    // flattened to be accessed as:
    // xi + gs * (yi + gs * zi)
    Eigen::MatrixXd sampleLocations(nsamples, 3);
    {
        size_t i = 0;
        for (size_t zi = 0; zi < gs; ++zi) {
            for (size_t yi = 0; yi < gs; ++yi) {
                for (size_t xi = 0; xi < gs; ++xi) {
                    sampleLocations.row(i) = bbox.interpolatePoint(
                            Point3D(xi / Real(gs - 1.0),
                                    yi / Real(gs - 1.0),
                                    zi / Real(gs - 1.0)));
                    ++i;
                }
            }
        }
    }

    // Evaluate signed distances at each grid point
    Eigen::VectorXd signedDistances(nsamples);
#if MICRO_WITH_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples),
            [&](const tbb::blocked_range<size_t> &r) {
                for (size_t i = r.begin(); i < r.end(); ++i)
                    signedDistances(i) = sdf.signedDistance(sampleLocations.row(i)) - isolevel;
            }
        );
#else
    for (size_t i = 0; i < nsamples; ++i)
        signedDistances(i) = sdf.signedDistance(sampleLocations.row(i)) - isolevel;
#endif
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::copyleft::marching_cubes(signedDistances, sampleLocations,
            gs, gs, gs, V, F);

    vertices.clear(), elements.clear();
    vertices.reserve(V.rows());
    elements.reserve(F.rows());
    for (int i = 0; i < V.rows(); ++i)
        vertices.emplace_back(Point3D(V.row(i)));
    for (int i = 0; i < F.rows(); ++i)
        elements.emplace_back(F(i, 0), F(i, 1), F(i, 2));
}

#else // !HAS_LIBIGL

#include <stdexcept>

void IGLSurfaceMesherMC::
mesh(const SignedDistanceRegion<3> &/* sdf */,
     std::vector<MeshIO::IOVertex> &/* vertices */,
     std::vector<MeshIO::IOElement> &/* elements */, const double /* isolevel */) const
{
    throw std::runtime_error("LIBIGL unavailable");
}
#endif // HAS_LIBIGL
