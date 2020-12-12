#include "BoxIntersectionMesher.hh"

#include <vector>
#include <utility>

#include "BoxIntersection1DFeatures.hh"

using namespace std;

void BoxIntersectionMesher::
mesh(const SignedDistanceRegion<3> &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements) const
{
    vector<vector<pair<size_t, size_t>>> polygons;
    vertices.clear();
    boxIntersection1DFeatures(sdf, meshingOptions.marchingSquaresGridSize,
                              meshingOptions.marchingSquaresCoarsening,
                              vertices, polygons);
    elements.clear(), elements.reserve(polygons.size());
    for (const auto &p : polygons)
        for (const auto &e : p)
            elements.emplace_back(e.first, e.second);
}
