#include "GraphLineMesher.hh"

#include "WireMesh.hh"

using namespace std;

template<class SignedDistanceFunction>
void BoxIntersectionMesher<SignedDistanceFunction>::
mesh(const SignedDistanceFunction &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements)
{
    vector<vector<pair<size_t, size_t>>> polygons;
    vertices.clear();
    boxIntersection1DFeatures(sdf, meshingOptions.marchingSquaresGridSize,
                              vertices, polygons);
    elements.clear(), elements.reserve(polygons.size());
    for (const auto &p : polygons)
        for (const auto &e : p)
            elements.emplace_back(e.first, e.second);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template class GraphLineMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>>>;
// Enable for slower builds...
// template class GraphLineMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>>>;
// template class GraphLineMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>>>;
