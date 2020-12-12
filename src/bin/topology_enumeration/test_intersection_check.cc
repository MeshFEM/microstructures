#include "intersection_check.hh"

using namespace std;

using Edge = std::pair<size_t, size_t>;

int main(int argc, char *argv[])
{
    string path = argv[1];

    vector<MeshIO::IOVertex> inVertices;
    vector<MeshIO::IOElement> elements;
    MeshIO::load(path, inVertices, elements);

    vector<Point3d> vertices;
    vector<Edge> edges;
    vertices.reserve(inVertices.size());
    edges.reserve(elements.size());
    for (auto &v : inVertices) vertices.emplace_back(v.point);
    for (auto &e : elements) {
        if (e.size() != 2) throw runtime_error("Non-edge element encountered.");
        edges.push_back({e[0], e[1]});
    }

    cout << "Has self intersection: " << hasSelfIntersection(vertices, edges) << endl;

    return 0;
}
