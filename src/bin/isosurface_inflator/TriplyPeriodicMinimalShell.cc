#include <isosurface_inflator/InflatorTypes.hh>
#include <isosurface_inflator/CGALClippedVolumeMesher.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <isosurface_inflator/TriplyPeriodicMinimalShell.hh>
#include <isosurface_inflator/SnapAndReflect.hh>
#include <MeshFEM/filters/remove_dangling_vertices.hh>
#include <MeshFEM/TetMesh.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/OrthotropicHomogenization.hh>

int main(int argc, const char *argv[]) {
    if (argc != 5) {
        std::cerr << "usage:   ./TriplyPeriodicMinimalShell mesher out_path t facet_distance" << std::endl;
        std::cerr << "example: ./TriplyPeriodicMinimalShell cgal out.msh 1.0 1e-3" << std::endl;
        exit(-1);
    }

    std::string mesherName(argv[1]),
                outPath(argv[2]);
    double t = std::stod(argv[3]);
    double facet_distance = std::stod(argv[4]);

    std::vector<Real> A = {        1.0,        1.0,        1.0 },
                 lambda = {    2.0 / t,    2.0 / t,    2.0 / t },
                      P = {          0,          0,          0 };
    std::vector<Vector3d> h = { Vector3d(1, 0, 0), Vector3d(0, 1, 0), Vector3d(0, 0, 1) };
    Real c = 0.25;
    TriplyPeriodicMinimalShell sdfunc(A, h, lambda, P, c);

    std::unique_ptr<MesherBase> mesher;
    if      (mesherName == "cgal") mesher = Future::make_unique<CGALClippedVolumeMesher>();
    else if (mesherName ==  "igl") mesher = Future::make_unique<IGLSurfaceMesherMC>();
    else if (mesherName ==  "midplane") mesher = Future::make_unique<MidplaneMesher>();
    else throw std::runtime_error("Unknown mesher; must be cgal, midplane, or igl");

    mesher->meshingOptions.facetDistance = facet_distance;
    mesher->meshingOptions.cellSize = 1.0;
    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;
    mesher->mesh(sdfunc, vertices, elements);

    remove_dangling_vertices(vertices, elements);

    MeshIO::save("pre_snap.msh", vertices, elements);

    TetMesh<> mesh(elements, vertices.size());
    smartSnap3D(vertices, mesh, sdfunc.boundingBox());

#if 0
    {
        std::vector<MeshIO::IOVertex > reflectedVertices;
        std::vector<MeshIO::IOElement> reflectedElements;
        reflectXYZ(3, vertices, elements, reflectedVertices, reflectedElements);
        std::swap(vertices, reflectedVertices);
        std::swap(elements, reflectedElements);
    }
#endif

    MeshIO::save(outPath, vertices, elements);

    LinearElasticity::HomogenousSimulator<3, 2> sim(elements, vertices);
    std::cout.precision(19);
    std::cout << "Volume:\t" << sim.mesh().volume() << std::endl;

    std::vector<VectorField<Real, 3>> w_ij;
    PeriodicHomogenization::Orthotropic::solveCellProblems(w_ij, sim, 1e-7);
    auto Eh = PeriodicHomogenization::Orthotropic::homogenizedElasticityTensorDisplacementForm(w_ij, sim);
    std::cout << "Homogenized elasticity tensor:" << std::endl;

    auto moduli = Eh.getOrthotropicParameters();

    std::cout << "Moduli:";
    for (Real val : moduli) std::cout << "\t" << val;
    std::cout << std::endl;

    return 0;
}
