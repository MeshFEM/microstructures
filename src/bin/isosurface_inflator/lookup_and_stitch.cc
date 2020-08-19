////////////////////////////////////////////////////////////////////////////////
// lookup_and_stitch.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Using the output of MaterialOptimization_cli, convert each cell's material
//  into a microstructure from an input database and produce a stitched mesh.
//
//  This is currently a quick hack that only works in 3D.
//
//  The input mesh is a tet mesh with cell indices marked (all tets  sharing a
//  cell_index form a hex).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  12/04/2017 16:46:17
////////////////////////////////////////////////////////////////////////////////

#include <isosurface_inflator/StitchedWireMesh.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <MeshFEM/MSHFieldParser.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/TetMesh.hh>
#include <MeshFEM/TriMesh.hh>
#include <MeshFEM/DenseCollisionGrid.hh>
#include <MeshFEM/Utilities/NDArray.hh>
#include <MeshFEM/filters/merge_duplicate_vertices.hh>
#include <MeshFEM/filters/extrude.hh>
#include <MeshFEM/filters/quad_tri_subdiv.hh>
#include <MeshFEM/filters/gen_cursor.hh>
#include <MeshFEM/Triangulate.h>
#include <MeshFEM/util.h>
#include <MeshFEM/Future.hh>
#include <boost/algorithm/string.hpp>
#include <igl/decimate.h>
#include <set>
#include <map>
#include <list>


struct EmbeddedVertex {
    Point3D p;
};

struct DatabaseEntry {
    Real E, nu;
    // ElasticityTensor<Real, 3> complianceTensor;
    Real linfStress;
    size_t topology;
    std::vector<Real> params;
};

static constexpr size_t NONE = std::numeric_limits<size_t>::max();

struct Database {
    // Database lines look like:
    //      file, reduction, init_stress, opt_stress, Y, nu, params
    Database(const std::string &file) {
        std::ifstream dbFile(file);
        if (!dbFile.is_open())
            throw std::runtime_error("Couldn't open microstructure database");

        std::string line;
        std::vector<std::string> lineComponents;
        while (getDataLine(dbFile, line)) {
            boost::trim(line);
            boost::split(lineComponents, line, boost::is_any_of("\t "), boost::token_compress_on);

            DatabaseEntry e;
            e.linfStress = std::stod(lineComponents[3]),
            e.E          = std::stod(lineComponents[4]),
            e.nu         = std::stod(lineComponents[5]);

            std::vector<std::string> resultPathComponents;
            boost::split(resultPathComponents, lineComponents[0], boost::is_any_of("_"));

            e.topology = std::stoi(resultPathComponents[2]);
            for (size_t i = 6; i < lineComponents.size(); ++i)
                e.params.push_back(std::stod(lineComponents[i]));
            entries.emplace_back(std::move(e));
        }
        std::cout << "Using database with " << entries.size() << " entries" << std::endl;
    }

    size_t bestEntry(Real targetE, Real targetNu, Real EThreshold = 0.20, Real nuThreshold = 0.05) {
        size_t bestIdx = NONE;
        Real bestStress = std::numeric_limits<double>::max();
        Real bestDistE = std::numeric_limits<double>::max();
        Real bestDistNu = std::numeric_limits<double>::max();
        for (size_t i = 0; i < entries.size(); ++i) {
            const auto &e = entries[i];
            Real distE  = std::abs((e.E -   targetE) / targetE);
            Real distNu = std::abs((e.nu - targetNu));
            bestDistE = std::min(bestDistE, distE);
            bestDistNu = std::min(bestDistNu, distNu);
            if ((distE  > EThreshold) || (distNu > nuThreshold)) continue;
            if (e.linfStress < bestStress) {
                bestIdx = i;
                bestStress = e.linfStress;
            }
        }
        if (bestIdx == NONE) {
            std::cout << "bestDistE: " << bestDistE << ", " << "bestDistNu: " << bestDistNu << std::endl;
            throw std::runtime_error("Couldn't find material: " + std::to_string(targetE) + ", " + std::to_string(targetNu) + ")");
        }
        return bestIdx;
    }

    std::vector<DatabaseEntry> entries;
};

int main(int argc, const char *argv[]) {
    // TODO: relative moduli error threshold (we'll take the lowest stress design meeting the threshold)
    // TODO: meshing options
    if (argc != 5) {
        std::cerr << "usage: lookup_and_stitch matopt_result.msh database.txt meshingOptions.json out.msh" << std::endl;
        exit(-1);
    }
    const std::string matoptPath(argv[1]),
                      databasePath(argv[2]),
                      meshingOptions(argv[3]),
                      outPath(argv[4]);

    MSHFieldParser<3> fieldParser(matoptPath);

    auto cellIdxSF = fieldParser.scalarField("cell_index", DomainType::PER_ELEMENT);
    auto young     = fieldParser.scalarField("Final E"   , DomainType::PER_ELEMENT);
    auto poisson   = fieldParser.scalarField("Final nu"  , DomainType::PER_ELEMENT);

    using TMesh = TetMesh<EmbeddedVertex>;
    TMesh mesh(fieldParser.elements(), fieldParser.vertices().size());
    for (auto v : mesh.vertices()) v->p = fieldParser.vertices()[v.index()].point;

    size_t numCells = cellIdxSF.max() + 1;
    const size_t nv = mesh.numVertices();
    const size_t numTets = mesh.numTets();
    using BB = BBox<Point3D>;
    std::vector<BB> hexes(numCells);

    using SHandle = TMesh::SHandle<TMesh>;
    auto cellIndex = [&](const SHandle &s) -> size_t { return size_t(cellIdxSF(s.index())); };
    // auto barycenter = [&](const SHandle &s) -> Point3D {
    //     Point3D result(Point3D::Zero());
    //     for (auto v : s.vertices()) result += v->p;
    //     result *= 1.0 / s.numVertices();
    //     return result;
    // };

    for (auto s : mesh.simplices()) hexes.at(cellIndex(s)) = BB(s.vertex(0)->p);
    for (auto s : mesh.simplices()) {
        for (auto v : s.vertices())
            hexes.at(cellIdxSF[s.index()]).unionPoint(v->p);
    }

    using NeigborGrid = NDArray<size_t, 3, 3, 3, 3>;
    std::vector<NeigborGrid> adj(numCells);
    for (auto &a : adj) a.fill(NONE);

    std::vector<std::set<size_t>> cellsAdjVertex(nv);
    std::vector<std::set<size_t>> cellsAdjCell(numCells);

    // Record that v is adjacent to u
    auto addAdjacency = [&](size_t u, size_t v) {
        Point3D cu = hexes[u].center(),
                cv = hexes[v].center();
        size_t i = (cv[0] < cu[0]) ? 0 : ((cv[0] == cu[0]) ? 1 : 2),
               j = (cv[1] < cu[1]) ? 0 : ((cv[1] == cu[1]) ? 1 : 2),
               k = (cv[2] < cu[2]) ? 0 : ((cv[2] == cu[2]) ? 1 : 2);
        size_t old = adj[u](i, j, k);
        assert((old == NONE) || (v == old));
        adj[u](i, j, k) = v;
    };

    // Extract hex grid adjacencies from cell->vertex incidence (since we care
    // about corner neighbors)
    for (auto s : mesh.simplices()) {
        size_t ci = cellIndex(s);
        for (auto v : s.vertices())
            cellsAdjVertex[v.index()].insert(ci);
    }
    for (size_t i = 0; i < nv; ++i) {
        for (size_t u : cellsAdjVertex[i]) {
            for (size_t v : cellsAdjVertex[i]) {
                addAdjacency(u, v);
            }
        }
    }

    using SField = ScalarField<Real>;
    auto tetFieldFromCellField = [&](SField &cellField) {
        SField result(numTets);
        for (auto s : mesh.simplices())
            result(s.index()) = Real(cellField(cellIndex(s)));
        return result;
    };

    auto cellFieldFromTetField = [&](SField &tetField) {
        SField result(numCells);
        for (auto s : mesh.simplices())
            result[cellIndex(s)] = tetField(s.index());
        return result;
    };

    young = cellFieldFromTetField(young);
    poisson = cellFieldFromTetField(poisson);

    // output neighbor counts as scalar field for debugging
    SField valence(numCells);
    for (size_t i = 0; i < numCells; ++i) {
        size_t nn = 0;
        adj[i].visit([&](size_t n) { nn += (n != NONE); });
        valence(i) = nn;
    }

    MSHFieldWriter writer("debug.msh", fieldParser.vertices(), fieldParser.elements());
    writer.addField("valence", tetFieldFromCellField(valence), DomainType::PER_ELEMENT);

    // read pattern database
    Database db(databasePath);

    // lookup cells' E, nu
    SField lutE(numCells), lutNu(numCells), wcs(numCells), topology(numCells);
    std::vector<size_t> dbEntryForCell(numCells);
    for (size_t i = 0; i < numCells; ++i) {
        dbEntryForCell[i] = db.bestEntry(young(i), poisson(i));
        const auto &e = db.entries.at(dbEntryForCell[i]);
        lutE [i] = e.E;
        lutNu[i] = e.nu;
        wcs  [i] = e.linfStress;
        topology[i] = e.topology;
    }

    // output looked up pattern's topology, E, nu / error as scalar field
    writer.addField("E", tetFieldFromCellField(young), DomainType::PER_ELEMENT);
    writer.addField("nu", tetFieldFromCellField(poisson), DomainType::PER_ELEMENT);
    writer.addField("lut E", tetFieldFromCellField(lutE), DomainType::PER_ELEMENT);
    writer.addField("lut nu", tetFieldFromCellField(lutNu), DomainType::PER_ELEMENT);
    writer.addField("wcs", tetFieldFromCellField(wcs), DomainType::PER_ELEMENT);
    writer.addField("topology", tetFieldFromCellField(topology), DomainType::PER_ELEMENT);

    using WMP = std::shared_ptr<WireMeshBase>;
    std::map<size_t, WMP> topologies;
    topologies.emplace(  24, std::make_shared<WireMesh<Symmetry::Orthotropic<>>>("/home/jpanetta/Research/microstructures/patterns/3D/reference_wires/pattern0024.wire"));
    topologies.emplace(  77, std::make_shared<WireMesh<Symmetry::Orthotropic<>>>("/home/jpanetta/Research/microstructures/patterns/3D/reference_wires/pattern0077.wire"));
    topologies.emplace( 646, std::make_shared<WireMesh<Symmetry::Orthotropic<>>>("/home/jpanetta/Research/microstructures/patterns/3D/reference_wires/pattern0646.wire"));
    topologies.emplace( 746, std::make_shared<WireMesh<Symmetry::Orthotropic<>>>("/home/jpanetta/Research/microstructures/patterns/3D/reference_wires/pattern0746.wire"));
    topologies.emplace(1053, std::make_shared<WireMesh<Symmetry::Orthotropic<>>>("/home/jpanetta/Research/microstructures/patterns/3D/reference_wires/pattern1053.wire"));
    topologies.emplace(1065, std::make_shared<WireMesh<Symmetry::Orthotropic<>>>("/home/jpanetta/Research/microstructures/patterns/3D/reference_wires/pattern1065.wire"));

    // Create geometry for each cell
    auto mesher = Future::make_unique<IGLSurfaceMesherMC>();
    mesher->meshingOptions.load(meshingOptions);
    std::vector<MeshIO::IOVertex > globalVertices;
    std::vector<MeshIO::IOElement> globalElements;
    for (size_t ci = 0; ci < numCells; ++ci) {
        std::cout << "processing cell " << ci << " (topology " << topology[ci] << ")" <<  std::endl;
        NDCubeArray<WMP, 3, 3> topologyGrid;
        NDCubeArray<std::vector<Real>, 3, 3> parameterGrid;
        // Fill in the topology and parameter grids for the cell + any present neighbors
        adj[ci].visit([&](size_t cell, size_t i, size_t j, size_t k) {
            if (cell == NONE) return;
            const auto &e = db.entries.at(dbEntryForCell[cell]);
            topologyGrid(i, j, k) = topologies.at(e.topology);
            parameterGrid(i, j, k) = e.params;
        });
        std::cout << "building swm" << std::endl;
        auto swm = make_stitched_wire_mesh<3, false>(topologyGrid);
        std::cout << "params from param grid" << std::endl;
        auto params = swm.paramsFromParamGrid(parameterGrid);
        std::cout << "build sdf" << std::endl;
        PatternSignedDistance<double, StitchedWireMesh<3, false>> sdf(swm);
        std::cout << "set parameters" << std::endl;
        sdf.setParameters(params, Eigen::Matrix3d::Identity(), mesher->meshingOptions.jointBlendingMode);
#if 1
        // Mesh each cell geometry inside the period cell
        std::vector<MeshIO::IOVertex > cellVertices;
        std::vector<MeshIO::IOElement> cellElements;
        mesher->mesh(sdf, cellVertices, cellElements);
        MeshIO::save("cell_" + std::to_string(ci) + ".msh", cellVertices, cellElements);
#else
        std::vector<MeshIO::IOVertex > cellVertices;
        std::vector<MeshIO::IOElement> cellElements;
        MeshIO::load("cell_" + std::to_string(ci) + ".msh", cellVertices, cellElements);
#endif

        // Bilinarly map each meshed geometry to fill the hex
        const size_t vtxOffset = globalVertices.size();
        for (auto &v : cellVertices) globalVertices.push_back(hexes[ci].interpolatePoint(0.5 * (v.point + Point3D(1, 1, 1))));
        for (auto e : cellElements) {
            for (size_t &vi : e)  vi += vtxOffset;
            globalElements.emplace_back(std::move(e));
        }
    }

    // Stitch the mesh together
    std::cout << "merging" << std::endl;
    std::vector<MeshIO::IOVertex > stitchedVertices;
    std::vector<MeshIO::IOElement> stitchedElements;
    merge_duplicate_vertices(globalVertices, globalElements, stitchedVertices, stitchedElements, 1e-8);
    // stitchedVertices = globalVertices;
    // stitchedElements = globalElements;

    TriMesh<EmbeddedVertex> surface(stitchedElements, stitchedVertices.size());
    for (auto v : surface.vertices()) v->p = stitchedVertices[v.index()].point;

    std::cout << "filling holes" << std::endl;
    // Fill holes/create the support structure.
    // Extract the boundary loops
    std::vector<bool> visited(surface.numBoundaryVertices(), false);
    std::vector<MeshIO::IOVertex > supportColVertices;
    std::vector<MeshIO::IOElement> supportColElements;
    size_t numLoops = 0;
    for (auto bv : surface.boundaryVertices()) {
        if (visited.at(bv.index())) continue;
        ++numLoops;
        std::vector<Point3D> bdryPts;
        do {
            visited[bv.index()] = true;
            bdryPts.push_back(bv.volumeVertex()->p);
            bv = bv.outEdge().tip();
        } while (!visited[bv.index()]);
        BB bb(bdryPts);
        auto dim = bb.dimensions();
        assert(std::abs(dim[0] * dim[1] * dim[2]) == 0);

        // Map to plane without caring about orientation for now
        size_t zeroComp;
        for (zeroComp = 0; zeroComp < 3; ++ zeroComp)
            if (dim[zeroComp] == 0.0) break;

        size_t xComp = (zeroComp + 1) % 3;
        size_t yComp = (zeroComp + 2) % 3;

        std::list<std::list<Point2D>> polygons(1);
        auto &poly = polygons.back();
        for (const auto &bp : bdryPts)
            poly.emplace_back(bp[xComp], bp[yComp]);

        // Note: boundary is traversed ccw, so signed area will
        // be positive if the mapping to the plane preserved orientation,
        // negative if it flipped orientation.
        Real signedArea = area(poly);
        if (std::abs(signedArea) < 1e-6) {
            std::cout << "signedArea: " << signedArea << std::endl;
            std::cout << "poly.size(): " << poly.size() << std::endl;
            for (auto &p : poly) {
                std::cout << "\t\t" << p.transpose() << std::endl;
            }

        }
        // assert(std::abs(signedArea) > 1e-6);
        if (signedArea < 0) {
            std::swap(xComp, yComp);
            for (auto &p : poly)
                std::swap(p[0], p[1]);
        }

        std::vector<MeshIO::IOVertex > fillVertices;
        std::vector<MeshIO::IOElement> fillElements;
        // Fill the hole with Triangle
        triangulatePSLG(polygons, std::vector<Point2D>(), fillVertices, fillElements, 100, "QY");

        for (size_t i = 0; i < bdryPts.size(); ++i) {
            assert(bdryPts[i][xComp] == fillVertices[i].point[0]);
            assert(bdryPts[i][yComp] == fillVertices[i].point[1]);
            assert(bdryPts[i][zeroComp] == bb.minCorner[zeroComp]);
        }

        size_t offset = stitchedVertices.size();
        for (const auto &v : fillVertices) {
            Point3D pt(Point3D::Zero());
            pt[xComp] = v[0];
            pt[yComp] = v[1];
            pt[zeroComp] = bb.minCorner[zeroComp];
            stitchedVertices.emplace_back(pt);
        }
        for (auto e : fillElements) {
            for (size_t &vi : e)  vi += offset;
            stitchedElements.emplace_back(std::move(e));
        }

        // Detect the parts needing support columns.
        // (support geometry could be pre-computed and then transformed to the desired loc)
        if ((zeroComp == 2) && std::abs(bb.minCorner[2] - BB(fieldParser.vertices()).minCorner[2]) < 1e-14) {
            Point3D center(Point3D::Zero());
            for (const auto &v : fillVertices) {
                Point3D pt;
                pt[xComp] = v[0];
                pt[yComp] = v[1];
                pt[zeroComp] = bb.minCorner[zeroComp];
                center += pt;
            }
            center *= 1.0 / fillVertices.size();

            const Real diameter = 0.4;
            const size_t nSubdiv = 15;
            const Real height = 5.0;

            // Generate a circle
            std::vector<MeshIO::IOVertex > colVertices;
            std::vector<MeshIO::IOElement> colElements;
            colVertices.push_back(center);
            for (size_t i = 0; i < nSubdiv; ++i) {
                Real theta = 2 * M_PI * i / (nSubdiv + 1);
                colVertices.push_back((center + Point3D(diameter/2 * cos(theta), diameter/2 * sin(theta), 0.0)).eval());
                colElements.emplace_back(0, 1 + i, 1 + (i + 1) % nSubdiv);
            }
            TriMesh<EmbeddedVertex> circMesh(colElements, colVertices.size());
            for (auto v : circMesh.vertices()) v->p = colVertices[v.index()].point;
            extrude(circMesh, height, colVertices, colElements);
            std::vector<MeshIO::IOVertex > colVerticesTri;
            std::vector<MeshIO::IOElement> colElementsTri;
            std::vector<size_t> quadIdx;
            quad_tri_subdiv(colVertices, colElements, colVerticesTri, colElementsTri, quadIdx);
            const size_t offset = supportColVertices.size();
            for (auto &v : colVerticesTri) supportColVertices.emplace_back(v);
            for (auto e : colElementsTri) {
                for (size_t &vi : e)  vi += offset;
                supportColElements.emplace_back(std::move(e));
            }
        }

    }

    std::cout << "filled " << numLoops << " holes" << std::endl;

    // Stitch the hole-closing patches on.
    std::cout << "merging hole infills" << std::endl;
    merge_duplicate_vertices(stitchedVertices, stitchedElements, stitchedVertices, stitchedElements, 0);

    TriMesh<EmbeddedVertex> closedMesh(stitchedElements, stitchedVertices.size());
    for (auto v : closedMesh.vertices()) v->p = stitchedVertices[v.index()].point;

    if (closedMesh.numBoundaryVertices() > 0) {
        std::cout << "ERROR: mesh not closed." << std::endl;
        std::vector<MeshIO::IOVertex > cursorVertices;
        std::vector<MeshIO::IOElement> cursorElements;
        for (auto bv: closedMesh.boundaryVertices()) {
            gen_cursor(1.0, bv.volumeVertex()->p, cursorVertices, cursorElements);
        }
        MeshIO::save("cursors.msh", cursorVertices, cursorElements);
    }

    MeshIO::save(outPath, stitchedVertices, stitchedElements);
    MeshIO::save("support.off", supportColVertices, supportColElements);

    // Run mesh decimation
    Eigen::MatrixXd V(stitchedVertices.size(), 3), U;
    Eigen::MatrixXi F(stitchedElements.size(), 3), G;
    Eigen::VectorXi J, I;
    for (size_t i = 0; i < stitchedVertices.size(); ++i)
        V.row(i) = stitchedVertices[i].point;
    for (size_t i = 0; i < stitchedElements.size(); ++i) {
        F(i, 0) = stitchedElements[i][0];
        F(i, 1) = stitchedElements[i][1];
        F(i, 2) = stitchedElements[i][2];
    }

    const size_t maxFaces = 4000000;;
    std::cout << "Decimating..." << std::endl;
    igl::decimate(V, F, maxFaces, U, G, J, I);

    std::cout << "Outputting decimated mesh" << std::endl;
    std::vector<MeshIO::IOVertex > decimatedVertices(U.rows());
    std::vector<MeshIO::IOElement> decimatedElements(G.rows());
    for (size_t i = 0; i < decimatedVertices.size(); ++i)
        decimatedVertices[i].point = U.row(i);
    for (size_t i = 0; i < decimatedElements.size(); ++i)
        decimatedElements[i] = MeshIO::IOElement(G(i, 0), G(i, 1), G(i, 2));
    MeshIO::save("decimated4M.msh", decimatedVertices, decimatedElements);

    V = U;
    F = G;
    igl::decimate(V, F, 2000000, U, G, J, I);

    decimatedVertices.resize(U.rows());
    decimatedElements.resize(G.rows());
    for (size_t i = 0; i < decimatedVertices.size(); ++i)
        decimatedVertices[i].point = U.row(i);
    for (size_t i = 0; i < decimatedElements.size(); ++i)
        decimatedElements[i] = MeshIO::IOElement(G(i, 0), G(i, 1), G(i, 2));
    MeshIO::save("decimated2M.msh", decimatedVertices, decimatedElements);


    return 0;
}
