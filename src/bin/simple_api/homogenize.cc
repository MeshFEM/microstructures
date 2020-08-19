////////////////////////////////////////////////////////////////////////////////
// homogenize.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Generates and homogenizes cells that have been deformed linearly.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  05/12/2015 15:22:14
////////////////////////////////////////////////////////////////////////////////

#include <MeshFEM/SymmetricMatrix.hh>
#include <MeshFEM/util.h>
#include <MeshFEM/Types.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/FEMMesh.hh>
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>
#include <MeshFEM/Materials.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/filters/remove_dangling_vertices.hh>
#include <CLI/CLI.hpp>
#include <json.hpp>
#include <string>
#include <vector>

using json = nlohmann::json;
using namespace PeriodicHomogenization;

////////////////////////////////////////////////////////////////////////////////

// Default arguments
struct Args {
	std::string mesh;
	bool homogenize = false;
	bool transformVersion = false;
	std::string material;
	std::string jacobian;
	std::string displacedMesh;
	double displacementScale = 0.1;
	bool parametrizedTransform = false;
	int degree = 2;
	std::string tile;
	std::string out;
	std::string dumpJson;
};

////////////////////////////////////////////////////////////////////////////////

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<typename T>
void dumpJson(const T &EhDefo, const std::string &filename) {
	json data;
	data["elasticity_tensor"] = json::array();
	auto E = EhDefo.getCoefficients();
	for (auto x : E) {
		data["elasticity_tensor"].push_back(x);
	}
	data["homogenized_moduli"] = json::array();
	std::vector<double> P;
	EhDefo.getOrthotropicParameters(P);
	for (auto x : P) {
		data["homogenized_moduli"].push_back(x);
	}

	// Write json
	std::ofstream out(filename);
	out << data;
}

template<size_t _N, size_t _FEMDegree>
void execute(const Args &args,
			 const std::vector<MeshIO::IOVertex> &inVertices,
			 const std::vector<MeshIO::IOElement> &inElements)
{
	auto &mat = HMG<_N>::material;
	if (!args.material.empty()) {
		mat.setFromFile(args.material);
	}
	typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
	typedef LinearElasticity::Simulator<Mesh> Simulator;
	typedef typename Simulator::VField VField;
	Simulator sim(inElements, inVertices);
	auto &mesh = sim.mesh();

	std::cout << std::setprecision(16);

	if (args.parametrizedTransform) {
		if (!args.homogenize || !args.tile.empty()) {
			throw std::runtime_error("parametrizedTransform only supports homogenization");
		}
		if (!args.transformVersion) {
			std::cerr << "WARNING: running transformVersion" << std::endl;
		}
		if (_N != 2) {
			throw std::runtime_error("parametrizedTransform only supports 2D");
		}
		std::string line;
		Real theta, lambda;

		auto EBase = mat.getTensor();

		while (getDataLine(std::cin, line)) {
			std::istringstream ss(line);
			if (!(ss >> theta >> lambda)) {
				throw std::runtime_error("invalid input transformation: " + line);
			}

			Eigen::Matrix<Real, _N, _N> rot, stretch, jacobian;
			rot << cos(theta), -sin(theta),
				   sin(theta),  cos(theta);
			stretch << lambda, 0,
					   0,      1;
			jacobian = rot * stretch * rot.transpose();
			std::cout << jacobian << std::endl;
			mat.setTensor(EBase.transform(jacobian.inverse()));

			std::vector<VField> w_ij;
			solveCellProblems(w_ij, sim);
			auto EhDefo = homogenizedElasticityTensorDisplacementForm(w_ij, sim).transform(jacobian);
			auto ShDefo = EhDefo.inverse();

			std::cout << theta << '\t' << lambda << '\t'
				 << EhDefo.D(0, 0) << '\t' << EhDefo.D(0, 1) << '\t' << EhDefo.D(0, 2) << '\t'
										   << EhDefo.D(1, 1) << '\t' << EhDefo.D(1, 2) << '\t'
																	 << EhDefo.D(2, 2) << '\t'
				 << ShDefo.D(0, 0) << '\t' << ShDefo.D(0, 1) << '\t' << ShDefo.D(0, 2) << '\t'
										   << ShDefo.D(1, 1) << '\t' << ShDefo.D(1, 2) << '\t'
																	 << ShDefo.D(2, 2) << std::endl;
			if (!args.dumpJson.empty()) {
				dumpJson(EhDefo, args.dumpJson);
			}
		}

		// BENCHMARK_REPORT();
		return;
	}

	// Parse jacobian.
	Eigen::Matrix<Real, _N, _N> jacobian;

	std::vector<double> jacobianComponents;
	{
		std::istringstream ss(args.jacobian);
		double x;
		while (ss >> x) {
			jacobianComponents.push_back(x);
		}
	}
	if (jacobianComponents.size() != _N * _N) {
		throw std::runtime_error("Invalid deformation jacobian");
	}
	for (size_t i = 0; i < _N; ++i) {
		for (size_t j = 0; j < _N; ++j) {
			jacobian(i, j) = jacobianComponents[_N * i + j];
		}
	}

	auto bbox = mesh.boundingBox();
	VectorND<_N> center = 0.5 * (bbox.minCorner + bbox.maxCorner);
	std::vector<MeshIO::IOVertex> deformedVertices;

	for (size_t vi = 0; vi < mesh.numVertices(); ++vi) {
		VectorND<_N> p = mesh.vertex(vi).node()->p;
		deformedVertices.emplace_back((jacobian * (p - center)).eval());
	}

	Real deformedCellVolume = bbox.volume() * jacobian.determinant();

	if (args.homogenize && args.transformVersion) {
		std::vector<VField> w_ij;
		// Morteza's transformation formulas
		mat.setTensor(mat.getTensor().transform(jacobian.inverse()));
		solveCellProblems(w_ij, sim);
		auto EhDefo = homogenizedElasticityTensorDisplacementForm(w_ij, sim).transform(jacobian);
		std::cout << "Elasticity tensor:" << std::endl;
		std::cout << EhDefo << std::endl << std::endl;
		std::cout << "Homogenized Moduli: ";
		EhDefo.printOrthotropic(std::cout);
		if (!args.dumpJson.empty()) {
			dumpJson(EhDefo, args.dumpJson);
		}
	} else if (args.homogenize) {
		sim.applyPeriodicConditions();
		sim.applyNoRigidMotionConstraint();
		sim.setUsePinNoRigidTranslationConstraint(true);
		sim.updateMeshNodePositions(deformedVertices);
		std::shared_ptr<MSHFieldWriter> writer;
		if (!args.out.empty()) {
			writer = std::make_shared<MSHFieldWriter>(args.out, mesh);
		}
		std::vector<VField> w_ij;
		typedef typename Simulator::SMatrix SMatrix;
		constexpr size_t numStrains = SMatrix::flatSize();
		for (size_t i = 0; i < numStrains; ++i) {
			VField rhs(sim.constantStrainLoad(-SMatrix::CanonicalBasis(i)));
			w_ij.push_back(sim.solve(rhs));
			if (writer) {
				writer->addField("load_ij " + std::to_string(i), sim.dofToNodeField(rhs), DomainType::PER_NODE);
				writer->addField("w_ij" + std::to_string(i), w_ij.back(), DomainType::PER_NODE);
				writer->addField("strain w_ij " + std::to_string(i), sim.averageStrainField(w_ij[i]),
								 DomainType::PER_ELEMENT);
			}
		}
		auto EhDefo = homogenizedElasticityTensorDisplacementForm(w_ij, sim, deformedCellVolume);
		std::cout << "Elasticity tensor:" << std::endl;
		std::cout << EhDefo << std::endl << std::endl;
		std::cout << "Homogenized Moduli: ";
		EhDefo.printOrthotropic(std::cout);

		if (!args.displacedMesh.empty()) {
			std::string out_mesh = args.displacedMesh;
			SymmetricMatrixValue<Real, _N> strain;
			auto bbox = mesh.boundingBox();
			VectorND<_N> center = bbox.center();

			std::vector<VField> cstrainDisp_ij;

			for (unsigned index = 0; index < 3; index++) {
				VField cstrainDisp(mesh.numNodes());
				strain = SMatrix::CanonicalBasis(index);
				for (auto n : mesh.nodes()) {
					cstrainDisp(n.index()) = strain.contract(n->p - center);
				}

				cstrainDisp_ij.push_back(cstrainDisp);
			}

			for (unsigned index = 0; index < 3; index++) {
				for (auto n : mesh.nodes()) {
					cstrainDisp_ij[index](n.index()) += w_ij[index](n.index());
				}

				writer->addField("u_cstrain_ij" + std::to_string(index), cstrainDisp_ij[index], DomainType::PER_NODE);
			}

			// Adding displacement to each node:
			for (unsigned index = 0; index < 3; index++) {
				std::vector<MeshIO::IOVertex> displacedVertices;
				VField cstrainDisp = cstrainDisp_ij[index];
				Eigen::Matrix<Real, _N, Eigen::Dynamic> cstrainDisp_v = sim.template nodeToVertexField<VField>(
						cstrainDisp).data();
				assert((long int) mesh.numVertices() == cstrainDisp_v.cols());
				for (size_t vi = 0; vi < mesh.numVertices(); ++vi) {
					VectorND<_N> p = mesh.vertex(vi).node()->p;
					VectorND<_N> cstrainDisp_vi = cstrainDisp_v.col(vi);
					displacedVertices.emplace_back((p + args.displacementScale * cstrainDisp_vi).eval());
				}

				sim.updateMeshNodePositions(displacedVertices);
				std::shared_ptr<MSHFieldWriter> writerDisplacedMesh;

				writerDisplacedMesh = std::make_shared<MSHFieldWriter>(out_mesh + std::to_string(index) + ".msh", mesh);
			}
		}

		if (!args.dumpJson.empty()) {
			dumpJson(EhDefo, args.dumpJson);
		}
	} else if (!args.tile.empty()) {
		mesh.setNodePositions(deformedVertices);
		std::vector<int> tileComponents;
		{
			std::istringstream ss(args.tile);
			int x;
			while (ss >> x) {
				tileComponents.push_back(x);
			}
		}
		if (tileComponents.size() != _N) {
			throw std::runtime_error("Invalid number of tiling dimensions");
		}
		std::vector<size_t> tilings;
		int numCells = 1;
		for (const auto &ci : tileComponents) {
			if (ci <= 0) { throw std::runtime_error("Invalid number of tilings"); }
			tilings.push_back((size_t) ci);
			numCells *= ci;
		}

		std::vector<MeshIO::IOVertex>  tiledVertices;
		std::vector<MeshIO::IOElement> tiledElements;
		tiledVertices.reserve(numCells * inVertices.size());
		tiledElements.reserve(numCells * inElements.size());
		size_t numCellVertices = inVertices.size();

		if (_N == 2) tilings.push_back(1);

		// Glue in (i, j, k) scanline order, attaching a copy of the base cell to the
		// existing structure by merging duplicated vertices. These duplicated
		// vertices for neighbors along dimension d are specified by the entries in
		// identifiedFacePairs[d]. In this scanline order, only neighbors
		// (i - 1, j, k), (i, j - 1, k), (i, j, k - 1) >= (0, 0, 0) exist.
		//      *-------*
		//     /   /   /|
		//    /---+---/ |
		//   /   /   /|/|     ^ j (y)
		//  *---+---* + *     |
		//  |   |   |/|/      |     i (x)
		//  |---+---| /       *----->
		//  |   |   |/       /
		//  *-------*       v k (z)
		// Gluing is done by overwriting tets' "tiled vertex indices" (i.e. index
		// before merging duplicates) with the corresponding "glued vertex indices"
		// (index after merging duplicates). For vertices that end up merging, the
		// final glued index is stored in a dictionary mapping the vertices'
		// original tiled vertex indices to the index they are merged to.
		//
		// Note: This scanline approach only works on geometry that is actually
		// triply periodic. For instance, it failed on James' zigzag shape that is
		// missing geometry in one of its corners.

		// Update a global tiled vertex index -> glued vertex index map.
		// Tiled vertex indieces are of the form:
		//   base_vertex_index + vtxOffset(i, j, k)
		//
		auto vtxIdxOffset = [tilings,numCellVertices](int i, int j, int k) -> size_t {
			return numCellVertices * (i * tilings[2] * tilings[1] + j * tilings[2] + k);
		};

		VectorND<_N> delta;
		for (size_t i = 0; i < tilings[0]; ++i) {
			delta[0] = i * bbox.dimensions()[0];
			for (size_t j = 0; j < tilings[1]; ++j) {
				delta[1] = j * bbox.dimensions()[1];
				for (size_t k = 0; k < tilings[2]; ++k) {
					if (_N > 2) delta[2] = k * bbox.dimensions()[2];
					auto vtxPosOffset = (jacobian * delta).eval();
					for (const auto &v : deformedVertices)
						tiledVertices.emplace_back((truncateFrom3D<VectorND<_N>>(v.point) + vtxPosOffset).eval());
					for (auto e : inElements) {
						for (size_t ei = 0; ei < e.size(); ++ei)
							e[ei] = e[ei] + vtxIdxOffset(i, j, k);
						tiledElements.push_back(e);
					}
				}
			}
		}

		// TODO: merge duplicated vertices.
		remove_dangling_vertices(tiledVertices, tiledElements);
		save(args.out, tiledVertices, tiledElements);
	}
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
	Args args;

	// Parse arguments
	CLI::App app{"homogenize"};

	app.add_option("mesh", args.mesh, "input mesh")->required()->check(CLI::ExistingFile);
	auto homogenize_opt = app.add_flag("--homogenize", args.homogenize, "run homogenization");
	app.add_flag("--transformVersion", args.homogenize, "use transform version of homogenization");
	app.add_option("-m,--material", args.material, "base material");
	app.add_option("-j,--jacobian", args.jacobian, "linear deformation jacobian");
	app.add_option("--displacedMesh", args.displacedMesh, "file prefix containing displaced mesh");
	app.add_option("--displacementScale", args.displacementScale, "used to scale displacements obtained from constant plus periodic strain");
	app.add_flag("-p,--parametrizedTransform", "read a list of parameterized deformations from stdin");
	app.add_option("-d,--degree", args.degree, "degree of finite elements");
	app.add_option("-t,--tile", args.tile, "tilings 'nx ny nz' (default: 1)")->excludes(homogenize_opt);
	app.add_option("-o,--out", args.out, "output file of deformed geometry (and w_ij fields if homogenization is run)");
	app.add_option("--dumpJson", args.dumpJson, "dump info into a json file)");

	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		return app.exit(e);
	}

	if (app.count("--out") + app.count("--homogenize") == 0) {
		return app.exit(CLI::ValidationError("No operation requested"));
	}

	if (app.count("--jacobian") + app.count("--parametrizedTransform") != 1) {
		return app.exit(CLI::ValidationError("Must specify either deformation jacobian or parametrizedTransform"));
	}

	std::vector<MeshIO::IOVertex>  inVertices;
	std::vector<MeshIO::IOElement> inElements;
	auto type = load(args.mesh, inVertices, inElements, MeshIO::FMT_GUESS, MeshIO::MESH_GUESS);

	// Infer dimension from mesh type.
	size_t dim;
	if (type == MeshIO::MESH_TET) dim = 3;
	else if (type == MeshIO::MESH_TRI) dim = 2;
	else throw std::runtime_error("Mesh must be triangle or tet.");

	// Look up and run appropriate instantiation.
	auto exec = (dim == 3) ? ((args.degree == 2) ? execute<3, 2> : execute<3, 1>)
						   : ((args.degree == 2) ? execute<2, 2> : execute<2, 1>);
	exec(args, inVertices, inElements);

	return 0;
}
