#ifndef WIREMESHEMBEDDING_H
#define WIREMESHEMBEDDING_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cstdlib>
#include <math.h>

template <class EMesh, class PolyMesh>
class WireMeshEmbedding
{
public:
	typedef WireMeshEmbedding<EMesh, PolyMesh> ThisType;

	typedef typename EMesh::CoordType          ECoordType;

	typedef typename PolyMesh::CoordType       PCoordType;
	typedef typename PolyMesh::FaceType        PFaceType;
	typedef typename PolyMesh::FacePointer     PFacePointer;\


	struct QuadParametrization
	{
		// these represents the indices of the vertices
		// defining the positive X-axis (V(index0) ->V(index1))
		char index0;

		PCoordType interpolate(const vcg::Point2d & t, const PFaceType & face) const
		{
			char index1 = (index0+1)%4;
			char index2 = (index0+2)%4;
			char index3 = (index0+3)%4;
			      PCoordType   p0 = face.cP(index0);
			      PCoordType   p1 = face.cP(index1);
			const PCoordType & p2 = face.cP(index2);
			const PCoordType & p3 = face.cP(index3);

			// bilinear interpolation
			p0 = p0 * (1 - t.X()) + p1 * (t.X());
			p1 = p3 * (1 - t.X()) + p2 * (t.X());
			return p0 * (1 - t.Y()) + p1 * (t.Y());
		}
	};

	typedef typename PolyMesh::template PerFaceAttributeHandle<QuadParametrization> QuadParamHandle;

	static QuadParamHandle getQuadParametrizationHandle(PolyMesh & pmesh)
	{
		QuadParamHandle h = vcg::tri::Allocator<PolyMesh>::template FindPerFaceAttribute<QuadParametrization>(pmesh, ThisType::ParametrizationAttributeName());
		assert(vcg::tri::Allocator<PolyMesh>::IsValidHandle(pmesh, h));
		return h;
	}

	// used to planarize and obtain a coherent mapping of edge meshes among quads
	static void preprocessQuadMesh(PolyMesh & pmesh)
	{
		if (pmesh.VN() == 0 || pmesh.FN() == 0)
			return;

		// planarize
		planarizeQuadMesh(pmesh);

		// compute adjacencies
		vcg::tri::UpdateTopology<PolyMesh>::FaceFace(pmesh);
		// vf adj does not exist for polygonal meshes
		auto vf = vcg::tri::Allocator<PolyMesh>::template GetPerVertexAttribute<std::vector<PFacePointer> >(pmesh, VFAttributeName());
		for (auto fi=pmesh.face.begin(); fi!=pmesh.face.end(); ++fi)
		{
			for (int i=0; i<fi->VN(); ++i)
				vf[fi->V(i)].push_back(&(*fi));
		}

		// compute border
		vcg::tri::UpdateFlags<PolyMesh>::FaceBorderFromFF(pmesh);

		// create a coherent (if possibile) parametrization for the quad mesh
		createParametrization(pmesh);

	}

	// MHS on Feb 8 2016
	// for viewing the parametrization sequence
	static void dumpParametrizationSequence(PolyMesh & pmesh)
	{
		cout << "the parametrization sequence is: ";
		for (auto f = pmesh.face.begin(); f != pmesh.face.end(); ++f){
			QuadParametrization qpar = getQuadParametrizationHandle(pmesh)[f];
			cout << (int) qpar.index0;
		}
		cout << endl;
	}

	// assume em is a normalized edge mesh
	static void embedWireMesh(EMesh & em, const PFaceType & f, PolyMesh & pmesh)
	{
		QuadParametrization qpar = getQuadParametrizationHandle(pmesh)[&f];
		for (size_t i=0; i<em.vert.size(); ++i)
		{
			ECoordType & p = em.vert[i].P();
			vcg::Point2d t(p[0], p[1]);
			p = ECoordType::Construct(qpar.interpolate(t, f));
		}
	}

	static std::vector<PFacePointer> adjacentFaceEdge(PolyMesh & pmesh, const PFaceType & f, unsigned char edgeIdx)
	{
		std::vector<PFacePointer> ret;
		if (edgeIdx >= 0 && edgeIdx < f.VN())
		{
			QuadParamHandle qpar = getQuadParametrizationHandle(pmesh);
			unsigned char actualIndex = (qpar[&f].index0 + edgeIdx) % 4;
			PFacePointer af = f.FFp(actualIndex);
			if (af != &f)
				ret.push_back(af);
		}
		return ret;
	}

	static std::vector<char> idxToElementInAdjacentFaceEdge(PolyMesh & pmesh, const PFaceType & f, unsigned char edgeIdx)
	{
		std::vector<char> ret;
		if (edgeIdx >= 0 && edgeIdx < f.VN())
		{
			QuadParamHandle qpar = getQuadParametrizationHandle(pmesh);
			unsigned char actualIndex = (qpar[&f].index0 + edgeIdx) % 4;
			PFacePointer af = f.FFp(actualIndex);
			char		 adjFaceEdgeIdx = f.cFFi(actualIndex);

			if (af != &f)
			{
				char actualAdjFaceEdgeIdx = (adjFaceEdgeIdx + qpar[af].index0) % 4;
				ret.push_back(actualAdjFaceEdgeIdx);
			}
		}
		return ret;
	}
	
	static std::vector<PFacePointer> adjacentFaceVertex(PolyMesh & pmesh, const PFaceType & f, unsigned char vertIdx)
	{
		std::vector<PFacePointer> ret;
		if (vertIdx >= 0 && vertIdx < f.VN())
		{
			QuadParamHandle qpar = getQuadParametrizationHandle(pmesh);
			unsigned char actualIndex = (qpar[&f].index0 + vertIdx) % 4;

			auto vf = vcg::tri::Allocator<PolyMesh>::template GetPerVertexAttribute<std::vector<PFacePointer> >(pmesh, VFAttributeName());
			std::vector<PFacePointer> & adjFaces = vf[f.V(actualIndex)];
			for (PFacePointer af : adjFaces)
			{
				if (af != &f)
					ret.push_back(af);
			}
		}
		return ret;
	}

	// MHS on AUG14, 2015:
	// this is not being used, right now. It must be updated in case it is needed in the future!
	static std::vector<char> idxToElementInAdjacentFaceVertex(PolyMesh & /* pmesh */, const PFaceType & /* f */, unsigned char /* vertIdx */)
	{
		std::vector<char> ret;
		return ret;
	}

	// MHS on FEB4, 2016:
	// this function returns the idxs vertex coordinate in the equivalent parallelogram obtained form the
	// quad [p1 p2 p3 p4]
	static PCoordType getEquivalentCoords(const PCoordType p1, // coordinate of the first  vertex in the quad
										  const PCoordType p2, // coordinate of the second vertex in the quad
										  const PCoordType p3, // coordinate of the third  vertex in the quad
										  const PCoordType p4, // coordinate of the fourth vertex in the quad
										  char idx) // index 0..3
	{
		if (idx < 0 || idx > 3)
			throw("index out of range");

		switch (idx){
		case 0:
			return (p1 + p1 + p1) / 4.0 + (p2 - p3 + p4) / 4.0;
			break;
		case 1:
			return (p2 + p2 + p2) / 4.0 + (p3 - p4 + p1) / 4.0;
			break;
		case 2:
			return (p3 + p3 + p3) / 4.0 + (p4 - p1 + p2) / 4.0;
			break;
		case 3:
			return (p4 + p4 + p4) / 4.0 + (p1 - p2 + p3) / 4.0;
			break;
		}	
	}

	// MHS on JUL14, 2015:
	// This method returns the deformation (F) corresponding to the equivalent parallelogram for each
	// bilinear quad in the quad mesh.
	static void  getJacobians(PolyMesh 						& pmesh, 
			                  std::vector<Eigen::Matrix2d>	& defs,
			                  std::vector<double> 			& angles)
	{
		defs.clear();
		angles.clear();


		for (auto fc = pmesh.face.begin(); fc != pmesh.face.end(); fc++)
		{
			QuadParametrization qpar = getQuadParametrizationHandle(pmesh)[fc];

			char index0 = qpar.index0;
			char index1 = (index0 + 1) % 4;
			char index2 = (index0 + 2) % 4;
			char index3 = (index0 + 3) % 4;

			PCoordType p1 = fc->cP(index0);
			PCoordType p2 = fc->cP(index1);
			PCoordType p3 = fc->cP(index2);
			PCoordType p4 = fc->cP(index3);

			PCoordType direction = p2 - p1;

			double cosAngle = direction[0] / sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
			double angle    = acos(cosAngle);

			PCoordType p1_equivalent = getEquivalentCoords(p1, p2, p3, p4, 0); 
			PCoordType p2_equivalent = getEquivalentCoords(p1, p2, p3, p4, 1);
			PCoordType p3_equivalent = getEquivalentCoords(p1, p2, p3, p4, 2);
			PCoordType p4_equivalent = getEquivalentCoords(p1, p2, p3, p4, 3);

			PCoordType alpha = (p2_equivalent - p1_equivalent); 
			PCoordType beta  = (p4_equivalent - p1_equivalent); 

			Eigen::Matrix2d jacobian;
			jacobian(0, 0) = alpha[0];
			jacobian(0, 1) = beta[0];
			jacobian(1, 0) = alpha[1];
			jacobian(1, 1) = beta[1];
			
			jacobian = jacobian / std::sqrt(std::abs(jacobian.determinant()));

			defs.push_back(jacobian); 
			angles.push_back(angle);
		}
	}

	// MHS on AUG10 2015:
	// This method returns the symmetric right stretc (U)
	static void getStretches(PolyMesh & pmesh,
			                 std::vector<Eigen::Matrix2d>	& stretches,
			                 std::vector<double> 			& angles)
	{
		std::vector<Eigen::Matrix2d> defs;
		getJacobians(pmesh, defs, angles);
		
		
		stretches.clear();
		for (size_t i = 0; i < defs.size(); ++i)
		{
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(defs[i] * defs[i].transpose());
			stretches.push_back(es.operatorSqrt());
		}
	}

	// MHS on Feb 16, 2016:
	// This is to overwrite Luigi's coherent parametrization
	static void createLocalParametrization(PolyMesh & pmesh)
	{
		vcg::tri::RequireFFAdjacency(pmesh);
		vcg::tri::RequirePerFaceFlags(pmesh);

		// clear visited flag
		vcg::tri::UpdateFlags<PolyMesh>::FaceClearV(pmesh);

		QuadParamHandle handle =
		        vcg::tri::Allocator<PolyMesh>::template GetPerFaceAttribute<QuadParametrization>(pmesh, ThisType::ParametrizationAttributeName());

		// each face is parameterized independently, i.e., index0 is set such that
		// ||P_eq(index1) - P_eq(index0)|| > ||P_eq(index2) - P_eq(index1)
		// where P_eq[i] (i=0..3) are the coordinates of the vertexes of the equivalent parallelogram
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			PFacePointer fp = &pmesh.face[i];
			PCoordType p1, p2, p3, p4,
					   p1_eq, p2_eq, p3_eq;

			p1 = fp->cP(0);
			p2 = fp->cP(1);
			p3 = fp->cP(2);
			p4 = fp->cP(3);

			p1_eq = getEquivalentCoords(p1, p2, p3, p4, 0);
			p2_eq = getEquivalentCoords(p1, p2, p3, p4, 1);
			p3_eq = getEquivalentCoords(p1, p2, p3, p4, 2);

			PCoordType dir_1 = p2_eq - p1_eq;
			PCoordType dir_2 = p3_eq - p2_eq;
		
			QuadParametrization & qp = handle[fp];


			if (dir_1.Norm() < dir_2.Norm())
				qp.index0 = 1;
			else
				qp.index0 = 0;

			fp->SetV();
		}
	}

private:
	static std::string ParametrizationAttributeName(void)
	{
		return "parametrization";
	}

	static std::string VFAttributeName(void)
	{
		return "vf_adj";
	}

	static void planarizeQuadMesh(PolyMesh & pmesh)
	{
		for (auto vi=pmesh.vert.begin(); vi!=pmesh.vert.end(); vi++)
		{
			vi->P()[2] = 0;
		}
	}




	static void createParametrization(PolyMesh & pmesh)
	{
		vcg::tri::RequireFFAdjacency(pmesh);
		vcg::tri::RequirePerFaceFlags(pmesh);

		// clear visited flag
		vcg::tri::UpdateFlags<PolyMesh>::FaceClearV(pmesh);

		QuadParamHandle handle =
		        vcg::tri::Allocator<PolyMesh>::template GetPerFaceAttribute<QuadParametrization>(pmesh, ThisType::ParametrizationAttributeName());

		// propagate if we encounter faces not visited
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			PFacePointer fp = &pmesh.face[i];
			if (!fp->IsV())
			{
				// set the parametrization for the first face
				setOptimalParametrization(handle, fp);

				// propagate to all other quads in the same connected component
				std::vector<PFacePointer> toPropagate;
				toPropagate.push_back(fp);

				for (size_t k=0; k<toPropagate.size(); ++k)
				{
					PFacePointer f = toPropagate[k];
					for (char e=0; e<4; ++e)
					{
						PFacePointer fAdj = setAdjacentParametrization(handle, f, e);
						if (fAdj)
							toPropagate.push_back(fAdj);
					}
				}
			}
		}
	}

	static void setOptimalParametrization(QuadParamHandle handle, PFacePointer f)
	{
		// TODO CHANGE TO OPTIMAL PARAMETRIZATION
		// right now, it gets the first quad and assumes a configuration like this:
		//
		//        Y
		//        ^
		//        |
		//  3-----------2
		//  |     |     |
		//  |     |     |
		//  |     +---- | --> X
		//  |           |
		//  |           |
		//  0-----------1

		assert(f);

		QuadParametrization & qp = handle[f];
		qp.index0 = 0;

		f->SetV();
	}

	static PFacePointer setAdjacentParametrization(QuadParamHandle & handle, const PFacePointer face, char adjIdx)
	{
		assert(adjIdx>=0 && adjIdx<4);

		PFaceType * adjFace = face->FFp(adjIdx);
		if (adjFace != NULL && !adjFace->IsV())
		{
			char edgeAdj = face->cFFi(adjIdx);
			handle[adjFace].index0 = (handle[face].index0 - adjIdx + (edgeAdj+2) + 4)%4;

			adjFace->SetV();
			return adjFace;
		}
		return NULL;
	}
};

#endif // WIREMESHEMBEDDING_H
