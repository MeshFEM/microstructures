#ifndef WIREINFLATOR2D_H
#define WIREINFLATOR2D_H

#include "EdgeMeshPattern.h"
#include "EdgeMeshType.h"
#include "PolyMeshType.h"
#include "PolyMeshUtils.h"
#include "TriMeshType.h"
#include "mshLoader.h"
#include "WireMesh2D.h"
#include "OutMesh.h"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <string>
#include <strings.h>
#include <memory>

template<template<class> class WMesh>
class WireInflator2DImpl;

// This is now an abstract base class so that we can dynamically choose between
// different WMesh template parameters (see WireInflator2DImpl) in the calling
// code.
// 
// This change was made so that Morteza could create an inflator with one WMesh
// type when his deformed cells experienced small shearing and a different WMesh
// type in the large shearing regime.
class WireInflator2D {
public:
    typedef std::shared_ptr<WireInflator2D>		 Ptr;
	typedef OutMesh<2, 3>                        OutMeshType;

    // Constructor passes edgeMeshPath and forwards any additional arguments to
    // WireInflator2DImpl<WMesh>'s constructor.
    template<template<class> class WMesh, typename... Args>
    static Ptr construct(const std::string &edgeMeshPath, Args&&... args) {
        return std::make_shared<WireInflator2DImpl<WMesh>>(edgeMeshPath, std::forward<Args>(args)...);
    }

	virtual void generatePattern(const CellParameters & inP,
	                             const TessellationParameters & inT,
	                             OutMeshType & out,
	                             bool genVelocityField = true) = 0;
	virtual void generateTiledPattern(const Array2D<CellParameters *> & inP,
                                      const TessellationParameters & inT,
                                      OutMeshType & out) = 0;

	/*!
	 * \brief generateQuadsPattern generates a wiremesh pattern over a quad mesh.
	 * \param nodes vector of nodes (coordinates)
	 * \param elements vector of elements (node indices)
	 * \param quadParameters the vector of cell parameters, one per quad.
	 * \param inT the tessellation parameters
	 * \param out the output mesh structure
	 */
	void generateQuadsPattern(const VectorF & nodes,
	                          const VectorI & elements,
	                          const std::vector<CellParameters> & quadParameters,
	                          const TessellationParameters & inT,
	                          OutMeshType & out,
	                          bool averageThicknessOnBoundary = false)
	{
		PolyMesh pmesh;
		bool ok = vectorsToPolyMesh(nodes, elements, pmesh);
		if (!ok)
		{
			std::cout << "Unable to process."<< std::endl
			          << "Input mesh is invalid." << std::endl;
			return;
		}

		generateQuadsPattern(pmesh, quadParameters, inT, out, averageThicknessOnBoundary);
	}

	/*!
	 * \brief generateQuadsPattern generates a wiremesh pattern over a quad mesh.
	 * \param quadMeshPath file path string (can be either an OBJ or MSH file).
	 * \param quadParameters the vector of cell parameters, one per quad.
	 * \param inT the tessellation parameters
	 * \param out the output mesh structure
	 */
	void generateQuadsPattern(const std::string & quadMeshPath,
	                          const std::vector<CellParameters> & quadParameters,
	                          const TessellationParameters & inT,
	                          OutMeshType & out,
	                          bool averageThicknessOnBoundary = false)
	{
		typedef PolyMeshUtils<PolyMesh> PMU;
		PolyMesh pmesh;
		bool ok = false;
		if (checkFileExt(quadMeshPath, ".OBJ"))
		{
			ok = PMU::importFromOBJ(quadMeshPath, pmesh);
			if (ok)
				generateQuadsPattern(pmesh, quadParameters, inT, out, averageThicknessOnBoundary);
		}
		else if (checkFileExt(quadMeshPath, ".MSH"))
		{
			VectorF nodes;
			VectorI elements;

			ok = loadQuadMsh(quadMeshPath, nodes, elements);
			if (ok)
				generateQuadsPattern(nodes, elements, quadParameters, inT, out, averageThicknessOnBoundary);
		}

		if (!ok)
		{
			std::cout << "Unable to process."<< std::endl
			          << "Input mesh is invalid." << std::endl;
		}
	}

	virtual void generateQuadsPattern(PolyMesh & pmesh,
	                                  const std::vector<CellParameters> & quadParameters,
	                                  const TessellationParameters & inT,
	                                  OutMeshType & out,
	                                  bool averageThicknessOnBoundary = false) = 0;

	CellParameters createParameters(void) const { return CellParameters(this->numberOfParameters()); }


	virtual bool parametersValid(const CellParameters & params) const = 0;
    virtual size_t numberOfParameters() const = 0;
	virtual const std::vector<ParameterOperation> &getParameterOperations(void) const = 0;
	virtual std::pair<double,double> getParameterRange(int i) const = 0;
    
    virtual ~WireInflator2D() { }

protected:
	static bool vectorsToPolyMesh(const VectorF & nodes, const VectorI & elements, PolyMesh & pmesh)
	{
		pmesh.Clear();

		static const int NodeDim   = 3;
		static const int FaceArity = 4;

		if ((nodes.size()    % NodeDim   != 0) ||
		    (elements.size() % FaceArity != 0))
			return false;

		// fill nodes
		vcg::tri::Allocator<PolyMesh>::AddVertices(pmesh, nodes.size() / NodeDim);
		for (size_t i=0; i<pmesh.vert.size(); ++i)
		{
			PolyMesh::CoordType & p = pmesh.vert[i].P();
			for (int k=0; k<NodeDim; k++)
			{
				p[k] = nodes[i*NodeDim + k];
			}
		}

		// fill faces
		vcg::tri::Allocator<PolyMesh>::AddFaces(pmesh, elements.size() / FaceArity);
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			PolyMesh::FaceType & f = pmesh.face[i];
			f.Alloc(FaceArity);
			for (int k=0; k<FaceArity; k++)
			{
				int nodeIdx = elements[i*FaceArity + k];
				assert(nodeIdx >= 0 && nodeIdx < pmesh.VN());
				f.V(k) = &pmesh.vert[nodeIdx];
			}
		}

		vcg::tri::UpdateBounding<PolyMesh>::Box(pmesh);

		return true;
	}

	static bool vectorsToEdgeMesh(const VectorF & nodes, const VectorI & elements, EMesh & em)
	{
		em.Clear();

		static const int NodeDim   = 3;
		static const int EdgeArity = 2;

		if ((nodes.size() % NodeDim != 0) ||
		    (elements.size() % EdgeArity != 0))
			return false;

		// fill nodes
		vcg::tri::Allocator<EMesh>::AddVertices(em, nodes.size() / NodeDim);
		for (size_t i=0; i<em.vert.size(); ++i)
		{
			EMesh::CoordType & p = em.vert[i].P();
			for (int k=0; k<NodeDim; k++)
			{
				p[k] = nodes[i*NodeDim + k];
			}
		}

		// fill edges
		vcg::tri::Allocator<EMesh>::AddEdges(em, elements.size() / EdgeArity);
		for (size_t i=0; i<em.edge.size(); ++i)
		{
			EMesh::EdgeType & e = em.edge[i];
			for (int k=0; k<EdgeArity; k++)
			{
				int nodeIdx = elements[i*EdgeArity + k];
				assert(nodeIdx >= 0 && nodeIdx < em.VN());
				e.V(k) = &em.vert[nodeIdx];
			}
		}

		vcg::tri::UpdateBounding<EMesh>::Box(em);

		return true;
	}

	// convert vcg mesh to outmesh structure
	static void vcgMeshToOutMesh(OutMeshType & out, TMesh & m)
	{
		typedef typename OutMeshType::IndexType Index;

		// resize vectors
		out.nodes.resize(m.vert.size());
		out.elements.resize(m.face.size());
		out.edge_fields.clear();

		// fill in nodes
		for (size_t i=0; i<m.vert.size(); ++i)
		{
			// node position
			const TMesh::VertexType & v = m.vert[i];
			out.nodes[i] = { {double(v.cP()[0]), double(v.cP()[1])} };
		}

		// fill in elements
		for (size_t i=0; i<m.face.size(); ++i)
		{
			const TMesh::FaceType & f = m.face[i];
			out.elements[i] = {
			    { Index(vcg::tri::Index(m, f.cV(0))),
			      Index(vcg::tri::Index(m, f.cV(1))),
			      Index(vcg::tri::Index(m, f.cV(2))) }
			};
		}
	}

	static bool checkFileExt(const std::string & filePath, const std::string & ext)
	{
		int diff = (filePath.length() - ext.length());
		if (diff < 0)
			return false;

		const char * cfile = filePath.c_str() + diff;
		const char * cext  = ext.c_str();
		return (strncasecmp(cfile , cext, ext.length()) == 0);
	}

	static bool loadQuadMsh(const std::string & file_name, VectorF & nodes, VectorI & elements)
	{
		try
		{
			MshLoader loader(file_name);
			if (loader.get_nodes_per_element() != 4)
				return false;

			nodes    = loader.get_nodes();
			elements = loader.get_elements();

		} catch (MshLoader::ErrorCode e) { (void)e; return false; }

		return true;
	}

};

template<template<class> class WMesh>
class WireInflator2DImpl : public WireInflator2D {
public:
	typedef EdgeMeshPattern<TMesh, EMesh, WMesh> PatternGen;

    // Constructor passes edgeMeshPath and forwards any additional arguments to
    // PatternGen's constructor.
    template<typename... Args>
    WireInflator2DImpl(const std::string &edgeMeshPath, Args&&... args)
        : m_pattern(edgeMeshPath, std::forward<Args>(args)...) { } 

	WireInflator2DImpl(EMesh & edgeMesh)
	    : m_pattern(edgeMesh) { }

	WireInflator2DImpl( const VectorF & vertices, const VectorI & edges)
	{
		EMesh em;

		bool ok = vectorsToEdgeMesh(vertices, edges, em);
		if (!ok)
		{
			std::cout << "Invalid input wire-mesh."<< std::endl;
			return;
		}

		PatternGen pg(em);
		m_pattern = pg;
	}

	void generatePattern(const CellParameters & inP,
	                     const TessellationParameters & inT,
	                     OutMeshType & out,
	                     bool genVelocityField = true)
	{
		TMesh m;

		m_pattern.params() = { inP, inT };

		m_pattern.generate();
		m_pattern.tessellate(m);

		vcg::tri::Allocator<TMesh>::CompactEveryVector(m);

		vcgMeshToOutMesh(out, m);

		// compute edge velocity field
		// JP: also extract the per-vertex velocity fields
		if (genVelocityField)
			m_pattern.computeVelocityField(m, out.edge_fields, &out.vertex_velocities);
	}

	void generateTiledPattern(const Array2D<CellParameters *> & inP, const TessellationParameters & inT, OutMeshType & out)
	{
		TMesh m;

		m_pattern.params().tessellationParams = inT;

		m_pattern.tile(inP);
		m_pattern.tessellate(m);

		vcg::tri::Allocator<TMesh>::CompactEveryVector(m);

		vcgMeshToOutMesh(out, m);
	}

	void generateQuadsPattern(PolyMesh & pmesh,
	                          const std::vector<CellParameters> & quadParameters,
	                          const TessellationParameters & inT,
	                          OutMeshType & out,
	                          bool averageThicknessOnBoundary = false)
	{
		typedef PolyMeshUtils<PolyMesh> PMU;

		if (!PMU::isQuadMesh(pmesh))
		{
			std::cout << "Unable to process."<< std::endl
			          << "Input mesh is not a quad mesh." << std::endl;
			return;
		}

		if ((pmesh.FN() != int(quadParameters.size())))
		{
			std::cout << "Unable to process."<< std::endl
			          << "Cells and parameters number mismatch." << std::endl;
			return;
		}

		for (size_t i=0; i<quadParameters.size(); ++i)
		{
			if (quadParameters[i].numberOfParameters() != m_pattern.numberOfParameters())
			{
				std::cout << "Unable to process."<< std::endl
				          << "Some quad has an invalid number of parameters." << std::endl;
				return;
			}
		}

		// copy parameters into per face attributes
		auto faceParams = vcg::tri::Allocator<PolyMesh>::GetPerFaceAttribute<CellParameters>(pmesh, m_pattern.CellParametersAttributeName());
		for (size_t i=0; i<pmesh.face.size(); ++i)
		{
			faceParams[i] = quadParameters[i];
		}

		// generate the triangulated mesh
		TMesh m;

		m_pattern.params().tessellationParams = inT;

		m_pattern.generateFromQuads(pmesh, averageThicknessOnBoundary);
		m_pattern.tessellate(m);

		vcg::tri::Allocator<TMesh>::CompactEveryVector(m);

		vcgMeshToOutMesh(out, m);
	}


	const PatternGen & patternGenerator(void)
	{
		return m_pattern;
	}

	bool parametersValid(const CellParameters & params) const { return m_pattern.parametersValid(params); }
    size_t numberOfParameters() const { return m_pattern.numberOfParameters(); }
	const std::vector<ParameterOperation> &getParameterOperations(void) const { return m_pattern.getParameterOperations(); }
	std::pair<double,double> getParameterRange(int i) const { return m_pattern.getParameterRange(i); }

    virtual ~WireInflator2DImpl() { }

private:
	PatternGen m_pattern;
};

#endif // WIREINFLATOR2D_H
