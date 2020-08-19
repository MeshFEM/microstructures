#ifndef WIREMESH2D_H
#define WIREMESH2D_H

#include "EdgeMeshUtils.h"
#include "InflatorParameters.h"

#include <string>
#include <utility>
#include <cmath>
#include <map>
#include <unordered_map>
#include <cassert>

#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>

template <class EMesh>
class WireMesh2D
{
public:
	typedef typename EMesh::VertexType      VertexType;
	typedef typename EMesh::EdgeType        EdgeType;
	typedef typename VertexType::ScalarType ScalarType;
	typedef typename VertexType::CoordType  CoordType;

	WireMesh2D(void) {;}

	WireMesh2D(EMesh & em)
	{
		this->setMesh(em);
	}

	WireMesh2D(const std::string & edgeMeshPath)
	{
		this->setMesh(edgeMeshPath);
	}

	WireMesh2D & operator = (WireMesh2D & wm)
	{
		m_operations      = wm.m_operations;
		m_symmetry_orbits = wm.m_symmetry_orbits;
		m_orbits_idx      = wm.m_orbits_idx;
		m_radius_range    = wm.m_radius_range;
		m_transl_range    = wm.m_transl_range;
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;

		return *this;
	}

	WireMesh2D(const WireMesh2D & wm)
	    : m_operations(wm.m_operations)
	    , m_symmetry_orbits(wm.m_symmetry_orbits)
	    , m_orbits_idx(wm.m_orbits_idx)
	    , m_radius_range(wm.m_radius_range)
	    , m_transl_range(wm.m_transl_range)
	{
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;
	}

	void setMesh(const std::string & edgeMeshPath)
	{
		m_operations.clear();
		bool done = EdgeMeshUtils<EMesh>::importObj(m_em, edgeMeshPath);

		if (!done)
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
			return;
		}

		this->setup();
	}

	void setMesh(EMesh & em)
	{
		m_operations.clear();
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em, em, false, true);
		this->setup();
	}

	bool isValid(void) const
	{
		return (m_em.VN() > 0) && (m_em.EN() > 0);
	}

	size_t numberOfParameters(void) const
	{
		return m_operations.size();
	}

	CellParameters createCellParameters(void) const
	{
		CellParameters params(this->numberOfParameters());
		for (size_t i=0; i<params.numberOfParameters(); ++i)
		{
			std::pair<double,double> range = this->getParameterRange(i);
			params.parameter(i) = (range.first + range.second) / 2;
		}

		return params;
	}

	/*!
	 * \brief parametersValid checks wheter the current set of parameters is valid for the current wire mesh.
	 * \return true if the parameters are valid, false otherwise.
	 */
	bool parametersValid(const CellParameters & params) const
	{
		if (params.numberOfParameters() != m_operations.size())
			return false;

		std::pair<double,double> range_radius = this->getRadiusParameterRange();
		std::pair<double,double> range_transl = this->getTranslationParameterRange();
		for (size_t i=0; i<m_operations.size(); ++i)
		{
			switch (m_operations.at(i).type)
			{
			case ParameterOperation::Radius:
				if (params.cParameter(i) < range_radius.first ||
				    params.cParameter(i) > range_radius.second)
					return false;
				break;
			case ParameterOperation::Translation:
				if (params.cParameter(i) < range_transl.first ||
				    params.cParameter(i) > range_transl.second)
					return false;
				break;
			default:
				assert(0);
			}
		}
		return true;
	}

	std::pair<double,double> getParameterRange(int i) const
	{
		assert(i>=0 && i<int(m_operations.size()));

		switch (m_operations.at(i).type)
		{
		case ParameterOperation::Radius:
			return this->getRadiusParameterRange();
		case ParameterOperation::Translation:
			return this->getTranslationParameterRange();
		default:
			assert(0);
		}
		return std::make_pair(0,0);
	}

	void getUnmodifiedNormalizedEdgeMesh(EMesh & em)
	{
		if (!this->isValid())
			return;

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(em, m_normalized_em, false, true);
	}

	void getUnmodifiedEdgeMesh(EMesh & em)
	{
		if (!this->isValid())
			return;

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(em, m_em, false, true);
	}

	void getEdgeMesh(EMesh & em, const CellParameters & params)
	{
		if (!this->isValid())
			return;

		assert(this->parametersValid(params));

		this->getUnmodifiedEdgeMesh(em);
		this->applyParameterOperations(em, params);
	}

	void getNormalizedMesh(EMesh & em, const CellParameters & params)
	{
		if (!this->isValid())
			return;

		assert(this->parametersValid(params));

		this->getUnmodifiedNormalizedEdgeMesh(em);
		this->applyParameterOperations(em, params);
	}

	vcg::Point2<ScalarType> getOriginalScale(void) const
	{
		return vcg::Point2<ScalarType>(m_bbox.DimX(), m_bbox.DimY());
	}

	const std::vector<ParameterOperation> & getParameterOperations(void) const
	{
		// all symmetry orbit radius parameters come first
		return m_operations;
	}

	int getOrbitIndexForNode(int node)
	{
		if (!this->isValid() || node < 0 || node >= m_em.VN())
			return -1;

		auto tag = vcg::tri::Allocator<EMesh>::template FindPerVertexAttribute<int>(m_em, SymmetryOrbitAttributeName());
		if (!vcg::tri::Allocator<EMesh>::IsValidHandle(m_em, tag))
			return -1;

		return m_orbits_idx[tag[node]];
	}

    virtual ~WireMesh2D() { }

protected:
	EMesh                           m_em;
	EMesh                           m_normalized_em;
	vcg::Box3<ScalarType>           m_bbox; // the original bounding box
	std::vector<ParameterOperation> m_operations;
	std::vector<std::vector<int> >  m_symmetry_orbits;
	std::unordered_map<int,int>     m_orbits_idx;
	std::pair<double,double>        m_radius_range;
	std::pair<double,double>        m_transl_range;

	static inline ScalarType Epsilon(void) { return 0.00001; }

	void setup(void)
	{
		this->m_radius_range = std::make_pair( 0.01, 5.0 );
		this->m_transl_range = std::make_pair(-0.18, 0.18);

		if (!this->isValid())
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
			return;
		}

		vcg::tri::Allocator<EMesh>::CompactEveryVector(m_em);
		vcg::tri::UpdateTopology<EMesh>::VertexEdge(m_em);
		vcg::tri::UpdateTopology<EMesh>::EdgeEdge(m_em);

		this->normalizeMesh();

		this->computeSymmetryOrbits();
		this->generateParameterOperations(m_symmetry_orbits);
	}

	virtual std::pair<double,double> getRadiusParameterRange(void) const
	{
		return m_radius_range;
	}

	virtual std::pair<double,double> getTranslationParameterRange(void) const
	{
		return m_transl_range;
	}

	static std::string SymmetryOrbitAttributeName(void)
	{
		return "symmetry_orbit";
	}

	void normalizeMesh(void)
	{
		if (!this->isValid())
			return;

		// store bounding box for later use
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);
		m_bbox = m_em.bbox;

		vcg::tri::UpdatePosition<EMesh>::Translate(m_em, -m_em.bbox.min);
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);

		assert(m_em.bbox.Dim()[0] != 0 && m_em.bbox.Dim()[1] != 0);

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, m_em, false, true);

		// real normalization
		for (size_t i=0; i<m_normalized_em.vert.size(); ++i)
		{
			CoordType & p = m_normalized_em.vert[i].P();
			p[0] /= m_normalized_em.bbox.max[0];
			p[1] /= m_normalized_em.bbox.max[1];
			p[2] = 0;
		}
		vcg::tri::UpdateBounding<EMesh>::Box(m_normalized_em);
	}

	// compute symmetry orbits clusters
	void computeSymmetryOrbits()
	{
		if (!isValid())
			return;

		// tag each vertex
		auto tag = vcg::tri::Allocator<EMesh>::template AddPerVertexAttribute<int>(m_em, SymmetryOrbitAttributeName());
		for (int i=0; i<int(m_em.vert.size()); ++i)
		{
			tag[i] = i;
		}

		// symmetrize wrt to x and y and find matching vertices
		std::map<int, std::vector<int> > orbits;
		for (int ax=0; ax<2; ++ax) // do the matching for x and then y axis
		{
			for (int i=0; i<int(m_em.vert.size()); ++i)
			{
				CoordType flipped = m_em.vert[i].cP();
				flipped[ax] = (m_em.bbox.Dim()[ax] - flipped[ax]);

				for (int j=0; j<int(m_em.vert.size()); ++j)
				{
					if ( (i == j) || (tag[i] == tag[j]) )
						continue;

					CoordType delta = vcg::Abs(flipped - m_em.vert[j].cP());

					if (delta[0] <= Epsilon() &&
					    delta[1] <= Epsilon() &&
					    delta[2] <= Epsilon())
					{
						tag[j] = tag[i] = std::min(tag[i], tag[j]);
					}
				}
			}
		}
		for (int i=0; i<int(m_em.vert.size()); ++i)
		{
			orbits[tag[i]].push_back(i);
		}

		m_symmetry_orbits.clear();
		m_orbits_idx.clear();
		int index = 0;
		for (auto & o : orbits)
		{
			m_symmetry_orbits.push_back(o.second);
			m_orbits_idx[o.first] = index++;
		}
	}

	void generateParameterOperations(const std::vector<std::vector<int> > & symmetryOrbits)
	{
		m_operations.clear();

		// generate radius change for each symmetry orbit
		{
			ParameterOperation op;
			op.type = ParameterOperation::Radius;
			for (auto & s : symmetryOrbits)
			{
				op.nodes = s;
				m_operations.push_back(op);
			}
		}

		// generate displacement change
		{
			ParameterOperation op;
			op.type = ParameterOperation::Translation;
			const auto & meshBbox = m_em.bbox;
			for (auto & s : symmetryOrbits)
			{
				vcg::Box3<ScalarType> bbox;
				for (int i : s)
				{
					bbox.Add(m_em.vert[i].P());
				}

				// add x and y displacement operations
				for (int ax=0; ax<2; ++ax)
				{
					if (bbox.Dim()[ax] < meshBbox.Dim()[ax] &&
					    bbox.Dim()[ax] > Epsilon())
					{
						vcg::Point2d displ(0,0);
						for (int i : s)
						{
							displ[ax] = double(this->sign(m_em.vert[i].P()[ax] - meshBbox.Center()[ax]));
							op.nodes_displ[i] = displ;
						}
						m_operations.push_back(op);
						op.nodes_displ.clear();
					}
				}
			}
		}
	}

	void applyParameterOperations(EMesh & em, const CellParameters & params)
	{
		assert(m_operations.size() == params.numberOfParameters());

		for (size_t i=0; i<m_operations.size(); ++i)
		{
			const ParameterOperation & p = m_operations[i];
			switch (p.type)
			{
			// store vertex radius into quality
			case ParameterOperation::Radius:
				for (int node : p.nodes)
				{
					em.vert[node].Q() = params.cParameter(i);
				}
				break;
			// displace the vertex
			case ParameterOperation::Translation:
				for (auto & displ : p.nodes_displ)
				{
					em.vert[displ.first].P() +=
					        CoordType(displ.second[0] * em.bbox.Dim()[0], displ.second[1] * em.bbox.Dim()[1], 0) * params.cParameter(i);
				}
				break;
			default:
				assert(0);
			}
		}
	}

	template <typename T>
	static int sign(T x)
	{
		return (T(0) < x) - (x < T(0));
	}
};


// MHS on JUL29, 2015:
// this is to streamline the symmetry operation calculations
template <class EMesh>
class metaEMesh
{
public:
	typedef typename EMesh::VertexType::CoordType CoordType;
	typedef typename EMesh::VertexType::ScalarType ScalarType;

	metaEMesh(EMesh & em)
	{
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em, em, false, true);
		// initialize other memebers here

		m_radius_ops.clear();
		m_transl_ops.clear();
		ParameterOperation op;
		vcg::Point2d dx(1.0, 0.0);
		vcg::Point2d dy(0.0, 1.0);

		for (size_t i = 0; i < em.vert.size(); ++i)
		{
			op.type = ParameterOperation::Radius;
			op.nodes.clear();
			op.nodes_displ.clear();
			op.nodes.push_back(i);
			m_radius_ops.push_back(op);

			op.type = ParameterOperation::Translation;
			op.nodes.clear();
			op.nodes_displ.clear();
			op.nodes_displ[i] = dx;
			m_transl_ops.push_back(op);

			op.type = ParameterOperation::Translation;
			op.nodes.clear();
			op.nodes_displ.clear();
			op.nodes_displ[i] = dy;
			m_transl_ops.push_back(op);
		}

		//addUpDownSymmetry();
		//addLeftRightSymmetry();
		//addPosDiagSymmetry();
		//addNegDiagSymmetry();
		//addBaseSymmetry();

		//displayRadiusOps();
		//displayTranslOps();
	}

	// symmetry about y = 0
	void addUpDownSymmetry()
	{
		addAboutAxisSymmetry(0);
	}

	// symmetry about x = 0
	void addLeftRightSymmetry()
	{
		addAboutAxisSymmetry(90);
	}

	// symmetry about y = x
	void addPosDiagSymmetry()
	{
		addAboutAxisSymmetry(45);
	}

	// symmetry about y = -x
	void addNegDiagSymmetry()
	{
		addAboutAxisSymmetry(135);
	}

	std::vector<ParameterOperation> getOperations() // radius ops first and then translation ones
	{
		std::vector<ParameterOperation> ops;
		ops.clear();

		for (auto op : m_radius_ops)
			if (op.nodes.size() > 0)
				ops.push_back(op);

		for (auto op : m_transl_ops)
			if (op.nodes_displ.size() > 0)
				ops.push_back(op);

		return ops;
	}

	void displayRadiusOps()
	{
		std::cout << std::endl << "---------- displaying radius ops ----------" << std::endl;
		int counter = 1;
		for (auto op : m_radius_ops)
		{
			if (op.nodes.size() > 0)
			{
				std::cout << "there are " << op.nodes.size() << " nodes in operation " << counter << ", and they are:" << std::endl;
				for (int i = 0; i < op.nodes.size(); ++i)
					std::cout << op.nodes[i] << "\t";
				std::cout << std::endl << endl;
			}
			++counter;
		}

	}

	void displayTranslOps()
	{
		std::cout << std::endl << "---------- displaying translation ops ----------" << std::endl;
		int counter = 1;
		for (auto op : m_transl_ops)
		{
			if (op.nodes_displ.size() > 0)
			{
				std::cout << "there are " << op.nodes_displ.size() << " nodes_displs in operation " << counter << ", and they are:" << std::endl;
				for (auto nd = op.nodes_displ.begin(); nd != op.nodes_displ.end(); ++nd)
					std::cout << nd->first << ":(" << nd->second[0] << "," << nd->second[1] << ")\t";
				std::cout << std::endl << endl;
			}
			++counter;
		}

	}


	void addBaseSymmetry()
	{
		for (size_t i = 0; i < m_em.vert.size(); ++i)
		{
			for (int flipAxis = 0; flipAxis < 2; ++flipAxis) // 0 corresponds to the x-axis and 1 corresponds to the y-axis
			{
			    vcg::Point2d  noDispDirection(0.0, 0.0);
				vcg::Point2d simDispDirection(0.0, 0.0);
				simDispDirection[(flipAxis+1)%2] = 1.0;
				 noDispDirection[(flipAxis+2)%2] = 1.0;

				int  firstAxisDegree = (flipAxis == 0) ? 0 : 90;
				int secondAxisDegree = (flipAxis == 0) ? 90 : 0;

				size_t j = findNodeIdx(flipPoint(m_em.vert[i].cP(), firstAxisDegree));
				if (j == i)
				{
					int k = findNodeIdx(flipPoint(m_em.vert[j].cP(), secondAxisDegree));

					// correct the radius symmetries
					int currentOp = findRadiusOperationIdx(j);
					int flippedOp = findRadiusOperationIdx(k);
					if (currentOp != flippedOp)
					while(!m_radius_ops[flippedOp].nodes.empty())
					{
						m_radius_ops[currentOp].nodes.push_back(m_radius_ops[flippedOp].nodes.back());
						m_radius_ops[flippedOp].nodes.pop_back();
					}

					vcg::Box3<ScalarType> bBox = m_em.bbox;
					if (vcg::math::Abs(m_em.vert[i].cP()[flipAxis]) <= Epsilon() ||
						vcg::math::Abs(m_em.vert[i].cP()[flipAxis] - bBox.Dim()[flipAxis]) <= Epsilon())
					{
						//std::cout<< std::endl << "F I AM HERE  " << i << "   REALLY HERE" << std::endl << std::endl;

						// correct the translation symmetries
						std::pair<int, ScalarType> currentOpSignPair = findTranslationOperationIdx(j, simDispDirection);
						int        currentOpIdx = std::get<0>(currentOpSignPair);
						//ScalarType currentSign  = std::get<1>(currentOpSignPair);
						if (currentOpIdx >= 0)
							m_transl_ops[currentOpIdx].nodes_displ.clear();
						// correct the translation symmetries
						currentOpSignPair = findTranslationOperationIdx(j, noDispDirection);
						currentOpIdx = std::get<0>(currentOpSignPair);
						//currentSign  = std::get<1>(currentOpSignPair);
						if (currentOpIdx >= 0)
							m_transl_ops[currentOpIdx].nodes_displ.clear();
					}

				}
			}
			if (i == 8 || i == 10)
			{
				std::pair<int, ScalarType> currentOpSignPair = findTranslationOperationIdx(i, vcg::Point2d(0.0, 1.0));
				int currentOpIdx = std::get<0>(currentOpSignPair);
				if (currentOpIdx > -1)
					m_transl_ops[currentOpIdx].nodes_displ.clear();
				continue;
			}
			if (i == 9 || i == 11)
			{
				std::pair<int, ScalarType> currentOpSignPair = findTranslationOperationIdx(i, vcg::Point2d(1.0, 0.0));
				int currentOpIdx = std::get<0>(currentOpSignPair);
				if (currentOpIdx > -1)
					m_transl_ops[currentOpIdx].nodes_displ.clear();
				continue;
			}
		}
	}

	static size_t findOperationIdx(const int nodeID, const std::vector<ParameterOperation> & allOps)
	{
		// this method assumes that the raduis operations are first in the list

		// extract raduis ops
		std::vector<ParameterOperation> radiusOps;
		radiusOps.clear();
		for (auto op : allOps)
		{
			if (op.type == ParameterOperation::Radius)
				radiusOps.push_back(op);
		}

		int idx = -1;
		int counter = 0;
		for (auto op : radiusOps)
		{
			std::vector<int>::iterator it;
			it = find(op.nodes.begin(), op.nodes.end(), nodeID);
			if (it != op.nodes.end())
			{
				idx = counter;
				break;
			}
			++counter;
		}
		return idx;
	}

protected:
	EMesh 								m_em;
	std::vector<ParameterOperation>		m_radius_ops;
	std::vector<ParameterOperation> 	m_transl_ops;

	static inline ScalarType Epsilon(void) { return 0.00001; }

	void addAboutAxisSymmetry(const int axisDegree)
	{
		// add the radius operations
		addRadiusSymmetries(axisDegree);

		// add the translation operationsn
		// swip to add "dx" symmetries
		addAboutCoordAxisSymmetries(axisDegree, "dx");
		// swip to add "dy" symmetries
		addAboutCoordAxisSymmetries(axisDegree, "dy");
	}




	void addRadiusSymmetries(const int axisDegree)
	{
		// add the radius symmetries
		for (size_t i = 0; i < m_em.vert.size(); ++i)
		{
			size_t currentNode = i;
			size_t flippedNode = findNodeIdx(flipPoint(m_em.vert[i].cP(), axisDegree));

			if (currentNode != flippedNode)
			{
				int currentOp = findRadiusOperationIdx(currentNode);
				int flippedOp = findRadiusOperationIdx(flippedNode);

				if (currentOp != flippedOp)
					while(!m_radius_ops[flippedOp].nodes.empty())
					{
						m_radius_ops[currentOp].nodes.push_back(m_radius_ops[flippedOp].nodes.back());
						m_radius_ops[flippedOp].nodes.pop_back();
					}
			}
		}
	}

	void addAboutCoordAxisSymmetries(const int axisDegree, const std::string axis)
	{
		for (size_t i = 0; i < m_em.vert.size(); ++i)
		{
			std::pair<int, ScalarType>	 currentOpSignPair, flippedOpSignPair;
			int							 currentOpIdx,      flippedOpIdx;
			ScalarType					 currentSign,       flippedSign;

			vcg::Point2d currentDisp(0.0, 0.0);
			if (axis == "dx")
				currentDisp[0] = 1.0;
			else
				currentDisp[1] = 1.0;


			// ------------------------
			// for the current node and disp find:
			// a) the "index" in the operation list
			// b) the "sign"
			// such that:
			// i:sign*disp is in m_transl_ops[index]
			// vcg::Point2d currentDisp(1.0, 0.0);
			currentOpSignPair = findTranslationOperationIdx(i, currentDisp);
			currentOpIdx = std::get<0>(currentOpSignPair);
			currentSign  = std::get<1>(currentOpSignPair);

			// if didn't find move to next node
			if (currentOpIdx < 0)
				continue;

			// find the flipped node and flipped disp
			int j = findNodeIdx(flipPoint(m_em.vert[i].cP(), axisDegree));
			vcg::Point2d flippedDisp = flipDisp(currentSign * currentDisp, axisDegree);

			// for node j and flipped disp find:
			// a) the "index" in the operation list
			// b) the "sign"
			// such that:
			// j:sign*flippedDisp is in m_transl_ops[index]
			flippedOpSignPair = findTranslationOperationIdx(j, flippedDisp);
			flippedOpIdx = std::get<0>(flippedOpSignPair);
			flippedSign  = std::get<1>(flippedOpSignPair);

			std::map<int, vcg::Point2d>::iterator it;
			if (flippedOpIdx < 0)
			{
				// clear all the map items located in currentOpIdx
				m_transl_ops[currentOpIdx].nodes_displ.clear();
				continue;
			}

			else if (flippedOpIdx == currentOpIdx && flippedSign > 0)
			{
				if (flippedSign <= Epsilon())
					m_transl_ops[currentOpIdx].nodes_displ[i] += flippedSign * flippedDisp;
				else
					continue;
			}
			else if (flippedOpIdx == currentOpIdx && flippedSign < 0)
			{
				// clear the map located in flippedIpIdx
				m_transl_ops[flippedOpIdx].nodes_displ.clear();
				continue;
			}
			else
			{
				// append the map items located in flippedOpIdx to the map items located in currentOpIcs
				for (auto nd = m_transl_ops[flippedOpIdx].nodes_displ.begin(); nd != m_transl_ops[flippedOpIdx].nodes_displ.end(); ++nd)
				{
					std::map<int, vcg::Point2d>::iterator it;
					it = m_transl_ops[currentOpIdx].nodes_displ.find(nd->first);
					if (it != m_transl_ops[currentOpIdx].nodes_displ.end())
						m_transl_ops[currentOpIdx].nodes_displ[nd->first] += flippedSign * nd->second;
					else
						m_transl_ops[currentOpIdx].nodes_displ[nd->first]  = flippedSign * nd->second;
				}
				// clear the maps located in flippedOpIdx
				m_transl_ops[flippedOpIdx].nodes_displ.clear();
				continue;
			}
		}
	}

	// Helper Functions
	CoordType flipPoint(CoordType inPoint, const int axisDegree)
	{
		vcg::Box3<ScalarType> bBox = m_em.bbox;
		CoordType flipped;
		switch(axisDegree)
		{
			case 0:
				flipped[0] = inPoint[0];
				flipped[1] = bBox.Dim()[1] - inPoint[1];
				break;
			case 90:
				flipped[0] = bBox.Dim()[0] - inPoint[0];
				flipped[1] = inPoint[1];
				break;
			case 45:
				flipped[0] = inPoint[1];
				flipped[1] = inPoint[0];
				break;
			case 135:
				flipped[0] = bBox.Dim()[0] - inPoint[1];
				flipped[1] = bBox.Dim()[1] - inPoint[0];
				break;
		} // end of switch on axisDegree
		return flipped;
	}

	vcg::Point2d flipDisp(vcg::Point2d inDisp, const int axisDegree)
	{
		vcg::Point2d flipped;
		switch(axisDegree)
		{
			case 0:
				flipped[0] =   inDisp[0];
				flipped[1] = - inDisp[1];
				break;
			case 90:
				flipped[0] = - inDisp[0];
				flipped[1] =   inDisp[1];
				break;
			case 45:
				flipped[0] =   inDisp[1];
				flipped[1] =   inDisp[0];
				break;
			case 135:
				flipped[0] = - inDisp[1];
				flipped[1] = - inDisp[0];
				break;
		} // end of switch on axisDegree
		return flipped;
	}


	int findNodeIdx(CoordType inPoint)
	{
		int idx, i;
		bool flag = false;

		for (i = 0; i < int(m_em.vert.size()); ++i)
		{
			CoordType delta = vcg::Abs(m_em.vert[i].cP() - inPoint);
			if (delta[0] <= Epsilon() &&
				delta[1] <= Epsilon() &&
				delta[2] <= Epsilon())
			{
				flag = true;
				break;
			}
		}

		if (flag)
			idx = i;
		else
			return -1;

		return idx;
	}




	int findRadiusOperationIdx(const int nodeID)
	{
		int idx = -1;
		int counter = 0;
		for (auto op : m_radius_ops)
		{
			std::vector<int>::iterator it;
			it = find(op.nodes.begin(), op.nodes.end(), nodeID);
			if (it != op.nodes.end())
			{
				idx = counter;
				break;
			}
			++counter;
		}
		return idx;
	}

	std::pair<int, ScalarType> findTranslationOperationIdx(const int nodeID, vcg::Point2d disp)
	{
		int idx = -1;
		ScalarType sign = 0;
		int counter = 0;
		for (auto op : m_transl_ops)
		{
			int flag = 0;
			for (auto nd = op.nodes_displ.begin(); nd != op.nodes_displ.end(); ++nd)
			{
				if (nd->first == nodeID &&
					vcg::math::Abs(disp.dot(nd->second)) > Epsilon())
				{
					flag = 1;
					idx = counter;
					sign = (disp.dot(nd->second) > 0) ? 1.0 : -1.0;
					break;
				}
			}
			if (flag == 1)
				break;

			++counter;
		}

		return std::make_pair(idx, sign);
	}




};


// Morteza: use this to override symmetry orbits.
template <class EMesh>
class WireMesh2DMorteza : public WireMesh2D<EMesh>
{

public:
	typedef typename EMesh::VertexType      VertexType;
	typedef typename EMesh::EdgeType        EdgeType;
	typedef typename VertexType::ScalarType ScalarType;
	typedef typename VertexType::CoordType  CoordType;

	WireMesh2DMorteza(void) {;}

	WireMesh2DMorteza(EMesh & em)
	{
		this->setMesh(em);
	}

	WireMesh2DMorteza(const std::string & edgeMeshPath)
	{
		this->setMesh(edgeMeshPath);
	}

	// MHS on JUL14, 2015:
	// A new constructor with inSymmetryMode as the input ...
	WireMesh2DMorteza(const std::string & edgeMeshPath, const int inSymmetryMode)
	{
		this->m_symmetry_mode = inSymmetryMode;
		this->setMesh(edgeMeshPath);
	}


	WireMesh2DMorteza & operator = (WireMesh2DMorteza & wm)
	{
		m_operations      = wm.m_operations;
		m_symmetry_orbits = wm.m_symmetry_orbits;
		m_orbits_idx      = wm.m_orbits_idx;
		m_radius_range    = wm.m_radius_range;
		m_transl_range    = wm.m_transl_range;
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;

		return *this;
	}

	WireMesh2DMorteza(const WireMesh2DMorteza & wm)
	    : m_operations(wm.m_operations)
	    , m_symmetry_orbits(wm.m_symmetry_orbits)
	    , m_orbits_idx(wm.m_orbits_idx)
	    , m_radius_range(wm.m_radius_range)
	    , m_transl_range(wm.m_transl_range)
	{
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em,            wm.m_em,            false, true);
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, wm.m_normalized_em, false, true);
		m_bbox = wm.m_bbox;
	}

	void setMesh(const std::string & edgeMeshPath)
	{
		m_operations.clear();
		bool done = EdgeMeshUtils<EMesh>::importObj(m_em, edgeMeshPath);

		if (!done)
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
			return;
		}

		this->setup();
	}

	void setMesh(EMesh & em)
	{
		m_operations.clear();
		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_em, em, false, true);
		this->setup();
	}

	bool isValid(void) const
	{
		return (m_em.VN() > 0) && (m_em.EN() > 0);
	}

	size_t numberOfParameters(void) const
	{
		return m_operations.size();
	}

	CellParameters createCellParameters(void) const
	{
		CellParameters params(this->numberOfParameters());
		for (size_t i=0; i<params.numberOfParameters(); ++i)
		{
			std::pair<double,double> range = this->getParameterRange(i);
			params.parameter(i) = (range.first + range.second) / 2;
		}

		return params;
	}

	/*!
	 * \brief parametersValid checks wheter the current set of parameters is valid for the current wire mesh.
	 * \return true if the parameters are valid, false otherwise.
	 */
	bool parametersValid(const CellParameters & params) const
	{
		if (params.numberOfParameters() != m_operations.size())
			return false;

		std::pair<double,double> range_radius = this->getRadiusParameterRange();
		std::pair<double,double> range_transl = this->getTranslationParameterRange();
		for (size_t i=0; i<m_operations.size(); ++i)
		{
			switch (m_operations.at(i).type)
			{
			case ParameterOperation::Radius:
				if (params.cParameter(i) < range_radius.first ||
				    params.cParameter(i) > range_radius.second)
					return false;
				break;
			case ParameterOperation::Translation:
				if (params.cParameter(i) < range_transl.first ||
				    params.cParameter(i) > range_transl.second)
					return false;
				break;
			default:
				assert(0);
			}
		}
		return true;
	}

	std::pair<double,double> getParameterRange(int i) const
	{
		assert(i>=0 && i<int(m_operations.size()));

		switch (m_operations.at(i).type)
		{
		case ParameterOperation::Radius:
			return this->getRadiusParameterRange();
		case ParameterOperation::Translation:
			return this->getTranslationParameterRange();
		default:
			assert(0);
		}
		return std::make_pair(0,0);
	}

	void getUnmodifiedNormalizedEdgeMesh(EMesh & em)
	{
		if (!this->isValid())
			return;

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(em, m_normalized_em, false, true);
	}

	void getUnmodifiedEdgeMesh(EMesh & em)
	{
		if (!this->isValid())
			return;

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(em, m_em, false, true);
	}

	void getEdgeMesh(EMesh & em, const CellParameters & params)
	{
		if (!this->isValid())
			return;

		assert(this->parametersValid(params));

		this->getUnmodifiedEdgeMesh(em);
		this->applyParameterOperations(em, params);
	}

	void getNormalizedMesh(EMesh & em, const CellParameters & params)
	{
		if (!this->isValid())
			return;

		assert(this->parametersValid(params));

		this->getUnmodifiedNormalizedEdgeMesh(em);
		this->applyParameterOperations(em, params);
	}

	vcg::Point2<ScalarType> getOriginalScale(void) const
	{
		return vcg::Point2<ScalarType>(m_bbox.DimX(), m_bbox.DimY());
	}

	const std::vector<ParameterOperation> & getParameterOperations(void) const
	{
		// all symmetry orbit radius parameters come first
		return m_operations;
	}

	int getOrbitIndexForNode(int node)
	{
		if (!this->isValid() || node < 0 || node >= m_em.VN())
			return -1;

		return metaEMesh<EMesh>::findOperationIdx(node, m_operations);
		/* auto tag = vcg::tri::Allocator<EMesh>::template FindPerVertexAttribute<int>(m_em, SymmetryOrbitAttributeName()); */
		/* if (!vcg::tri::Allocator<EMesh>::IsValidHandle(m_em, tag)) */
		/* { */
		/* 	std::cout << "THIS WAS CALLED!" << std::endl; */
		/* 	return -1; */
		/* } */

		/* return m_orbits_idx[tag[node]]; */
	}

	static std::vector<double> generateNewParameters(const int currentSymmetry, const int targetSymmetry, std::vector<double> inParams, const std::string & edgeMeshPath)
	{

		// read in the mesh and do the initialization stuff see setMesh() and setup() above
		EMesh em;
		bool done = EdgeMeshUtils<EMesh>::importObj(em, edgeMeshPath);

		if (!done)
			throw("Could not read the input Mesh!");

		bool isMeshValid = (em.VN() > 0) && (em.EN() > 0);

		if (!isMeshValid)
			throw("Input Mesh is not valid!");

		vcg::tri::Allocator<EMesh>::CompactEveryVector(em);
		vcg::tri::UpdateTopology<EMesh>::VertexEdge(em);
		vcg::tri::UpdateTopology<EMesh>::EdgeEdge(em);


		std::vector<double> outParams;
		std::vector<ParameterOperation> currentOperations = generateParameterOperations(currentSymmetry, em);

		std::vector<std::pair<int, bool>> cMap = conversionMap(currentSymmetry, targetSymmetry, em);

		if (inParams.size() != currentOperations.size())
			throw("the number of input Parameters does not match the corrent symmetry mode!");

		for (int i = 0; i < cMap.size(); ++i)
			outParams.push_back(inParams[std::get<0>(cMap[i])] * (2.0 * double(std::get<1>(cMap[i])) - 1.0) );
		return outParams;
	}

    virtual ~WireMesh2DMorteza() { }

protected:
	EMesh                           m_em;
	EMesh                           m_normalized_em;
	vcg::Box3<ScalarType>           m_bbox; // the original bounding box
	std::vector<ParameterOperation> m_operations;
	std::vector<std::vector<int> >  m_symmetry_orbits;
	std::unordered_map<int,int>     m_orbits_idx;
	std::pair<double,double>        m_radius_range;
	std::pair<double,double>        m_transl_range;
	int 							m_symmetry_mode;

	static inline ScalarType Epsilon(void) { return 0.00001; }
	static inline ScalarType OneForth(void) { return 0.25; }
	static inline ScalarType ThreeForth(void) { return 0.75; }

	void setup(void)
	{
		this->m_radius_range = std::make_pair( 0.01, 5.0 );
		this->m_transl_range = std::make_pair(-0.18, 0.18);

		if (!this->isValid())
		{
			m_em.Clear();
			m_normalized_em.Clear();
			m_bbox.SetNull();
			return;
		}

		vcg::tri::Allocator<EMesh>::CompactEveryVector(m_em);
		vcg::tri::UpdateTopology<EMesh>::VertexEdge(m_em);
		vcg::tri::UpdateTopology<EMesh>::EdgeEdge(m_em);

		this->normalizeMesh();

		this->setParameterOperations();
		/* displayParameterOperations(m_operations); */
	}

	virtual std::pair<double,double> getRadiusParameterRange(void) const
	{
		return m_radius_range;
	}

	virtual std::pair<double,double> getTranslationParameterRange(void) const
	{
		return m_transl_range;
	}

	static std::string SymmetryOrbitAttributeName(void)
	{
		return "symmetry_orbit";
	}

	void normalizeMesh(void)
	{
		if (!this->isValid())
			return;

		// store bounding box for later use
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);
		m_bbox = m_em.bbox;

		vcg::tri::UpdatePosition<EMesh>::Translate(m_em, -m_em.bbox.min);
		vcg::tri::UpdateBounding<EMesh>::Box(m_em);

		assert(m_em.bbox.Dim()[0] != 0 && m_em.bbox.Dim()[1] != 0);

		vcg::tri::Append<EMesh,EMesh>::MeshCopy(m_normalized_em, m_em, false, true);

		// real normalization
		for (size_t i=0; i<m_normalized_em.vert.size(); ++i)
		{
			CoordType & p = m_normalized_em.vert[i].P();
			p[0] /= m_normalized_em.bbox.max[0];
			p[1] /= m_normalized_em.bbox.max[1];
			p[2] = 0;
		}
		vcg::tri::UpdateBounding<EMesh>::Box(m_normalized_em);
	}


	// MHS on AUG01, 2015:
	// this a much neater version using the new metaEMesh class where each case has more symmetry than the previous case
	// note that case 3 is equivalent to the Luigi's case (even though the ordering of the parameters might be different, like Luigi's case all the radius ops come first)
	static std::vector<ParameterOperation> generateParameterOperations(const int inSymmetryMode, EMesh & em)
	{
		// create the metaEMesh form em
		metaEMesh<EMesh> metaEM(em);

		// apply the symmetries based on the inSymmetryMode
		switch (inSymmetryMode)
		{
			case 0: // no symmetry
				metaEM.addBaseSymmetry();
				break;
			case 1: // symmetry about y = 0 (or x-axis) only
				metaEM.addUpDownSymmetry();
				metaEM.addBaseSymmetry();
				break;
			case 2: // symmetry about x = 0 (or y-axis) only
				metaEM.addLeftRightSymmetry();
				metaEM.addBaseSymmetry();
				break;
			case 3: // symmetry about x = 0 and y = 0 (this is Luigi's symmetry)
				metaEM.addUpDownSymmetry();
				metaEM.addLeftRightSymmetry();
				metaEM.addBaseSymmetry();
				break;
			case 4: // symmetry about x = y
				metaEM.addPosDiagSymmetry();
				metaEM.addBaseSymmetry();
				break;
			case 5: // symmetry about x = -y
				metaEM.addNegDiagSymmetry();
				metaEM.addBaseSymmetry();
				break;
			case 6: // symmetry about x = y and x = -y
				metaEM.addPosDiagSymmetry();
				metaEM.addNegDiagSymmetry();
				metaEM.addBaseSymmetry();
				break;
			case 7: // symmetry about x = 0, y = 0, x = y, and x = -y
				metaEM.addUpDownSymmetry();
				metaEM.addLeftRightSymmetry();
				metaEM.addPosDiagSymmetry();
				metaEM.addNegDiagSymmetry();
				metaEM.addBaseSymmetry();
				break;
			// add other cases here
		}
		return metaEM.getOperations(); // radius ops first and then displacements
	}


	// MHS on JUL24, 2015:
	// this method displays the parameter operations
	static void displayParameterOperations(const std::vector<ParameterOperation> & operations)
	{
		int counter = 0;
		std::cout << std::endl <<  "----------" << std::endl << "printing out the parameter operations" << std::endl << "----------" << std::endl;
		for (auto& op : operations)
		{
			std::cout << "operation " << counter << "'s type is " << ((op.type == 1) ? "radius" : "translation") << std::endl;
			if (op.type == ParameterOperation::Radius){
				std::cout << "nodes:" << std::endl;
				for (auto n : op.nodes)
					std::cout << n << "\t";
				std::cout << std::endl;
			}
			else {
				std::cout << "nodes_displ:" << std::endl;
				for (auto nd = op.nodes_displ.begin() ; nd != op.nodes_displ.end(); ++nd)
					std::cout << "node " << nd->first << "'s displ is: (" << nd->second[0] << "," << nd->second[1] << ")\t";
				std::cout << std::endl;
			}
			++counter;
		}
		std::cout << std::endl;
	}

	void setParameterOperations()
	{
		std::vector<ParameterOperation> operations = generateParameterOperations(m_symmetry_mode, m_em);
		m_operations = operations;
	}

	// MHS on JUL 23, 2015:
	// the method generates a conversion map from parameters of currentSymmetry mode to parameters of targetSymmetry mode
	// the map is a vector of <int, bool> pairs where
	// pair->first  points to the parameter in the current parameter list
	// pair->second specifies the direction 0: same 1: oposite
	static std::vector<std::pair<int, bool>> conversionMap(const int currentSymmetry, const int targetSymmetry, EMesh & em)
	{
		std::vector<ParameterOperation> currentOperations, targetOperations;
		std::vector<std::pair<int, bool>> conversion;
		currentOperations = generateParameterOperations(currentSymmetry, em);
		targetOperations  = generateParameterOperations(targetSymmetry,  em);

		int counter;
		bool direction;
		conversion.clear();
		for (auto& top : targetOperations)
		{
			counter = 0;
			direction = true;
			for (auto& cop : currentOperations)
			{
				if (top.type ==	cop.type)
				{
					if (top.type == ParameterOperation::Radius)
					{
						int targetNode = top.nodes[0];
						int flag = 0;
						for (int i = 0; i < cop.nodes.size(); ++i)
						{
							if (cop.nodes[i] == targetNode)
							{
								flag = 1;
								break;
							}
						}
						if (flag == 1)
							break;
						else
							++counter;
					}
					else
					{
						int targetNode = top.nodes_displ.begin()->first;
						vcg::Point2d targetDispl = top.nodes_displ.begin()->second;
						int flag = 0;
						for (auto it = cop.nodes_displ.begin(); it != cop.nodes_displ.end(); ++it)
						{
							if (it->first == targetNode &&
								vcg::math::Abs(targetDispl * it->second) > Epsilon())
							{
								flag = 1;
								direction = targetDispl * it->second > 0.0;
								break;
							}
						}
						if (flag == 1)
							break;
						else
							++counter;
					}
				}
				else
				{
					++counter;
				}
			}
			conversion.push_back(std::make_pair(counter, direction));
		}
		return conversion;
	}




	void applyParameterOperations(EMesh & em, const CellParameters & params)
	{
		assert(m_operations.size() == params.numberOfParameters());

		for (size_t i=0; i<m_operations.size(); ++i)
		{
			const ParameterOperation & p = m_operations[i];
			switch (p.type)
			{
			// store vertex radius into quality
			case ParameterOperation::Radius:
				for (int node : p.nodes)
				{
					em.vert[node].Q() = params.cParameter(i);
				}
				break;
			// displace the vertex
			case ParameterOperation::Translation:
				for (auto & displ : p.nodes_displ)
				{
					em.vert[displ.first].P() +=
					        CoordType(displ.second[0] * em.bbox.Dim()[0], displ.second[1] * em.bbox.Dim()[1], 0) * params.cParameter(i);
				}
				break;
			default:
				assert(0);
			}
		}
	}

	template <typename T>
	static int sign(T x)
	{
		return (T(0) < x) - (x < T(0));
	}
};

#endif // WIREMESH2D_H
