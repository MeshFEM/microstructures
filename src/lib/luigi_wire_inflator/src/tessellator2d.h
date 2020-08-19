#ifndef TESSELLATOR2D_H
#define TESSELLATOR2D_H

#include "clipperHelper.h"

extern "C" {
#define ANSI_DECLARATORS
#define REAL double
#define VOID int
#include <triangle.h>
}

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/clean.h>

#include <iostream>
#include <cstring>
#include <iomanip>
#include <unordered_map>

struct Tessellator2DSettings
{
	bool plsg;
	bool area_constrained;
	double max_area;
	bool angle_constrained;
	double min_angle;
	bool steiner_points;
	bool quiet;

	Tessellator2DSettings(void)
	    : plsg             (true)
	    , area_constrained (true)
	    , max_area         (0.001)
	    , angle_constrained(true)
	    , min_angle        (30)
	    , steiner_points   (true)
	    , quiet            (false)
	{
		;
	}

	std::string getTriangleSwitches(void) const
	{
		// z  -> uses 0 based index
		// p  -> triangulates a Planar Straight Line Graph
		// a  -> set maximum triangle area
		// q  -> quality tessellation (min angle constraint)
		// Y  -> no steiner points (AVOID)
		// YY -> steiner point added except on boundary segments
		// Q  -> quiet

		std::ostringstream ss;
		ss << 'z';
		if (plsg) ss << 'p';
		if (quiet) ss << 'Q';
		if (!steiner_points) ss << "YY";
		if (angle_constrained)
		{
			ss << 'q';
			if (min_angle > 0)
				ss << std::fixed << std::setprecision(9) << min_angle;
		}
		if (area_constrained)
		{
			ss << 'a';
			if (max_area > 0)
				ss << std::fixed << std::setprecision(9) << max_area;
		}
		return ss.str();
	}

	bool isValid() const
	{
		return (!area_constrained  || max_area > 0)
		    && (!angle_constrained || (min_angle > 0 && min_angle <= 33));
	}
};

template <class VcgMesh>
class Tessellator2D
{
public:
	typedef typename VcgMesh::ScalarType ScalarType;
	typedef typename VcgMesh::CoordType  CoordType;

	static bool execute(const Tessellator2DSettings & settings, const ClipperLib::Paths & polygons, const ScalarType clipperScaleFactor, VcgMesh & out_mesh)
	{
		// Do triangulation using only boundary segments (a Planar Straight Line Graph)
		if (polygons.size() == 0) return false;

		// create in and out structs for triangle
		triangulateio in, out;
		memset(&in , 0, sizeof(triangulateio));
		memset(&out, 0, sizeof(triangulateio));

		// count and allocate vertices and indices
		std::unordered_map<ClipperLib::IntPoint, size_t> v_idx;
		std::vector<ClipperLib::IntPoint>                vertices;
		size_t totalVertices = 0;

		for (size_t i=0; i<polygons.size(); ++i)
		{
			const ClipperLib::Path & path = polygons[i];
			totalVertices += path.size();
			for (size_t k=0; k<path.size(); ++k)
			{
				const ClipperLib::IntPoint & p = path[k];
				if (v_idx.count(p) == 0)
				{
					size_t idx = vertices.size();
					v_idx[p] = idx;
					vertices.push_back(p);
				}
			}
		}
		in.numberofpoints   = vertices.size();
		in.numberofsegments = totalVertices;

		// initialize lists
		in.pointlist         = (REAL *) malloc(in.numberofpoints   * 2 * sizeof(REAL));
		in.segmentlist       = (int *)  malloc(in.numberofsegments * 2 * sizeof(int));
		in.segmentmarkerlist = (int *)  malloc(in.numberofsegments * sizeof(int));

		// fill triangle input structure with points
		for (size_t i=0; i<vertices.size(); ++i)
		{
			const ClipperLib::IntPoint & p = vertices[i];
			in.pointlist[i*2]     = REAL(p.X) / clipperScaleFactor;
			in.pointlist[i*2 + 1] = REAL(p.Y) / clipperScaleFactor;
		}

		// fill triangle input structure with boundary segments
		size_t sidx = 0;
		for (size_t i=0; i<polygons.size(); ++i)
		{
			const ClipperLib::Path & path = polygons[i];
			for (size_t k=0; k<path.size(); ++k)
			{
				const ClipperLib::IntPoint & p      = path[k];
				const ClipperLib::IntPoint & next_p = path[(k+1)%path.size()];
				in.segmentlist[sidx*2]     = v_idx[p];
				in.segmentlist[sidx*2 + 1] = v_idx[next_p];
				in.segmentmarkerlist[sidx] = 1; // mark each segment as boundary
				sidx++;
			}
		}

		//find holes
		ClipperLib::Path holePoints = getPointsInHoles(polygons);
		if (holePoints.size() > 0)
		{
			in.numberofholes = holePoints.size();
			in.holelist = (REAL *) malloc(in.numberofholes * 2 * sizeof(REAL));
			for (int i=0; i<in.numberofholes; ++i)
			{
				in.holelist[i*2]     = REAL(holePoints[i].X) / clipperScaleFactor;
				in.holelist[i*2 + 1] = REAL(holePoints[i].Y) / clipperScaleFactor;
			}
		}

		// launch the triangulation with parameters from settings
		std::string cmd = settings.getTriangleSwitches();
		triangulate((char *)cmd.c_str(), &in, &out, (triangulateio *) NULL);

		// output check
		if (out.numberoftriangles <= 0)
		{
			std::cout << "unable to perform triangulation!" << std::endl << std::flush;
			deallocateTriangulation(in, out);
			return false;
		}

		if (!settings.quiet)
			std::cout << std::flush;

		// create the vcg mesh
		out_mesh.Clear();
		vcg::tri::Allocator<VcgMesh>::AddVertices(out_mesh, out.numberofpoints);
		vcg::tri::Allocator<VcgMesh>::AddFaces(out_mesh, out.numberoftriangles);
		for (int i=0; i<out.numberofpoints; i++)
		{
			// vertices
			out_mesh.vert[i].P() = CoordType(ScalarType(out.pointlist[i*2]), ScalarType(out.pointlist[(i*2)+1]), 0);
		}
		for (int i=0; i<out.numberoftriangles; i++)
		{
			// faces
			int index0 = out.trianglelist[i*3];
			int index1 = out.trianglelist[(i*3)+1];
			int index2 = out.trianglelist[(i*3)+2];
			out_mesh.face[i].V(0) = &out_mesh.vert[index0];
			out_mesh.face[i].V(1) = &out_mesh.vert[index1];
			out_mesh.face[i].V(2) = &out_mesh.vert[index2];
		}
		deallocateTriangulation(in, out);

		// Clean steps
		vcg::tri::Clean<VcgMesh>::RemoveUnreferencedVertex(out_mesh);
		vcg::tri::Clean<VcgMesh>::RemoveDuplicateVertex(out_mesh);
		vcg::tri::Allocator<VcgMesh>::CompactEveryVector(out_mesh);

		// Update normals and bbox
		vcg::tri::UpdateNormal<VcgMesh>::PerFaceNormalized(out_mesh);
		vcg::tri::UpdateNormal<VcgMesh>::PerVertexFromCurrentFaceNormal(out_mesh);
		vcg::tri::UpdateBounding<VcgMesh>::Box(out_mesh);

		return true;
	}

	static void deallocateTriangulation(triangulateio & in, triangulateio & out)
	{
		// deallocate the triangle library input
		if (in.edgelist)              trifree((VOID *)in.edgelist);
		if (in.edgemarkerlist)        trifree((VOID *)in.edgemarkerlist);
		if (in.holelist)              trifree((VOID *)in.holelist);
		if (in.neighborlist)          trifree((VOID *)in.neighborlist);
		if (in.normlist)              trifree((VOID *)in.normlist);
		if (in.pointattributelist)    trifree((VOID *)in.pointattributelist);
		if (in.pointlist)             trifree((VOID *)in.pointlist);
		if (in.pointmarkerlist)       trifree((VOID *)in.pointmarkerlist);
		if (in.regionlist)            trifree((VOID *)in.regionlist);
		if (in.segmentlist)           trifree((VOID *)in.segmentlist);
		if (in.segmentmarkerlist)     trifree((VOID *)in.segmentmarkerlist);
		if (in.trianglearealist)      trifree((VOID *)in.trianglearealist);
		if (in.triangleattributelist) trifree((VOID *)in.triangleattributelist);
		if (in.trianglelist)          trifree((VOID *)in.trianglelist);

		// deallocate the triangle library output (this is unbelievable!!)
		if (out.edgelist              && (out.edgelist              != in.edgelist)             ) trifree((VOID *)out.edgelist);
		if (out.edgemarkerlist        && (out.edgemarkerlist        != in.edgemarkerlist)       ) trifree((VOID *)out.edgemarkerlist);
		if (out.holelist              && (out.holelist              != in.holelist)             ) trifree((VOID *)out.holelist);
		if (out.neighborlist          && (out.neighborlist          != in.neighborlist)         ) trifree((VOID *)out.neighborlist);
		if (out.normlist              && (out.normlist              != in.normlist)             ) trifree((VOID *)out.normlist);
		if (out.pointattributelist    && (out.pointattributelist    != in.pointattributelist)   ) trifree((VOID *)out.pointattributelist);
		if (out.pointlist             && (out.pointlist             != in.pointlist)            ) trifree((VOID *)out.pointlist);
		if (out.pointmarkerlist       && (out.pointmarkerlist       != in.pointmarkerlist)      ) trifree((VOID *)out.pointmarkerlist);
		if (out.regionlist            && (out.regionlist            != in.regionlist)           ) trifree((VOID *)out.regionlist);
		if (out.segmentlist           && (out.segmentlist           != in.segmentlist)          ) trifree((VOID *)out.segmentlist);
		if (out.segmentmarkerlist     && (out.segmentmarkerlist     != in.segmentmarkerlist)    ) trifree((VOID *)out.segmentmarkerlist);
		if (out.trianglearealist      && (out.trianglearealist      != in.trianglearealist)     ) trifree((VOID *)out.trianglearealist);
		if (out.triangleattributelist && (out.triangleattributelist != in.triangleattributelist)) trifree((VOID *)out.triangleattributelist);
		if (out.trianglelist          && (out.trianglelist          != in.trianglelist)         ) trifree((VOID *)out.trianglelist);
	}
};

#endif // TESSELLATOR2D_H
