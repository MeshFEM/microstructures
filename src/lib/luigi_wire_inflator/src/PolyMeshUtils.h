#ifndef POLYMESHUTILS_H
#define POLYMESHUTILS_H

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <string>
#include <iostream>

template <class MESH>
class PolyMeshUtils
{
public:
	static bool isQuadMesh(MESH& mesh)
	{
		for (auto it = mesh.face.begin(); it != mesh.face.end(); it++)
			if (it->VN() != 4)
				return false;

		return true;
	}

	static bool importFromOBJ(const std::string & file_name, MESH& mesh)
	{
		typedef vcg::tri::io::ImporterOBJ<MESH> PolyMeshImporterOBJ;

		int loadmask;

		// Clear the input polygonal mesh and load from file
		mesh.Clear();
		int res = PolyMeshImporterOBJ::Open(mesh, file_name.c_str(), loadmask);
		if (res != PolyMeshImporterOBJ::E_NOERROR)
		{
			std::cerr << "error loading file: \"" << file_name << "\":" << std::endl
			          << PolyMeshImporterOBJ::ErrorMsg(res) << std::endl;

			if (PolyMeshImporterOBJ::ErrorCritical(res))
				return false;
		}

		vcg::tri::UpdateTopology<MESH>::FaceFace(mesh);
		return true;
	}

	static bool exportToOBJ(const std::string & file_name, MESH & mesh)
	{
		typedef vcg::tri::io::ExporterOBJ<MESH> PolyMeshExporterOBJ;

		int mask = 0;

		int res = PolyMeshExporterOBJ::Save(mesh, file_name.c_str(), mask);
		if (res != PolyMeshExporterOBJ::E_NOERROR)
		{
			std::cerr << "error saving file: \"" << file_name << "\":" << std::endl
			          << PolyMeshExporterOBJ::ErrorMsg(res) << std::endl;

			return false;
		}
		return true;
	}
};

#endif // POLYMESHUTILS_H
