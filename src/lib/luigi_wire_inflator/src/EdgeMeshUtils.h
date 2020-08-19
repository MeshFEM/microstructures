#ifndef EDGEMESHUTILS_H
#define EDGEMESHUTILS_H

#include <vcg/complex/algorithms/update/bounding.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

template <class MESH>
class EdgeMeshUtils
{
public:
	static void cleanMesh(MESH & mesh)
	{
		vcg::tri::Clean<MESH>::RemoveUnreferencedVertex(mesh);
		vcg::tri::Clean<MESH>::RemoveDuplicateVertex(mesh);
	}

	static void compactMesh(MESH & mesh)
	{
		vcg::tri::Allocator<MESH>::CompactEveryVector(mesh);
	}

	///
	/// \brief exportObj
	/// \param mesh the mesh to be exported.
	/// \param filePath the file path to which export to.
	/// \param export_attrib if false it does not export per vertex attributes (i.e. exports geometry only).
	/// \return true if the export was successful, false otherwise.
	///
	static bool exportObj(MESH & mesh, const std::string & filePath, bool export_attrib = true)
	{
		typedef typename vcg::tri::io::ExporterOBJ<MESH> Exporter;

		int cap = Exporter::GetExportMaskCapability();
		int mask = export_attrib ? createExportMask(cap, mesh) : 0;

		int res = Exporter::Save(mesh, filePath.c_str(), mask);

		return res == Exporter::E_NOERROR;
	}

	///
	/// \brief importObj
	/// \param mesh the mesh where to import.
	/// \param filePath the file path from which import to.
	/// \return true if the import was successful, false otherwise.
	///
	static bool importObj(MESH & mesh, const std::string & filePath)
	{
		typedef typename vcg::tri::io::ImporterOBJ<MESH> Importer;

		int load_mask;

		if (!Importer::LoadMask(filePath.c_str(), load_mask))
			return false;

		int result = Importer::Open(mesh, filePath.c_str(), load_mask);
		if (Importer::ErrorCritical(result))
			return false;

		vcg::tri::UpdateBounding<MESH>::Box(mesh);

		return true;
	}

	static int createExportMask(const int exportCapabilities, const MESH & mesh)
	{
		int mask = vcg::tri::io::Mask::IOM_NONE;

		if (vcg::tri::HasPerVertexNormal   (mesh)) mask |= (vcg::tri::io::Mask::IOM_VERTNORMAL & exportCapabilities);

		if (vcg::tri::HasPerVertexColor    (mesh)) mask |= (vcg::tri::io::Mask::IOM_VERTCOLOR & exportCapabilities);

		if (vcg::tri::HasPerVertexTexCoord (mesh)) mask |= (vcg::tri::io::Mask::IOM_VERTTEXCOORD & exportCapabilities);

		return mask;
	}
};


#endif // EDGEMESHUTILS_H
