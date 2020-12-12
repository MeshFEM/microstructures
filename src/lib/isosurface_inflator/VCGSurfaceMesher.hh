////////////////////////////////////////////////////////////////////////////////
// VCGSurfaceMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      VCG-based surface mesher (used for quick previews of the object).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  11/11/2015 20:29:11
////////////////////////////////////////////////////////////////////////////////
#ifndef VCGSURFACEMESHER_HH
#define VCGSURFACEMESHER_HH

#include "MesherBase.hh"
#include "SignedDistanceRegion.hh"
#include <MeshFEM/MeshIO.hh>

class VCGSurfaceMesher : public MesherBase {
public:
    using Real = SignedDistanceRegion<3>::Real;
    using MesherBase::MesherBase;
    using MesherBase::meshingOptions;

    virtual void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements) const override;
};

#endif /* end of include guard: VCGSURFACEMESHER_HH */
