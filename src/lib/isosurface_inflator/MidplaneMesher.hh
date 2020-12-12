////////////////////////////////////////////////////////////////////////////////
// MidplaneMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Mesh midplane slice through the 3D signed distance function. I.e. mesh
//      the 1D intersection of the 3D pattern surface with the z = 0 midplane.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  02/03/2016 11:41:48
////////////////////////////////////////////////////////////////////////////////
#ifndef MIDPLANEMESHER_HH
#define MIDPLANEMESHER_HH

#include "MesherBase.hh"
#include "SignedDistanceRegion.hh"
#include <MeshFEM/MeshIO.hh>
#include <string>

class MidplaneMesher : public MesherBase {
public:
    using MesherBase::MesherBase;
    using MesherBase::meshingOptions;

    virtual void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements) const override;
    virtual void meshSlice(const SignedDistanceRegion<2> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements) const;
    void dumpSDF(const SignedDistanceRegion<3> &sdf, const std::string &path);

    std::string msPolygonPath;
};

#endif /* end of include guard: MIDPLANEMESHER_HH */
