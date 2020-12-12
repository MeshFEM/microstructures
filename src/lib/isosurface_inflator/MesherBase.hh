////////////////////////////////////////////////////////////////////////////////
// MesherBase.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Abstract base class for all isosurface inflator meshers.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/06/2017 02:27:47
////////////////////////////////////////////////////////////////////////////////
#ifndef MESHERBASE_HH
#define MESHERBASE_HH

#include "MeshingOptions.hh"
#include "SignedDistanceRegion.hh"
#include <MeshFEM/MeshIO.hh>

class MesherBase {
public:
    virtual void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements) const = 0;

    MeshingOptions meshingOptions;

    virtual ~MesherBase() = default;

    // Whether the mesher should attempt to make the boundary mesh lying on the
    // symmetry cell interface depend only on the interface geometry, allowing
    // meshes to be stitched together when the geometry matches.
    // This is only supported by the MidplaneMesher, and for patterns with
    // orthotropic symmetry it must be used with
    // meshingOptions::curveSimplifier=RESAMPLE (the default).
    // This mode should only be enabled when actually needed since it prevents
    // Triangle from inserting Steiner points.
    bool meshInterfaceConsistently = false;
};

#endif /* end of include guard: MESHERBASE_HH */
