////////////////////////////////////////////////////////////////////////////////
// CGALClippedVolumeMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      CGAL-based volume mesher for meshing a smooth signed distance function
//      intersected with an axis-aligned box; assumes that the only sharp
//      features come from the boolean intersection, and those sharp feature
//      curves are resolved using marching squares on the box faces.
//
//      The SignedDistanceFunction template parameter should be a class
//      provding the following functions:
//          static BBox<Point3<Real>> boundingBox()
//          template<R> signedDistance(Point3<R> p) const
//          template<R> boundingSphere(Point3<R> &c, R &r) const
//      boundingBox() determines the clipping region.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/29/2015 15:42:57
////////////////////////////////////////////////////////////////////////////////
#ifndef CGALCLIPPEDVOLUMEMESHER_HH
#define CGALCLIPPEDVOLUMEMESHER_HH

#include "MesherBase.hh"
#include <MeshFEM/MeshIO.hh>

#include "SignedDistanceRegion.hh"

class CGALClippedVolumeMesher : public MesherBase {
public:
    using Real = typename SignedDistanceRegion<3>::Real;
    using MesherBase::MesherBase;
    using MesherBase::meshingOptions;

    virtual ~CGALClippedVolumeMesher() = default;

    virtual void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements) const override;
private:
    struct ClippedSignedDistanceFunction;
};

#endif /* end of include guard: CGALCLIPPEDVOLUMEMESHER_HH */
