////////////////////////////////////////////////////////////////////////////////
// SignedDistanceRegion.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Abstract base class defining the interface to a Nd region to be meshed
//      using a signed distance/inside-outside query. Also implements slicing
//      3D regions into 2D regions.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/31/2017 21:43:22
////////////////////////////////////////////////////////////////////////////////
#ifndef SIGNEDDISTANCEREGION_HH
#define SIGNEDDISTANCEREGION_HH

#include "InflatorTypes.hh"

template<size_t Dim>
struct SignedDistanceRegion {
    using Real = double;
    virtual const BBox<PointNd<Dim>> &boundingBox() const = 0;
    virtual double signedDistance(const PointNd<Dim> &p) const = 0;

    // Unless implemented more efficiently in derived class, implement
    // inside/outside test using the signed distance query
    virtual bool isInside(const PointNd<Dim> &p) const { return signedDistance(p) <= 1.0; }

    // When the CGAL mesher is used, the boundingSphere must contain the
    // geometry and be centered inside the object. This method should be
    // overriden by regions intended to be used for meshing with CGAL.
    virtual void boundingSphere(PointNd<Dim> &/* c */, double &/* r */) const {
        throw std::runtime_error("Unimplemented");
    }

    virtual ~SignedDistanceRegion() { }
};

#endif /* end of include guard: SIGNEDDISTANCEREGION_HH */
