////////////////////////////////////////////////////////////////////////////////
// GraphLineMesher.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Mesher of the 1D pattern "inflation graph" (intended for debugging).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  11/14/2015 12:44:30
////////////////////////////////////////////////////////////////////////////////
#ifndef GRAPHLINEMESHER_HH
#define GRAPHLINEMESHER_HH

#include "MeshingOptions.hh"
#include <MeshFEM/MeshIO.hh>

template<class SignedDistanceFunction>
class GraphLineMesher {
public:
    typedef typename SignedDistanceFunction::Real Real;
    GraphLineMesher(const MeshingOptions &opts = MeshingOptions())
        : meshingOptions(opts) { }

    void mesh(const SignedDistanceFunction &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements);

    MeshingOptions meshingOptions;
};

#endif /* end of include guard: GRAPHLINEMESHER_HH */
