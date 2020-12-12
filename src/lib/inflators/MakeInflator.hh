#ifndef MAKEINFLATOR_HH
#define MAKEINFLATOR_HH

#include <MeshFEM/Utilities/ci_string.hh>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

using po_vm = boost::program_options::variables_map;

////////////////////////////////////////////////////////////////////////////////
// Forward Declarations
////////////////////////////////////////////////////////////////////////////////
class InflatorBase;
template<size_t N> class Inflator;

// Construct and configure an inflator by name with the specified options and
// optional linear equality constraints.
// Names:
//      Isosurface[23]D: Isosurface-based 2 or 3D inflator
//      James:           James' 3D inflator
//      Luigi:           Luigi's 2D inflator
//      LpHole:          Punctured hole inflator
std::unique_ptr<InflatorBase> make_inflator(const std::string &name, po_vm options,
                                            const std::vector<std::string> &constraints = std::vector<std::string>());

// Construct and configure an dimension-specific inflator by name with the
// specified options and optional linear equality constraints.
// Names:
//      Isosurface:      Isosurface-based 2 or 3D inflator (inferred from N)
//      James:           James' 3D inflator
//      Luigi:           Luigi's 2D inflator
//      LpHole:          Punctured hole inflator
// Throws an error if the inflator specified by "name" is not of dimension N.
template<size_t N>
std::unique_ptr<Inflator<N>> make_inflator(std::string name, po_vm options,
                                           const std::vector<std::string> &constraints = std::vector<std::string>())
{
    if (ci_string("Isosurface") == name.c_str()) { name += std::to_string(N); name += "D"; }
    auto infl = make_inflator(name, options, constraints);

    InflatorBase *rawPtr = infl.release();
    if (auto inflND = dynamic_cast<Inflator<N> *>(rawPtr))
        return std::unique_ptr<Inflator<N>>(inflND);
    else throw std::runtime_error("Invalid inflator dimension.");
}

// Extract the options that are meant to be passed to the inflator's
// constructor/configuration:
//      "pattern"
//      "isotropicParameters"
//      "vertexThickness"
//      "cell_size"
//      "hole_segments"
//      "max_volume"
//      "meshingOptions"
//      "subdivide"
//      "sub_algorithm"
//      "ortho_cell"
//      "inflation_dump_path"
po_vm filterInflatorOptions(const po_vm &options);

#endif /* end of include guard: MAKEINFLATOR_HH */
