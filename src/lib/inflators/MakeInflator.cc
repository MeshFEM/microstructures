#include "MakeInflator.hh"

#include "Inflator.hh"
#include "wrappers/IsoinflatorWrapper.hh"
#include "wrappers/LpHoleInflator.hh"
#include "wrappers/EqualityConstrainedInflator.hh"
#include "wrappers/BoundaryPerturbationInflator.hh"
#include "wrappers/JamesInflatorWrapper.hh"
#include "wrappers/LuigiInflatorWrapper.hh"

#include <MeshFEM/Future.hh>
#include <MeshFEM/Utilities/ci_string.hh>

#include <stdexcept>

using namespace std;
namespace po = boost::program_options;

template<typename T>
T extract_required(po::variables_map &opts, const string &key) {
    if (opts.count(key)) {
        T val(opts[key].as<T>());
        opts.erase(key);
        return val;
    }
    throw runtime_error("Options passed to makeInflator missing required option '" + key + "'");
}

template<typename T>
T extract_defaulted(po::variables_map &opts, const string &key, T default_val) {
    if (opts.count(key)) {
        T val = opts[key].as<T>();
        opts.erase(key);
        return val;
    }
    return default_val;
}

bool extract_flag(po::variables_map &opts, const string &key) {
    if (opts.count(key)) { opts.erase(key); return true; }
    return false;
}

template<typename T>
void extract_notify(po::variables_map &opts, const string &key, const T &func) {
    static_assert(std::is_same<typename function_traits<T>::result_type, void>::value, "Notify function must return void.");
    static_assert(function_traits<T>::arity == 1, "Notify function must take a single argument.");
    using OptionType = typename function_traits<T>::template arg<0>::type;
    if (opts.count(key)) {
        func(opts[key].as<OptionType>());
        opts.erase(key);
    }
}

unique_ptr<InflatorBase> make_inflator(const string &name, po::variables_map opts,
                                       const vector<string> &constraints) {
    unique_ptr<InflatorBase> infl;
    if (ci_string("Isosurface2D") == name.c_str()) {
        infl = Future::make_unique<IsoinflatorWrapper<2>>(
                extract_required<string>(opts, "pattern"),
                extract_required<std::string>(opts, "symmetry"),
                extract_flag(opts, "vertexThickness"),
                extract_defaulted<size_t>(opts, "inflation_graph_radius", 2));
    }
    else if (ci_string("Isosurface3D") == name.c_str()) {
        infl = Future::make_unique<IsoinflatorWrapper<3>>(
                extract_required<string>(opts, "pattern"),
                extract_required<std::string>(opts, "symmetry"),
                extract_flag(opts, "vertexThickness"),
                extract_defaulted<size_t>(opts, "inflation_graph_radius", 2));
    }
    else if (ci_string("James") == name.c_str()) {
        infl = Future::make_unique<JamesInflatorWrapper>(
                extract_required<string>(opts, "pattern"),
                extract_defaulted<double>(opts, "cell_size", 5.0),
                0.5 * sqrt(2),
                extract_required<std::string>(opts, "symmetry") == "cubic",
                extract_flag(opts, "vertexThickness"));
        infl->configureSubdivision(extract_defaulted<string>(opts, "sub_algorithm", "simple"),
                                   extract_defaulted<size_t>(opts,     "subdivide",        0));
    }
    else if (ci_string("Luigi") == name.c_str()) {
        infl = Future::make_unique<LuigiInflatorWrapper>(extract_required<string>(opts, "pattern"));
    }
    else if (ci_string("LpHole") == name.c_str()) {
        auto lphole_infl = Future::make_unique<LpHoleInflator>();
        extract_notify(opts, "hole_segments", [&](size_t hs) { lphole_infl->setNumSubdiv(hs); });
        infl = std::move(lphole_infl);
    }
    else if (ci_string("BoundaryPerturbation") == name.c_str()) {
        std::vector<MeshIO::IOVertex>  inVertices;
        std::vector<MeshIO::IOElement> inElements;
        auto type = MeshIO::load(extract_required<string>(opts, "pattern"),
                                 inVertices, inElements);

        if      (type == MeshIO::MESH_TRI) infl = Future::make_unique<BoundaryPerturbationInflator<2>>(inVertices, inElements, 1e-10);
        else if (type == MeshIO::MESH_TET) infl = Future::make_unique<BoundaryPerturbationInflator<3>>(inVertices, inElements, 1e-10);
        else    throw std::runtime_error("Mesh must be triangle or tet.");
    }
    else throw runtime_error("Invalid inflator: " + name);

    bool orthoCell = extract_flag(opts, "ortho_cell");
    bool  fullCell = extract_flag(opts, "fullCellInflator");

    if (orthoCell && fullCell) throw std::runtime_error("ortho_cell and fullCellInflator are mutually exclusive.");
    infl->setOrthoBaseCell(orthoCell);
    infl->setReflectiveInflator(!fullCell);

    if (fullCell) infl->setReflectiveInflator(false);

    // Dump to inflation_dump_path, if specified
    infl->setInflationDumpPath(extract_defaulted<string>(opts, "inflation_dump_path", ""));

    extract_notify(opts, "meshingOptions", [&](const string &v) { infl->loadMeshingOptions(v); });
    extract_notify(opts, "max_volume",     [&](double        v) { infl->setMaxElementVolume(v); });

    // Report any ignored options.
    for (const auto &entry : opts)
        cerr << "WARNING: inflator " << name << " ignored option " << entry.first << endl;

    // Wrap the inflator in an EqualityConstrainedInflator if there are
    // constraints.
    if (constraints.size()) {
        InflatorBase *rawPtr = infl.release();
        if (auto infl2D = dynamic_cast<Inflator<2> *>(rawPtr))
            infl = Future::make_unique<EqualityConstrainedInflator<2>>(unique_ptr<Inflator<2>>(infl2D), constraints);
        else if (auto infl3D = dynamic_cast<Inflator<3> *>(rawPtr))
            infl = Future::make_unique<EqualityConstrainedInflator<3>>(unique_ptr<Inflator<3>>(infl3D), constraints);
        else throw runtime_error("Invalid inflator dimension.");
    }

    return infl;
}

// Extract the options that are meant to be passed to the inflator constructor:
po_vm filterInflatorOptions(const po_vm &opts) {
    auto keys = {
        "pattern",
        "symmetry",
        "vertexThickness",
        "cell_size",
        "hole_segments",
        "max_volume",
        "meshingOptions",
        "subdivide",
        "sub_algorithm",
        "ortho_cell",
        "inflation_dump_path",
        "inflation_graph_radius",
        "params",
        "metaParams",
    };
    po_vm filtered;
    for (const string &key : keys)
        if (opts.count(key)) filtered.emplace(key, opts[key]);

    return filtered;
}
