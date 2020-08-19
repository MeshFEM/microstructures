#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/Symmetry.hh>
#include <nlohmann/json.hpp>
////////////////////////////////////////////////////////////////////////////////

using json = nlohmann::json;

////////////////////////////////////////////////////////////////////////////////

template<typename Symmetry> struct SymmetryTraits { };

////////////////////////////////////////////////////////////////////////////////

#ifndef SYMMETRY_NAME
#define SYMMETRY_NAME(sym, name)                  \
    template<>                                    \
    struct SymmetryTraits<sym> {                  \
        static constexpr char value[] = name;     \
    };
#endif

SYMMETRY_NAME(Symmetry::Cubic<>, "cubic")
SYMMETRY_NAME(Symmetry::Orthotropic<>, "orthotropic")
SYMMETRY_NAME(Symmetry::Diagonal<>, "diagonal")
SYMMETRY_NAME(Symmetry::Square<>, "square")
SYMMETRY_NAME(Symmetry::TriplyPeriodic<>, "triply_periodic")
SYMMETRY_NAME(Symmetry::DoublyPeriodic<>, "doubly_periodic")
SYMMETRY_NAME(Symmetry::NonPeriodic<>, "non_periodic")
