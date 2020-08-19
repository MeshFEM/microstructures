////////////////////////////////////////////////////////////////////////////////
// AutomaticDifferentiation.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Includes and functions needed for automatic differentiation
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/23/2015 15:26:43
////////////////////////////////////////////////////////////////////////////////
#ifndef MICROAUTOMATICDIFFERENTIATION_HH
#define MICROAUTOMATICDIFFERENTIATION_HH

#include <MeshFEM/AutomaticDifferentiation.hh>
#include <unsupported/Eigen/AutoDiff>
#include <cmath>
#include <algorithm>
#include <limits>
#include <type_traits>

// using ADReal = adept::adouble;

template<typename T> struct IsAutoDiffType : public std::false_type { };
template<typename T> struct IsAutoDiffType<Eigen::AutoDiffScalar<T>> : public std::true_type { };
// template<>           struct IsAutoDiffType<          adept::adouble> : public std::true_type { };

// GetADTypeOfPair<T1, T2>
// Metafunction to return either T1 or T2 depending on which is an autodiff
// type. If both or neither is an autodiff type, we return T1, which is assumed
// to equal T2.
template<typename T1, typename T2, bool T1AD, bool T2AD>
struct GetADTypeOfPairImpl { using type = T1; static_assert(std::is_same<T1, T2>::value, "Types must equal if neither or both are autodiff."); };
template<typename T1, typename T2> struct GetADTypeOfPairImpl<T1, T2,  true, false> { using type = T1; };
template<typename T1, typename T2> struct GetADTypeOfPairImpl<T1, T2, false,  true> { using type = T2; };

template<typename T1, typename T2>
struct GetADTypeOfPair : public GetADTypeOfPairImpl<T1, T2, IsAutoDiffType<T1>::value, IsAutoDiffType<T2>::value> { };

#endif /* end of include guard: MICROAUTOMATICDIFFERENTIATION_HH */
