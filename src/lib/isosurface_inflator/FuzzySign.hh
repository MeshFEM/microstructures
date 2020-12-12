////////////////////////////////////////////////////////////////////////////////
// FuzzySign.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Implements sign queries with tolerance.
//      Note: with nonzero tolerance, isZero, isPositive, and isNegative all
//      return true for val = 0.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/20/2015 11:50:08
////////////////////////////////////////////////////////////////////////////////
#ifndef FUZZYSIGN_HH
#define FUZZYSIGN_HH

#include <ratio>
#include <cmath>

typedef std::ratio<1, long(1e12)> DEFAULT_TOL;
template<typename TOL = DEFAULT_TOL> bool isZero    (double val) { return std::abs(val) < (double(TOL::num) / double(TOL::den)); }
template<typename TOL = DEFAULT_TOL> bool isPositive(double val) { return val >          -(double(TOL::num) / double(TOL::den)); }
template<typename TOL = DEFAULT_TOL> bool isNegative(double val) { return val <           (double(TOL::num) / double(TOL::den)); }

#endif /* end of include guard: FUZZYSIGN_HH */
