////////////////////////////////////////////////////////////////////////////////
// ObjectiveTermNormalizations.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Organizes the normalizations for each term in a composite pattern
//  optimization objective.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  04/07/2016 22:07:09
////////////////////////////////////////////////////////////////////////////////
#ifndef OBJECTIVETERMNORMALIZATIONS_HH
#define OBJECTIVETERMNORMALIZATIONS_HH

#include <map>
#include <string>
#include <stdexcept>

namespace PatternOptimization {

struct ObjectiveTermNormalizations {
    Real operator[](const std::string &name) const { return m_normalizations.at(name); }

    void set(const std::string &name, Real val) {
        if (m_normalizations.count(name)) throw std::runtime_error(name + " normalization already set.");
        m_normalizations.emplace(name, val);
    }

    bool isSet(const std::string &name) const { return m_normalizations.count(name); }

private:
    std::map<std::string, Real> m_normalizations;
};

}

#endif /* end of include guard: OBJECTIVETERMNORMALIZATIONS_HH */
