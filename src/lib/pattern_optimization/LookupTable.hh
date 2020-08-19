////////////////////////////////////////////////////////////////////////////////
// LookupTable.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Reads, writes, and processes a pattern lookup table.
//      Format:
//      pattern_num\tE\tnu\tanisotropy\tnum_params\tparam_0....
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/30/2015 18:50:10
////////////////////////////////////////////////////////////////////////////////
#ifndef LOOKUPTABLE_HH
#define LOOKUPTABLE_HH

#include <MeshFEM/ElasticityTensor.hh>
#include <Eigen/Dense>

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>
#include <sstream>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
        elems.push_back(item);
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

template<typename Real>
struct IsotropicLookupTable {
    IsotropicLookupTable() { }
    IsotropicLookupTable(const std::string &path) { load(path); }

    void load(const std::string &path) {
        std::ifstream lutStream(path);
        if (!lutStream.is_open())
            throw std::runtime_error("Couldn't open " + path);

        std::string line;
        while (std::getline(lutStream, line)) {
            records.emplace_back(line);
        }
    }

    void write(const std::string &path) const {
        std::ofstream outFile(path);
        if (!outFile.is_open())
            throw std::runtime_error("Couldn't open " + path + " for writing");
        write(outFile);
    }

    void write(std::ostream &os) const {
        for (auto r : records)
            os << r << std::endl;
    }

    struct Record {
        Record(const std::string &lutLine) {
            std::runtime_error invalidColSize("Invalid number of columns in LUT.");
            auto tokens = split(lutLine, '\t');
            if (tokens.size() < 5) throw invalidColSize;
            size_t nParams = std::stoi(tokens[4]);
            if (tokens.size() != 5 + nParams) throw invalidColSize;
            pattern = std::stoi(tokens[0]);
            E = std::stod(tokens[1]);
            nu = std::stod(tokens[2]);
            anisotropy = std::stod(tokens[3]);

            for (size_t i = 0; i < nParams; ++i)
                patternParams.push_back(std::stod(tokens[5 + i]));
        }

        size_t numParams() const { return patternParams.size(); }
        size_t pattern;
        Real E, nu, anisotropy;

        template<size_t Dim>
        ElasticityTensor<Real, Dim> tensor() const { return ElasticityTensor<Real, Dim>(E, nu); }

        std::vector<Real> patternParams;

        // Write record without trailing newline.
        friend std::ostream &operator<<(std::ostream &os, const Record &r) {
            std::ios state(NULL);
            state.copyfmt(os);
            os << std::setfill('0') << std::setw(4) << r.pattern << '\t';
            os.copyfmt(state);
            os << std::setprecision(16);
            os << r.E << '\t' << r.nu << '\t' << r.anisotropy << '\t' << r.numParams();
            for (auto p : r.patternParams)
                os << '\t' << p;
            os.copyfmt(state);
            return os;
        }
    };

    // Add the patterns in this collection to set pats
    void patterns(std::set<size_t> &pats) const {
        for (const auto &r : records) pats.insert(r.pattern);
    }

    std::set<size_t> patterns() const {
        std::set<size_t> result;
        patterns(result);
        return result;
    }

    size_t  pattern(size_t r) const { return records.at(r).pattern; }
    Real        E(size_t r) const { return records.at(r).E; }
    Real       nu(size_t r) const { return records.at(r).nu; }
    Real        A(size_t r) const { return records.at(r).anisotropy; }
    size_t numParams() const {
        if (records.size() == 0) { return 0; }
        size_t result = records[0].numParams();
        for (const auto &r : records) {
            if (r.numParams() != result) throw std::runtime_error("All LUT entries must have the same number of paramters for this operation.");
        }
        return result;
    }

    // All of this is extremely inefficient for now!
    IsotropicLookupTable selectPattern(size_t pattern) {
        IsotropicLookupTable result;
        result.records.reserve(records.size());
        for (const auto &r : records) {
            if (r.pattern == pattern)
                result.records.push_back(r);
        }
        return result;
    }
    IsotropicLookupTable subset(const std::vector<size_t> idxs) {
        IsotropicLookupTable result;
        result.records.reserve(idxs.size());
        for (size_t i : idxs) result.records.push_back(records.at(i));
        return result;
    }

    size_t size() const { return records.size(); }

    // Allow type conversion of fields
    template<typename T> void  getEs(                   std::vector<T> &out) const { out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.E); }
    template<typename T> void getNus(                   std::vector<T> &out) const { out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.nu); }
    template<typename T> void  getAs(                   std::vector<T> &out) const { out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.anisotropy); }
    template<typename T> void  getParamValues(size_t p, std::vector<T> &out) const {
        if (p >= numParams()) {
            throw std::runtime_error("Invalid parameter " + std::to_string(p));
        }
                                                                                     out.clear(), out.reserve(size()); for (const auto &r : records) out.push_back(r.patternParams.at(p));
    }

    template<typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  getParamValues() const {
        size_t np = numParams();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> values(np, size());
        for (size_t i = 0; i < size(); ++i) {
            for (size_t p = 0; p < np; ++p) {
                values(p, i) = records[i].patternParams[p];
            }
        }
        return values;
    }

    template<typename T> std::vector<T>  getEs()                  const { std::vector<T> result;              getEs(result); return result; }
    template<typename T> std::vector<T> getNus()                  const { std::vector<T> result;             getNus(result); return result; }
    template<typename T> std::vector<T>  getAs()                  const { std::vector<T> result;              getAs(result); return result; }
    template<typename T> std::vector<T>  getParamValues(size_t p) const { std::vector<T> result;  getParamValues(p, result); return result; }

    std::vector<Record> records;
};


#endif /* end of include guard: LOOKUPTABLE_HH */
