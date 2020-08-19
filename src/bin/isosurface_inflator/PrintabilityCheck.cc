////////////////////////////////////////////////////////////////////////////////
// PrintabilityCheck.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Iterates through a table of patterns and divides them into two tables:
//      those printable and those unprintable.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/05/2017 00:29:14
////////////////////////////////////////////////////////////////////////////////

#include <isosurface_inflator/WireMesh.hh>
#include <pattern_optimization/LookupTable.hh>

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>

#include <sstream>
#include <iomanip>

int main(int argc, const char *argv[]) {
    if (argc != 4) {
        std::cerr << "usage: ./PrintabilityCheck table.txt printable_out.txt unprintable_out.txt" << std::endl;
        exit(-1);
    }

    const std::string patternTablePath = argv[1];
    const std::string    printablePath = argv[2];
    const std::string  unprintablePath = argv[3];

    using PSym = Symmetry::Orthotropic<>;
    using WMesh = WireMesh<PSym>;

    IsotropicLookupTable<double> table(patternTablePath);

    // Create wire mesh and printabilty graph for each pattern
    std::vector<std::unique_ptr<WMesh>> wireMeshes;

    for (size_t i = 0; i < table.size(); ++i) {
        size_t pat = table.pattern(i);
        if (pat >= wireMeshes.size())
            wireMeshes.resize(pat + 1);
        std::stringstream ss;
        ss << ".." // TODO: get $MICRO_DIR environment variable.
           << "/patterns/3D/reference_wires/pattern"
           << std::setfill('0') << std::setw(4) << pat
           << ".wire";
        if (!wireMeshes[pat])
            wireMeshes[pat] = Future::make_unique<WMesh>(ss.str());
    }

    IsotropicLookupTable<double> printableTable, unprintableTable;
    size_t numPrintable = 0, numUnprintable = 0;
    std::vector<bool> isPrintable(table.size(), false);
    for (size_t i = 0; i < table.size(); ++i) {
        size_t pat = table.pattern(i);
        const std::vector<double> &params = table.records[i].patternParams;
        isPrintable[i] = wireMeshes.at(pat)->isPrintable(params);
        if (isPrintable[i]) {
            printableTable.records.push_back(table.records[i]);
            ++numPrintable;
        }
        else {
            unprintableTable.records.push_back(table.records[i]);
            ++numUnprintable;
        }
    }

    std::cout << numPrintable << " printable, " << numUnprintable
              << " unprintable" << std::endl;

      printableTable.write(  printablePath);
    unprintableTable.write(unprintablePath);

    IsotropicLookupTable<double> printableDisagree, unprintableDisagree;

    // Validate constraints
    for (size_t i = 0; i < table.size(); ++i) {
        size_t pat = table.pattern(i);
        const std::vector<double> &params = table.records[i].patternParams;
        auto C = wireMeshes.at(pat)->selfSupportingConstraints(params);

        // Use position map to determine current z coords; we need to construct
        // homogeneous param vector.
        Eigen::VectorXd paramVec(params.size() + 1);
        for (size_t p = 0; p < params.size(); ++p) paramVec[p] = params[p];
        paramVec[params.size()] = 1.0;

        Eigen::VectorXd cval = C * paramVec;
        bool satisfied = true;
        for (size_t crow = 0; crow < size_t(C.rows()); ++crow) {
            if (cval[crow] < -1e-12) satisfied = false;
        }
        if (satisfied != isPrintable[i]) {
            if (isPrintable[i])   printableDisagree.records.push_back(table.records[i]);
            else                unprintableDisagree.records.push_back(table.records[i]);
            std::cout << "disagreement with cval = " << cval.transpose() << std::endl;
        }
    }

      printableDisagree.write(  printablePath + ".disagree");
    unprintableDisagree.write(unprintablePath + ".disagree");

    return 0;
}
