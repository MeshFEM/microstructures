////////////////////////////////////////////////////////////////////////////////
// OffsetBoundsToTranslationBounds.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Output a pattern optimization job with per-translation variable bounds
//      set based on bounded offset around the pattern's default translation
//      variables.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/30/2016 03:10:24
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cstdlib>

#include <isosurface_inflator/IsosurfaceInflator.hh>
#include <pattern_optimization/PatternOptimizationJob.hh>

int main(int argc, const char *argv[]) {
    if (argc != 6) {
        std::cerr << "usage: ./OffsetBoundsToTranslationBounds mesher_name pattern job_file offset_min offset_max"
                  << std::endl;
        exit(-1);
    }

    std::string mesher(argv[1]);
    std::string patternPath(argv[2]);
    std::string jobPath(argv[3]);

    Real offsetMin = std::stod(argv[4]),
         offsetMax = std::stod(argv[5]);

    const bool vertexThickness = true;
    IsosurfaceInflator inflator(mesher, vertexThickness, patternPath);

    auto paramJob = PatternOptimization::parseJobFile(jobPath);
    if ((paramJob->varLowerBounds.size() + paramJob->varUpperBounds.size()) > 0) {
        std::cerr << "WARNING: overwriting per-variable bounds." << std::endl;
        paramJob->varLowerBounds.clear();
        paramJob->varUpperBounds.clear();
    }

    auto defaultParams = inflator.defaultParameters();
    for (size_t p = 0; p < defaultParams.size(); ++p) {
        if (inflator.isPositionParam(p)) {
            paramJob->varLowerBounds.emplace(p, defaultParams[p] + offsetMin);
            paramJob->varUpperBounds.emplace(p, defaultParams[p] + offsetMax);
        }
    }

    paramJob->writeJobFile(std::cout);

    return 0;
}
