////////////////////////////////////////////////////////////////////////////////
// sort.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Sort the enumerated patterns by complexity. First by edge count, then by
//      max valence.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/14/2016 03:38:41
////////////////////////////////////////////////////////////////////////////////
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/Geometry.hh>
#include <MeshFEM/utils.hh>

#include <stdexcept>
#include <iostream>
#include <fstream>

using namespace std;

// num edges, valence
using PatternRecord = std::pair<size_t, size_t>;

int main(int argc, char *argv[])
{
    if (argc != 2) {
        cerr << "usage: ./match directory.txt" << endl;
        exit(-1);
    }

    string dirpath = argv[1];
    ifstream dir(dirpath);
    if (!dir) throw std::runtime_error("Couldn't open directory " + dirpath);

    vector<PatternRecord> patternRecords;
    vector<string>        paths;

    string path;
    while (getline(dir, path)) {
        paths.push_back(path);
        vector<MeshIO::IOVertex> inVertices;
        vector<MeshIO::IOElement> elements;
        MeshIO::load(path, inVertices, elements);

        size_t numEdges = elements.size();
        vector<size_t> valences(inVertices.size(), 0);
        for (const auto &e : elements) {
            assert(e.size() == 2);
            ++valences.at(e[0]);
            ++valences.at(e[1]);
        }

        size_t maxValence = *max_element(valences.begin(), valences.end());

        patternRecords.push_back({numEdges, maxValence});
    }

    // Sort patterns lexicographically
    vector<size_t> perm = sortPermutation(patternRecords);
    // for (size_t p : perm) { cout << p << ":\t" << patternRecords.at(p).first << ", " << patternRecords.at(p).second << endl; }
    for (size_t p : perm) { cout << paths.at(p) << endl; }

    return 0;
}
