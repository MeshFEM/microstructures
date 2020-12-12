////////////////////////////////////////////////////////////////////////////////
// pattern_info.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Gets useful information about a pattern topology
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/14/2016 15:32:28
////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/WireMesh.hh>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cassert>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2) {
        cerr << "usage: ./pattern_info pattern.obj" << endl;
        exit(-1);
    }

    string pat = argv[1];
    WireMesh<Symmetry::Orthotropic<>> wm(pat);

    vector<MeshIO::IOVertex> inVertices;
    vector<MeshIO::IOElement> elements;
    MeshIO::load(pat, inVertices, elements);

    vector<size_t> valences(inVertices.size(), 0);
    for (const auto &e : elements) {
        assert(e.size() == 2);
        ++valences.at(e[0]);
        ++valences.at(e[1]);
    }

    size_t maxValence = *max_element(valences.begin(), valences.end());
    cout << "num edges\t" << elements.size() << endl;
    cout << "max valence\t" << maxValence << endl;

    cout << "parameters\t" << wm.numParams() << endl;

    cout << "(under orthotropic symmetry constraints)" << endl;
    cout << "position parameters\t" << wm.numPositionParams() << endl;
    cout << "thickness parameters\t" << wm.numThicknessParams() << endl;
    cout << "blending parameters\t" << wm.numBlendingParameters() << endl;

    return 0;
}
