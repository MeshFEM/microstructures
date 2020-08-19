////////////////////////////////////////////////////////////////////////////////
// RemoveUnsupported.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Detect and remove unsupported components from sliced images.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/10/2017 00:39:59
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cassert>
#include <vector>

#include <CImg.h>
#include <Future.hh>

int main(int argc, const char *argv[]) {
    if (argc < 2) {
        std::cerr << "usage: RemoveUsupported image1 image2..." << std::endl;
        exit(-1);
    }

    using Img = cimg_library::CImg<unsigned char>;
    std::unique_ptr<Img> prev;
    for (size_t i = 1; i < argc; ++i) {
        std::string path(argv[i]);
        std::unique_ptr<Img> curr = Future::make_unique<Img>(path.c_str());

        curr->channel(0); // make bitmap
        if (prev) {
            const size_t width  = curr->width(),
                         height = curr->height();

            assert(width  == prev->width());
            assert(height == prev->height());
            auto origSlice = *curr;

            curr->label(true /* consider diagonal voxels supported */); // Overwrite image with connected component labels.
            const size_t nComponents = size_t(curr->max()) + 1;
            if (nComponents > 250) throw std::runtime_error("Too many components"); // make sure we're not near overflow
            
            // Determine which components are supported.
            // A component is supported if any of its pixels are.
            std::vector<bool> hasSupport(nComponents, false);
            for (size_t c = 0; c < width; ++c) {
                for (size_t r = 0; r < height; ++r) {
                    unsigned char label = (*curr)(c, r);
                    if (label != 0)
                        if ((*prev)(c, r) != 0) hasSupport.at(label) = true;
                }
            }

            // Clear out the unsupported components.
            size_t numRemoved = 0;
            for (size_t c = 0; c < width; ++c) {
                for (size_t r = 0; r < height; ++r) {
                    if (origSlice(c, r) == 0) {
                        // Prevent enclosed voids (counted as separate components)
                        // from filling in.
                        (*curr)(c, r) = 0;
                        continue;
                    }
                    unsigned char label = (*curr)(c, r);
                    unsigned char newLabel = hasSupport.at(label) ? 255 : 0;
                    (*curr)(c, r) = newLabel;
                    if ((label > 0) != (newLabel > 0)) ++numRemoved;
                }
            }
            if (numRemoved)
                std::cerr << "WARNING: removed " << numRemoved << " voxels from " << path << std::endl;

        }
        std::string outPath = path.substr(0, path.size() - 3) + "cleaned.png";
        // std::cerr << "writing " << outPath << std::endl;;
        curr->save(outPath.c_str());
        prev = std::move(curr);
    }
}
