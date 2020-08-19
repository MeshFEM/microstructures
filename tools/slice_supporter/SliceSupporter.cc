////////////////////////////////////////////////////////////////////////////////
// SliceSupporter.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Creates a sequence of slices holding an "easy peel" foundation layer and
//      support structure for a given first slice.
//      Note: must be run from this directory since it accesses images stored
//      here.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/07/2017 22:37:28
////////////////////////////////////////////////////////////////////////////////
#include <CImg.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <deque>

#define SUPPORT_WIDTH_2X 0

namespace ci = cimg_library;

struct Loc {
    Loc(size_t c, size_t r)
        : col(c), row(r) { }
    Loc() : col(0), row(0) { }
    size_t col, row;
};

struct VisitedArray {
    template<typename T>
    VisitedArray(const ci::CImg<T> &img)
        : data(img.height() * img.width()),
          m_width(img.width()) { }

    void     visit(size_t col, size_t row)       { data[row * m_width + col] = true; }
    bool isVisited(size_t col, size_t row) const { return data[row * m_width + col]; }

    void visit(const Loc &l) { visit(l.col, l.row); }
    bool isVisited(const Loc &l) const { return isVisited(l.col, l.row); }

    void reset() { data.assign(data.size(), false); }

    size_t m_width;
    std::vector<bool> data;
};


int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cerr << "usage: SliceSupporter path_to_first_slice.bmp output_dir" << std::endl;
        exit(-1);
    }
    const std::string firstSlicePath(argv[1]);
    const std::string slice_dir(argv[2]);

    ci::CImg<unsigned char> first_slice(firstSlicePath.c_str());
    first_slice.channel(0); // Slice to grayscale.

    auto labels = first_slice.get_label(false);
    // labels.save("components.png");

    std::vector<ci::CImg<unsigned char>> supports;
    supports.reserve(100);
    for (size_t i = 1; i <= 29; ++i) {
        std::stringstream ss;
        ss << "support_"
           << std::setw(4) << std::setfill('0') << i
           << ".png";
        supports.emplace_back(ss.str().c_str());
    }

    if (SUPPORT_WIDTH_2X) {
        // for (auto &s : supports)
        //     s.resize((s.width() * 3) / 2, (s.height() * 3) /2, 1, 1, 1 /* nearest neighbor interp */);
        for (auto &s : supports)
            s.resize_doubleXY();
    }

    const size_t nComponents = size_t(labels.max()) + 1;
    std::cout << "Supporting " << nComponents - 1 << " components..." << std::endl;

    // Zeroth component should be black.
    assert((first_slice(0, 0) == 0) && (labels(0, 0) == 0));
    
    struct BBox {
        BBox()
            : min_r(std::numeric_limits<size_t>::max()),
              min_c(std::numeric_limits<size_t>::max()),
              max_r(0),
              max_c(0) { }

        size_t min_r, min_c, max_r, max_c;
        size_t rows() const { return max_r - min_r + 1; }
        size_t cols() const { return max_c - min_c + 1; }

        void update(size_t c, size_t r) {
            min_r = std::min(min_r, r);
            min_c = std::min(min_c, c);
            max_r = std::max(max_r, r);
            max_c = std::max(max_c, c);
        }

        void unionWith(const BBox &b) {
            min_r = std::min(min_r, b.min_r);
            min_c = std::min(min_c, b.min_c);
            max_r = std::max(max_r, b.max_r);
            max_c = std::max(max_c, b.max_c);
        }
    };

    std::vector<BBox> componentBBox(nComponents);
    for (size_t c = 0; c < labels.width(); ++c)
        for (size_t r = 0; r < labels.height(); ++r)
            componentBBox.at(labels(c, r)).update(c, r);
    assert(componentBBox.at(0).rows() == first_slice.height());
    assert(componentBBox.at(0).cols() == first_slice.width());

    ////////////////////////////////////////////////////////////////////////////
    // Optimize support placement
    ////////////////////////////////////////////////////////////////////////////
    const auto &supportSlice = supports.back();
    
    // Determine optimal placement for the support for each component by brute force:
    //  trying all positions, maximizing the the distance from the support
    //  pixels to boundary (unprinted pixels)
    // First, determine the (pre-offset) support pixel locations
    std::vector<Loc> supportPixelLocs;
    for (size_t c = 0; c < supportSlice.width(); ++c) {
        for (size_t r = 0; r < supportSlice.height(); ++r)
            if (supportSlice(c, r) != 0) supportPixelLocs.emplace_back(c, r);
    }

    VisitedArray va(first_slice);
    std::vector<Loc> supportLocations;
    for (size_t i = 1; i < componentBBox.size(); ++i) {
        const auto &cb = componentBBox[i];

        int bestMargin = std::numeric_limits<int>::min();
        Loc bestMarginLoc;
        // Loop over placements of the support in the full image (staying inside
        // the current component box.)
        // Note: the placement is of the support's center. This tries more
        // candidates than necessary, but is still fast enough.
        Loc supportCenter;
        int supportLeftOffset = -int(supportSlice.width() / 2);
        int supportTopOffset  = -int(supportSlice.height() / 2);
        for (supportCenter.col = cb.min_c; supportCenter.col <= cb.max_c; ++supportCenter.col) {
        for (supportCenter.row = cb.min_r; supportCenter.row <= cb.max_r; ++supportCenter.row) {
            // Col and row of upper left corner of the support image's candidate
            // placement in first_slice
            int slc = supportLeftOffset + int(supportCenter.col),
                str =  supportTopOffset + int(supportCenter.row);
            assert((slc >= 0) && (str >= 0));
            Loc supportTopLeft(slc, str);
            // Compute distance from support pixels to the component's boundary
            // using a BFS.
            int dist = 0;
            std::deque<Loc> bfsQueue;
            va.reset(); // we could improve the efficiency here.
            for (const Loc &l : supportPixelLocs) {
                Loc ol(supportTopLeft.col + l.col, supportTopLeft.row + l.row);
                if (first_slice(ol.col, ol.row) == 0) {
                    --dist;
                }
                assert(!va.isVisited(ol));
                va.visit(ol);
                bfsQueue.push_back(ol);
            }
            if (dist < 0) goto hitBoundary; // negative distance measure support pixels outside object; we want to minimize this count.
            dist = 0;
            while (!bfsQueue.empty()) {
                ++dist;
                std::deque<Loc> nextLocations;
                while (!bfsQueue.empty()) {
                    Loc l = bfsQueue.front();
                    bfsQueue.pop_front();
                    // Loop over neighbors
                    // Simplification: geometry should be well inside the image.
                    assert((l.col > 0) && (l.col < first_slice.width()));
                    assert((l.row > 0) && (l.row < first_slice.height()));
                    for (int nco = -1; nco <= 1; ++nco) {
                    for (int nro = -1; nro <= 1; ++nro) {
                        // TODO: what about other diagonal?!?
                        if ((nco == 0) == (nro == 0)) continue; // don't consider center or diagonals.
                        Loc nl(l.col + nco, l.row + nro);
                        if (first_slice(nl.col, nl.row) == 0) { // reached the boundary!
                            bfsQueue.clear();
                            nextLocations.clear();
                            goto hitBoundary;
                        }
                        if (va.isVisited(nl)) continue;
                        va.visit(nl);
                        nextLocations.push_back(nl);
                    }
                    }   
                }
                std::swap(bfsQueue, nextLocations);
            }
hitBoundary:
            // std::cerr << dist;
            if (dist > bestMargin) {
                bestMargin = dist;
                bestMarginLoc = Loc(supportTopLeft.col, supportTopLeft.row);
            }
        }
            // std::cerr << std::endl;
        }
        supportLocations.push_back(bestMarginLoc);
        // std::cout << "margin: " << bestMargin << std::endl;
    }

    // Compute bounding box of all support pixels.
    // Assumes that the first support structure column image is the largest
    // (which is always the case for tapering support).
    BBox supportPixelBB;
    const auto &baseSupport = supports.at(0);
    for (const Loc &sl : supportLocations) {
        for (size_t c = 0; c < baseSupport.width(); ++c) {
            for (size_t r = 0; r < baseSupport.height(); ++r) {
                if (baseSupport(c, r) != 0)
                    supportPixelBB.update(sl.col + c, sl.row + r);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Create foundation and support structure.
    ////////////////////////////////////////////////////////////////////////////
    // We use the "easy peel" foundation from B9. This is essentially a
    // rectangle with rounded corners. The very first slice is a rectangle with
    // its corners chopped off:
    //      +----------...
    //     /
    //    /
    //   +
    //   |
    //   |
    //   .
    // (The diagonal edge makes a 45 degree angle).
    // We place this diagonal edge to be tangent to the support column bounding
    // box. This wastes some material, but less than the B9 layout software.
    // The corner geometry is in the "corners_0001.png".."corners_0014.png"
    // images.
    const size_t nCornerLayers = 14;
    std::vector<ci::CImg<unsigned char>> corners;
    corners.reserve(nCornerLayers);
    for (size_t i = 1; i <= nCornerLayers; ++i) {
        std::stringstream ss;
        ss << "corner_" << std::setw(4) << std::setfill('0') << i << ".png";
        corners.emplace_back(ss.str().c_str());
    }

    size_t  cwidth = corners.back(). width(),
           cheight = corners.back().height();
    for (const auto &c : corners) {
        if ((c.width() != cwidth) || (c.height() != cheight))
            throw std::runtime_error("Mismatched corner image sizes");
    }

    // Offset the corners to make the diagonal edge tangent to the first
    // slice's bounding box. This offset in turn defines a rectangle that is the
    // bounding box of the foundation image. We form the image by filling this
    // box and drawing the corner images atop.
    size_t tangentRowOffset = size_t(ceil(cheight / 2.0)),
           tangentColOffset = size_t(ceil( cwidth / 2.0));
    BBox fullBBox;
    fullBBox.min_c = supportPixelBB.min_c - tangentColOffset;
    fullBBox.min_r = supportPixelBB.min_r - tangentRowOffset;
    fullBBox.max_c = supportPixelBB.max_c + tangentColOffset;
    fullBBox.max_r = supportPixelBB.max_r + tangentRowOffset;

    // Draw the support layers that include the corners.
    unsigned char color = 255;
    for (size_t i = 0; i < supports.size(); ++i) {
        auto supportImg = first_slice;
        supportImg.fill(0);
        if (i < 5) {
            supportImg.draw_rectangle(fullBBox.min_c, fullBBox.min_r,
                                      fullBBox.max_c, fullBBox.max_r,
                                      &color);
        }
        if (i < nCornerLayers) {
            supportImg.draw_image(fullBBox.min_c               , fullBBox.min_r                , corners.at(i)                 ); // Top left
            supportImg.draw_image(fullBBox.max_c - (cwidth - 1), fullBBox.min_r                , corners.at(i).get_mirror( "x")); // Top right
            supportImg.draw_image(fullBBox.min_c               , fullBBox.max_r - (cheight - 1), corners.at(i).get_mirror( "y")); // Bottom left
            supportImg.draw_image(fullBBox.max_c - (cwidth - 1), fullBBox.max_r - (cheight - 1), corners.at(i).get_mirror("xy")); // Bottom right
        }

        // Fill in the support structure columns
        const auto &s = supports.at(i);
        for (const Loc &sl : supportLocations) {
            for (size_t c = 0; c < s.width(); ++c) {
                for (size_t r = 0; r < s.height(); ++r) {
                    if (s(c, r) != 0) {
                        supportImg(sl.col + c, sl.row + r) = 255;
                    }
                }
            }
        }

        std::stringstream ss;
        ss << slice_dir << "/out_" << i << ".png";
        supportImg.save(ss.str().c_str());
    }
    
    return 0;
}
