////////////////////////////////////////////////////////////////////////////////
// AdaptiveEvaluator.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Temporary hack to speed up signed distance evaluation: don't sample
//      regions away from the contour. Eventually this should be
//      replaced/supplemented with an truly adaptive mesher.
//
//      The approach ensures that, assuming all contour components are "seen"
//      by the coarsened grid, an identical fine contour is extracted to the
//      one given by evaluating every grid point.
//
//      The approach is to first evaluate on a coarse grid, then recursively
//      refine the cells that detect a contour. Finally, a BFS is run on the
//      fine grid along the contour, starting from each fine contour cell
//      found by recursive refinement.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  04/04/2016 23:09:34
////////////////////////////////////////////////////////////////////////////////
#ifndef ADAPTIVEEVALUATOR_HH
#define ADAPTIVEEVALUATOR_HH

#include <vector>
#include <array>
#include <bitset>
#include <queue>
#include "Grid.hh"

struct CellChunk {
        size_t offsetRow, offsetCol;
        size_t chunkRowSize, chunkColSize;

        CellChunk(size_t ofr, size_t ofc, size_t crs, size_t ccs)
            : offsetRow(ofr), offsetCol(ofc),
              chunkRowSize(crs), chunkColSize(ccs) {
            assert((chunkRowSize > 0) && (chunkColSize > 0));
        }

        void getCornerVertices(const Grid2D &fullGrid,
                               std::array<size_t, 4> &corners) const {
            // A vertex's 2D index coincides with the 2D index of the cell
            // having the vertex as a lower-left corner.
            corners[0] = fullGrid.get1DVertexIndex(offsetRow,                offsetCol);
            corners[1] = fullGrid.get1DVertexIndex(offsetRow,                offsetCol + chunkColSize);
            corners[2] = fullGrid.get1DVertexIndex(offsetRow + chunkRowSize, offsetCol + chunkColSize);
            corners[3] = fullGrid.get1DVertexIndex(offsetRow + chunkRowSize, offsetCol);
        }

        bool isFineCell() const { return (chunkRowSize == 1) && (chunkColSize == 1); }

        // Up to four sub-chunks
        void subchunks(std::vector<CellChunk> &chunks) const {
            chunks.clear();
            size_t numSubRows = (chunkRowSize > 1) ? 2 : 1,
                   numSubCols = (chunkColSize > 1) ? 2 : 1;

            // Nothing to subdivide
            if ((numSubRows == 1) && (numSubCols == 1)) return;

            // Size of a full sub-chunk (round up)
            // Warning: the last sub-chunk may be partial for non pow-of-2 sizes.
            size_t subRowSize = (chunkRowSize + 1) / 2,
                   subColSize = (chunkColSize + 1) / 2;

            chunks.reserve(numSubRows * numSubCols);
            for (size_t srOffset = 0; srOffset < chunkRowSize; srOffset += subRowSize) {
                // Compute the actual (possibly partial) subchunk row size
                size_t srSize = std::min(subRowSize, chunkRowSize - srOffset);
                for (size_t scOffset = 0; scOffset < chunkColSize; scOffset += subColSize) {
                    // Compute the actual (possibly partial) subchunk col size
                    size_t scSize = std::min(subColSize, chunkColSize - scOffset);
                    chunks.emplace_back(offsetRow + srOffset, offsetCol + scOffset,
                                        srSize, scSize);
                }
            }
        }
};

struct AdaptiveEvaluator {
    template<class Domain>
    AdaptiveEvaluator(const Domain &d, const Grid2D &fullGrid, size_t coarseningLevels) {
        // Initially, no cells have been visited/determined to overlap contour
        m_tmp_cellVisited.assign(fullGrid.numCells(), false);
        m_cellOverlapsContour.assign(fullGrid.numCells(), false);

        // No vertices are evaluated yet.
        m_vtxSD.resize(fullGrid.numVertices());
        m_vtxSDEvaluated.assign(fullGrid.numVertices(), false);

        // Evaluate and adaptively refine all coarse chunks
        size_t chunkSize = 1 << coarseningLevels;
        const size_t nrows = fullGrid.rows(), ncols = fullGrid.cols();
        for (size_t offsetRow = 0; offsetRow < nrows; offsetRow += chunkSize) {
            size_t rsize = std::min(chunkSize, nrows - offsetRow);
            for (size_t offsetCol = 0; offsetCol < ncols; offsetCol += chunkSize) {
                size_t csize = std::min(chunkSize, ncols - offsetCol);
                m_recursiveEvalAndRefine(d, fullGrid, CellChunk(offsetRow,
                    offsetCol, rsize, csize));
            }
        }

        // Run BFS of neighbors from all fine cells overlapping the contour
        // to find connected contour cells missed by the coarse grid.
        std::queue<size_t> bfsQueue;
        for (size_t ci = 0; ci < fullGrid.numCells(); ++ci)
            if (m_cellOverlapsContour[ci]) bfsQueue.push(ci);

        // Note: cells are marked visited *before* they are added to the bfs queue.
        Grid2D::AdjacencyVec cornerVertices;
        while (!bfsQueue.empty()) {
            const size_t ci = bfsQueue.front();
            bfsQueue.pop();
            assert(m_tmp_cellVisited.at(ci));

            // Sample the cell corners to determine if the cell is actually
            // inside the object (even though we may already know this)
            fullGrid.cellVertices(ci, cornerVertices);
            std::bitset<4> inside;
            for (size_t i = 0; i < 4; ++i)
                inside[i] = (m_sampleVtx(d, fullGrid, cornerVertices[i]) <= 0);
            bool isContour = !(inside.all() || inside.none());
            m_cellOverlapsContour.at(ci) = isContour;
            if (!isContour) continue;

            // For each edge cut by the contour, visit the corresponding
            // neighbor if it hasn't been visited already.
            // Follow Grid2D's vertex numbering:
            // 3--2--2   increasing row idx
            // |     |   ^
            // 3     1   |
            // |     |   |
            // 0--0--1   +-----> increasing col idx
            size_t row, col;
            fullGrid.get2DCellIndex(ci, row, col);
            for (size_t ei = 0; ei < 4; ++ei) {
                if (inside[ei] == inside[(ei + 1) % 4]) continue;
                int neighborRow = int(row), neighborCol = int(col);
                if      (ei == 0)  --neighborRow;
                else if (ei == 1)  ++neighborCol;
                else if (ei == 2)  ++neighborRow;
                else if (ei == 3)  --neighborCol;
                // Ignore traverals outside the grid.
                if ((neighborRow < 0) || (size_t(neighborRow) >= nrows)) continue;
                if ((neighborCol < 0) || (size_t(neighborCol) >= ncols)) continue;
                size_t ni = fullGrid.get1DCellIndex(neighborRow, neighborCol);

                // Mark unvisited neighbors and add them to the queue.
                if (m_tmp_cellVisited.at(ni)) continue;
                m_tmp_cellVisited[ni] = true;
                bfsQueue.push(ni);
            }
        }

        // We no longer need the cell visited flags.
        m_tmp_cellVisited.clear();
    }

    bool cellOverlapsContour(size_t cellIndex) const {
        return m_cellOverlapsContour.at(cellIndex);
    }

    bool vertexWasEvaluated(size_t vtxIndex) const { return m_vtxSDEvaluated.at(vtxIndex); }

    // Only call this on vertices whose signed distances were evaluated during
    // construction (Otherwise throws an assertion or returns zero if DNDEBUG)
    Real signedDistance(size_t vtxIndex) const {
        assert(vertexWasEvaluated(vtxIndex));
        return m_vtxSD.at(vtxIndex);
    }

private:
    template<class Domain>
    void m_recursiveEvalAndRefine(const Domain &d, const Grid2D &g, const CellChunk &cc) {
        // Sample corners
        std::array<size_t, 4> corners;
        cc.getCornerVertices(g, corners);
        std::bitset<4> inside;
        for (size_t i = 0; i < 4; ++i)
            inside[i] = (m_sampleVtx(d, g, corners[i]) <= 0);

        // Does the chunk contain the contour? (At least, as visible by
        // the coarse chunked grid).
        bool isContour = !(inside.all() || inside.none());

        // If we've refined down to a fine cell, mark it as visited and
        // mark its contour membership.
        if (cc.isFineCell()) {
            size_t ci = g.get1DCellIndex(cc.offsetRow, cc.offsetCol);
            m_tmp_cellVisited.at(ci) = true;
            m_cellOverlapsContour.at(ci) = isContour;
            return; // This refinement tree branch is done.
        }

        // Don't refine cells that are entirely inside or outside
        if (!isContour) return;

        std::vector<CellChunk> subchunks;
        cc.subchunks(subchunks);
        assert(subchunks.size() > 0); // Already should have bailed...
        for (const auto &sc : subchunks)
            m_recursiveEvalAndRefine(d, g, sc);
    }


    // Memoized vertex sdf sampling
    template<class Domain>
    Real m_sampleVtx(const Domain &d, const Grid2D &g, size_t vtxIndex) {
        if (m_vtxSDEvaluated.at(vtxIndex)) return m_vtxSD.at(vtxIndex);
        Real sd = d.signedDistance(g.vertexPosition(vtxIndex));
        m_vtxSD.at(vtxIndex) = sd;
        m_vtxSDEvaluated.at(vtxIndex) = true;
        return sd;
    }

    std::vector<Real> m_vtxSD;
    std::vector<bool> m_vtxSDEvaluated, m_cellOverlapsContour;
    std::vector<bool> m_tmp_cellVisited; // Temporary--freed after construction
};

#endif /* end of include guard: ADAPTIVEEVALUATOR_HH */
