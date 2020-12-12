#ifndef RREF_H
#define RREF_H

#include <float.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>

////////////////////////////////////////////////////////////////////////////////
/*! Compute reduced row echelon form.
//  Assumes by default we have infinite precision (no tolerance value).
//  Variables of lowest indices possible are chosen as dependent variables.
//  @param[inout] C          Matrix provided as a vector of rows. Overwritten
//                           with RREF. All rows must have equal lengths.
//  @param[out]   dep_idx    Dependent variable column indices (which column of
//                           RREF matrix corresponds to the variable)
//                           (listed in sorted order)
//  @param[out]   dep_row    Dependent variable row indices (which row of the
//                           RREF matrix holds the variable's dependencies)
//  @param[in]    forceIndep which variables are forced to be independent
//                           (defaults to none)
//  @param[in]    tol        Threshold at which a value is considered zero
//                           (defaults to 0)
*///////////////////////////////////////////////////////////////////////////////
template <typename CType> 
void rref(std::vector<std::vector<CType> > &C,
          std::vector<size_t> &dep_idx,
          std::vector<size_t> &dep_row,
          const std::vector<bool> &forceIndep = std::vector<bool>(),
          Real tol = 0.0)
{
    int m = C.size();
    if (m < 1) return;
    int n = C[0].size(); 
    if (n < 1) return;
    for (int l = 1; l < m; l++)
        assert(C[l].size() == (size_t) n);
    int i = 0, j = 0;
    dep_idx.clear();

    bool forcingIndependentVars = (forceIndep.size() == size_t(n));

    while ((i < m) && (j < n)) {
        // Skip over variables we force to be independent
        if (forcingIndependentVars && forceIndep[j]) {
            ++j;
            continue;
        }

        // Find value and index of the smallest nonzero element in the remainder
        // of column j.
        // We choose the smallest nonzero value (above the tolerance) to reduce
        // fractional entries.
        int k = -1; 
        Real minNonzeroVal = std::numeric_limits<Real>::max();
        for (int l = i; l < m; ++l) {
            Real tmp = std::abs(C[l][j]);
            if ((tmp > tol) && (tmp < minNonzeroVal)) {
                minNonzeroVal = tmp;
                k = l;
            }
        }
        if (k == -1) { 
            // all zero column, reset to exactly zero
            for (int l = i; l < m; ++l)
                C[l][j] = 0; 
            // Keep looking for a nonzero pivot element
            j++;
            continue;
        }
        else {
            // record pivot index (= dependent variable)
            dep_idx.push_back(j);
            // Record the row expressing the dependent variable's dependencies
            dep_row.push_back(i);

            for (int l = j; l < n; ++l)
                std::swap(C[i][l], C[k][l]);

            // Rescale pivot row so pivot element is 1
            CType scale = (CType) 1.0 / C[i][j];  
            C[i][j] = 1;
            for (int l = j + 1; l < n; ++l)
                C[i][l] *= scale;   

            // Remove dependent variable from all other constraints
            for (int s = 0; s < m; ++s) {
                if (s != i) {
                    CType row_scale = C[s][j];
                    for (int l = j; l < n; ++l)
                        C[s][l] -= row_scale * C[i][l];
                }
            }
            i++; j++; 
        }
    }

    if (forcingIndependentVars) {
        // Make sure there are no independent variables in rows without
        // dependent variables (i.e. the rows i .. n - 1)
        for (int row = i; row < n; ++row) {
            for (int l = 0; l < n; ++l) {
                if (std::abs(C[row][l]) > tol)  {
                    // The forced independent variable set was invalid!
                    assert(false);
                }
            }
        }
    }
}

#endif // RREF_H
