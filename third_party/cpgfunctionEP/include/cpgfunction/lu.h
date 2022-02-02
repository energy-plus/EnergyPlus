//
// Created by jackcook on 8/17/21.
//

#include <cmath>
#include <stdexcept>
#include <vector>

#ifndef LU_DECOMPOSITION_LU_H
#define LU_DECOMPOSITION_LU_H

namespace jcc {

void CroutDecomposition(std::vector<double> &A, int &n, std::vector<int> &indx,
                        int &n_threads);

void CroutSolve(std::vector<double> &LU, std::vector<double> &b, int &n,
                std::vector<int> &indx, int &n_threads);

double dot(int &n, std::vector<double> &A, int &begin_x, int &incx,
           int &begin_y, int &incy, int &n_threads);

} // namespace jcc

#endif // LU_DECOMPOSITION_LU_H
