//!
//! @file       matrix.h, header file
//! @brief      defining functions for matrix operations
//!
//! @author     Miran Kim
//! @date       Dec. 1, 2017
//! @copyright  GNU Pub License
//!

#ifndef Matrix_H_
#define Matrix_H_

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>

#include <NTL/BasicThreadPool.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include "NTL/RR.h"
#include <NTL/ZZX.h>
#include "NTL/mat_RR.h"
#include "NTL/vec_RR.h"

using namespace std;
using namespace NTL;
  
void printRvector(vec_RR& vec, const long print_size = 0);

void printRmatrix(Mat<RR>& mat, const long print_size = 0);

RR getError(mat_RR Amat, mat_RR Bmat, const long nrows = 0, const long ncols = 0);


#endif
