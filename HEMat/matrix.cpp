//!
//! @file       matrix.cpp
//! @brief      implementing functions for matrix operations in plaintext
//!
//! @author     Miran Kim
//! @date       Dec. 1, 2017
//!

#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>

#include <cmath>
#include <map>
#include <math.h>  // pow
#include <sys/time.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"
#include <cassert>
#include <random>
#include <string>
#include <iomanip>

#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include "NTL/RR.h"
#include <NTL/ZZX.h>
#include "NTL/mat_RR.h"
#include "NTL/vec_RR.h"

#include "matrix.h"

//!@ Input: vec_RR
//!@ Function: print the vector
//!@ If print_size = 0, then print out all the components of an input vector

void printRvector(vec_RR& vec, long print_size){
    long len;
    
    if(print_size == 0){
        len = vec.length();
    }
    else{
        len = print_size;
    }
    
    cout << "   [" ;
    for(int i = 0; i < len; ++i){
        cout << " " << vec[i] << ((i != len - 1) ? "\t" : "]\n");
    }
}


//!@ Input: RR-matrix
//!@ Function: print the matrix
void printRmatrix(Mat<RR>& mat, const long print_size){
    long rlen, clen;
    
    if(print_size == 0){
        rlen = mat.NumRows();
        clen = mat.NumCols();
    }
    else{
        rlen = print_size;
        clen = print_size;
    }
    
    for(int i = 0; i< rlen; ++i){
        cout << "   [";
        for(int j = 0; j < clen; ++j){
            cout << mat[i][j] << ((j != clen - 1) ? "\t" : "]\n");
        }
    }
}

//!@ Input: A and B
//!@ Function: return the maximum norm of the difference of two input matrices A and B
RR getError(mat_RR Amat, mat_RR Bmat, long nrows, long ncols){
    RR ret = to_RR("0");
    
    for(long i = 0; i < nrows; ++i){
        for(long j = 0; j < ncols; ++j){
            RR temp = abs(Amat[i][j]-Bmat[i][j]);
            if (ret < temp){
                ret = temp;
            }
            if(temp > 1e-2){
                cout << "(" << i << "," << j  << ") = " << Amat[i][j] << ", " << Bmat[i][j]<< endl;
            }
        }
    }
    return ret;
}

