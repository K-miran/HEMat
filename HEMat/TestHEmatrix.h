//!
//! @file       TestHEmatrix.h, header file
//! @brief      defining functions for testing homomorphic matrix computation
//!
//! @author     Miran Kim
//! @date       Dec. 1, 2017
//! @copyright  GNU Pub License
//!


#ifndef HEAANNTT_TESTHEMATRIX_H_
#define HEAANNTT_TESTHEMATRIX_H_

#include <iostream>

using namespace std;

class TestHEmatrix {
public:
	
    static void testHEAAN(long logN, long logQ);
    
    static void testEnc(long dim);
    
    static void testAdd(long dim);

    static void testTrans(long dim);
    
    static void testShift(long dim, long k);
    
    static void testMult(long dim);
    
    static void testRMult(long nrows, long subdim);
    
    static void testMult_preprocessing(long dim);
    
    static void testRMult_preprocessing(long nrows, long subdim);
    
    static void testSIMDAdd(long dim, long nbatching, const long nitr);
    
    static void testSIMDTrans(long dim, long nbatching, const long nitr);
    
    static void testSIMDMult(long dim, long nbatching, const long nitr);
    
};

#endif /* TESTHEMATRIX_H_ */
