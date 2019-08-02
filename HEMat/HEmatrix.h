//!
//! @file       HEmatrix.h, header file
//! @brief      defining functions for matrix operations
//!
//! @author     Miran Kim
//! @date       Dec. 1, 2017
//! @copyright  GNU Pub License
//!

#ifndef HEMatrix_H_
#define HEMatrix_H_


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

#include "../src/Scheme.h"
#include "../src/SecretKey.h"

using namespace std;
using namespace NTL;

//! structure for parameters
typedef struct HEMatpar{
    long nrows;      //! real dimension of data
    long ncols;
    
    long dim;        //! power-of-two integer larger than dim
    long dim1;
    long subdim;     //! it can be used for denoting the number of rows (or columns) in rectangular matrix multiplication
    
    long logdim; 
    long nslots;     //! the number of slots (Power of Two)
    long sqrdim;     //! sqrt(d)
    long nbatching;  //! number of matrices in a single ciphertext (default = 1)
    
    long pBits;      //! bits for scaling factor of message
    long cBits;      //! bits for scaling factor of a constant
    long logQ;       //! bitlength for fresh ciphertext
}HEMatpar;

//!@ Read parameters for matrix over HE
void readHEMatpar(HEMatpar& HEmatpar, long nrows, long ncols, long pBits, long cBits, long logQ, const long subdim = 0, const long nbatching = 1);


class HEmatrix {
    
public:
    Scheme& scheme;
    SecretKey& secretKey;
    HEMatpar& HEmatpar;
    
    //! constructor
    HEmatrix(Scheme& scheme, SecretKey& secretKey, HEMatpar& HEmatpar) : scheme(scheme), secretKey(secretKey), HEmatpar(HEmatpar) {}
    
    //! encrypt and decrypt
    void encryptRmat(Ciphertext& ctxt, mat_RR& mat, long logp);
    void decryptRmat(mat_RR& mat, Ciphertext& ctxt);
    
    void encryptParallelRmat(Ciphertext& ctxt, mat_RR*& mat, long logp, long nbatching);
    void decryptParallelRmat(mat_RR*& mat, Ciphertext& ctxt);

    //! optimization: used for linear transformation
    void msgleftRotate(complex<double>*& res, complex<double>* vals, long dim, long nrot);
    void msgrightRotate(complex<double>*& res, complex<double>* vals, long dim, long nrot);
    
    void msgleftRotateAndEqual(complex<double>*& vals, long dim, long nrot);
    void msgrightRotateAndEqual(complex<double>*& vals, long dim, long nrot);


    //! transpose
    void genTransPoly(ZZX*& transpoly);
    void transpose(Ciphertext& res, Ciphertext& ctxt, ZZX*& transpoly);
    
    void genTransPoly_Parallel(ZZX*& transpoly);
    void transpose_Parallel(Ciphertext& res, Ciphertext& ctxt, ZZX*& transpoly);
    
    //! shift by k-position columns
    void genShiftPoly(ZZX*& shiftpoly, const long num = 0);
    void shiftBycols(Ciphertext& res, Ciphertext& ctxt, long k, ZZX*& shiftpoly);
   
    //! parallel-multiplication
    void genShiftPoly_Parallel(ZZX*& shiftpoly);
    void shiftBycols_Parallel(Ciphertext& res, Ciphertext& ctxt, long k, ZZX*& shiftpoly);
    
    
    //-------------------------------------------
    // multiplication
    //-------------------------------------------
    
    //! generate polynomials for linear transformations of multiplication
    void genMultPoly(ZZX**& Initpoly);
    void genMultPoly_Parallel(ZZX**& Initpoly); //! parallel
    void genMultBPoly(ZZX*& Initpoly);  //! only for "B"
    
    //! generate the initial ciphertexts for multiplicaiton
    void genInitCtxt(Ciphertext& resA, Ciphertext& resB, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly);
    void genInitCtxt_Parallel(Ciphertext& resA, Ciphertext& resB, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly);
    
   
    void genInitActxt(Ciphertext*& Actxts, Mat<RR>& mat);    //! generate the encryptions of Amat (d ctxts)
    void genInitBctxt(Ciphertext& resB, Ciphertext& Bctxt, ZZX*& Initpoly); //!  generate the initial ciphertexts (only for "B")
    void genInitRecActxt(Ciphertext*& Actxts, Mat<RR>& mat);  //! multiplication when encryptions of rectangular Amat are given as fresh
    
    //! perform Hadamard multiplication
    void HEmatmul_Hadamard(Ciphertext& res, Ciphertext* Actxts, Ciphertext* Bctxts, long num);
    
    //! performt the matrix multiplication
    void HEmatmul(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly);
    void HEmatmul_Parallel(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly);
    void HErmatmul(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly);    //! rectangular matrix multiplication
    
    
    //! performt the matrix multiplication where "A" are given as fresh ciphertexts with a form of mat-mult algorithm
    void HEmatmul_preprocessing(Ciphertext& res, Ciphertext*& Actxts, Ciphertext& Bctxt, ZZX*& Initpoly);
    void HErmatmul_preprocessing(Ciphertext& res, Ciphertext*& Actxts, Ciphertext& Bctxt, ZZX*& Initpoly);
    
};


#endif
