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
    
    /**
     @param[in] mat, The input matrix
     @param[in] logp, The scale of messages
     @param[out] ctxt, The ciphertext encrypting a scaled matrix "(1 << logp) * mat"
     */
    void encryptRmat(Ciphertext& ctxt, mat_RR& mat, long logp);
    
    /**
     @param[in] ctxt, The input ciphertext encrypting "(1 << logp) * mat"
     @param[out] mat, The matrix
     */
    void decryptRmat(mat_RR& mat, Ciphertext& ctxt);
    
    /**
     @param[in] mat, The multiple input matrices
     @param[in] logp, The scale of messages
     @param[in] nbathcing, The number of multiple matrices in a single ciphertext
     @param[out] ctxt, The ciphertext encrypting multiple matrices with scale (1 << logp)
     */
    void encryptParallelRmat(Ciphertext& ctxt, mat_RR*& mat, long logp, long nbatching);
    
    /**
     @param[in] ctxt, The input ciphertext encrypting "(1 << logp) * mat"
     @param[out] mat, Multiple matrices
     */
    void decryptParallelRmat(mat_RR*& mat, Ciphertext& ctxt);

    /**
     @param[in] vals, The input vector, vals = (vals[0],vals[1],...,vals[d-1])
     @param[in] dim, The dim of of the input vector
     @param[in] nrot, The number of steps to rotate by left (positive)
     @param[out] res, The output vector s.t.
     res = (vals[nrot], vals[nrot+1],..., vals[d-1]), vals[0], vals[1],...., vals[nrot-1])
     */
    void msgleftRotate(complex<double>*& res, complex<double>* vals, long dim, long nrot);
    
    /**
     @param[in] vals, The input vector, vals = (vals[0],vals[1],...,vals[d-1])
     @param[in] dim, The dim of of the input vector
     @param[in] nrot, The number of steps to rotate by right (positive)
     @param[out] res, The output vector s.t
     res = (vals[d-nrot], vals[d-nrot+1],..., vals[d-1]), vals[0], vals[1],...., vals[d-nrot-1])
     */
    void msgrightRotate(complex<double>*& res, complex<double>* vals, long dim, long nrot);
    
    /**
     @param[in] vals, The input vector to rotate
     @param[in] dim, The dim of of the input vector
     @param[in] nrot, The number of steps to rotate by left (positive)
     */
    void msgleftRotateAndEqual(complex<double>*& vals, long dim, long nrot);
    
    /**
     @param[in] vals, The input vector to rotate
     @param[in] dim, The dim of of the input vector
     @param[in] nrot, The number of steps to rotate by right (positive)
     */
    void msgrightRotateAndEqual(complex<double>*& vals, long dim, long nrot);

    /**
     @param[out] transpoly, The polynomials needed for transposition
     */
    void genTransPoly(ZZX*& transpoly);
    
    /**
     @param[in] ctxt, The input ciphertext
     @param[in] transpoly, The polynomials
     @param[out] res, The output ciphertext that encrypts the transpose of the corresponding input matrix
     consumes a constant-multiplication level
     */
    void transpose(Ciphertext& res, Ciphertext& ctxt, ZZX*& transpoly);
    
    /**
     @param[out] transpoly, The polynomials needed for transpositions of multiple matrices
     */
    void genTransPoly_Parallel(ZZX*& transpoly);
    
    /**
     @param[in] ctxt, The input ciphertext encrypting multiple matrices
     @param[in] transpoly, The polynomials
     @param[out] res, The output ciphertext that encrypts the transpose results of the multiple input matrices
     */
    void transpose_Parallel(Ciphertext& res, Ciphertext& ctxt, ZZX*& transpoly);
    
    /**
     @param[in] num, The number of steps to shift by columns
     @param[out] shiftpoly, The polynomials needed for column shift by num position
     num = 0: shift by (d-1)
     */
    void genShiftPoly(ZZX*& shiftpoly, const long num = 0);
    
    /**
     @param[in] ctxt, The input ciphertext, Enc(m[0],...m[d-1] | m[d],...m[2d-1] | ... )
     @param[in] k, The number of steps to shift by columns
     @param[in] shiftpoly, The polynomials
     @param[out] res, The output ciphertext, Enc(m[k] ... m[k-1] | m[d+k]...m[d+k-1] | ... )
     consumes a constant-multiplication level
     */
    void shiftBycols(Ciphertext& res, Ciphertext& ctxt, long k, ZZX*& shiftpoly);
   
    /**
     @param[out] shiftpoly, The polynomials needed for column shift of multiple matrices
     */
    void genShiftPoly_Parallel(ZZX*& shiftpoly);
    
    /**
     @param[in] ctxt, The input ciphertext
     @param[in] k, The number of steps to shift by columns
     @param[in] shiftpoly, The polynomials for parallel shift-by-column operations
     @param[out] res, The output ciphertext
     */
    void shiftBycols_Parallel(Ciphertext& res, Ciphertext& ctxt, long k, ZZX*& shiftpoly);
    
    
    //-------------------------------------------
    // multiplication
    //-------------------------------------------
    
    /**
     @param[out] Initpoly, The polynomials needed for the initial linear transformations to get A[0] and B[0]
     */
    void genMultPoly(ZZX**& Initpoly);
    
    /**
     @param[out] Initpoly, The polynomials needed for the initial linear transformations  
     */
    void genMultPoly_Parallel(ZZX**& Initpoly);
    
    /**
     @param[out] Initpoly, The polynomials needed for the initial linear transformations to get B[0]
     */
    void genMultBPoly(ZZX*& Initpoly);
    
    /**
     @param[in] Actxt, The input ciphertext encrypting a matrix A
     @param[in] Bctxt, The input ciphertext encrypting a matrix B
     @param[in] Initpoly, The polynomials needed for linear transformations of multiplication
     @param[out] resA, The output ciphertext encrypting a permuated matrix A[0]
     @param[out] resB, The output ciphertext encrypting a permuated matrix B[0]
     */
    void genInitCtxt(Ciphertext& resA, Ciphertext& resB, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly);
    
    /**
     @param[in] Actxt, The input ciphertext encrypting multiple matrices As
     @param[in] Bctxt, The input ciphertext encrypting multiple matrices Bs
     @param[in] Initpoly, The polynomials needed for linear transformations of parallel multiplication
     @param[out] resA, The output ciphertext encrypting multiple permuated matrices A[0]'s
     @param[out] resB, The output ciphertext encrypting multiple permuated matrices B[0]'s
     */
    void genInitCtxt_Parallel(Ciphertext& resA, Ciphertext& resB, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly);
    
    /**
     @param[in] Actxt, The input ciphertext encrypting a matrix A
     @param[out] resA, The output ciphertexts encrypting the permuated matrices A[i] for 0 <= i < dim
     */
    void genInitActxt(Ciphertext*& Actxts, Mat<RR>& mat);
    
    /*
     @param[in] Bctxt, The input ciphertext encrypting a matrix B
     @param[in] Initpoly, The polynomials needed for linear transformations of multiplication
     @param[out] resB, The output ciphertext encrypting a permuted matrix B[0]
     */
    void genInitBctxt(Ciphertext& resB, Ciphertext& Bctxt, ZZX*& Initpoly);
    
    /**
     @param[in] Actxt, The input ciphertext encrypting a wide rectangular matrix A
     @param[out] resA, The output ciphertexts encrypting matrices Ai's
     */
    void genInitRecActxt(Ciphertext*& Actxts, Mat<RR>& mat);
    
    /**
     @param[in] Actxts, The input ciphertexts encrypting Ai's
     @param[in] Bctxts, The input ciphertexts encrypting Bi's
     @param[in] num, The number of ciphertexts to multiply
     @param[out] res, The output ciphertext s.t.
     Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1] followed by a rescaling operation
     */
    void HEmatmul_Hadamard(Ciphertext& res, Ciphertext* Actxts, Ciphertext* Bctxts, long num);
    
    /**
     @param[in] Actxt, The input ciphertext encrypting a matrix A
     @param[in] Bctxt, The input ciphertext encrypting a matrix B
     @param[in] Initpoly, The polynomials needed for linear transformations of multiplication
     @param[in] shiftpoly, The polynomials for shift-by-column operations
     @param[out] res, The output ciphertext encrypting a matrix (A * B)
     */
    void HEmatmul(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly);
    
    /**
     @param[in] Actxt, The input ciphertext encrypting multiple matrices As
     @param[in] Bctxt, The input ciphertext encrypting multiple matrices Bs
     @param[in] Initpoly, The polynomials needed for linear transformations of parallel multiplication
     @param[in]  shiftpoly, The polynomials for parallel shift-by-column operations
     @param[out] res, The output ciphertext encrypting matrices (As * Bs)
     */
    void HEmatmul_Parallel(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly);
    
    /**
     @param[in] Actxt, The input ciphertext encrypting a wide rectangular matrix A (an m*n matrix with m<n)
     @param[in] Bctxt, The input ciphertext encrypting a square matrix B (an n*n matrix)
     @param[in] Initpoly, The polynomials needed for linear transformations of multiplication
     @param[in] shiftpoly, The polynomials for shift-by-column operations
     @param[out] res, The output ciphertext encrypting a matrix (A * B)
     */
    void HErmatmul(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly);
    
    /**
     @param[in] Actxts, The input ciphertext encrypting permuted matrices A[i] ("A" are given as fresh ciphertexts)
     @param[in] Bctxt, The input ciphertext encrypting a matrix B
     @param[in] Initpoly, The polynomials needed for linear transformations of multiplication
     @param[out] res, The output ciphertext encrypting a matrix (A * B)
     */
    void HEmatmul_preprocessing(Ciphertext& res, Ciphertext*& Actxts, Ciphertext& Bctxt, ZZX*& Initpoly);
    
    /**
     @param[in] Actxts, The input ciphertext encrypting permuted matrices A[i] (a rectangular matrix "A" are given as fresh ciphertexts)
     @param[in] Bctxt, The input ciphertext encrypting a matrix B
     @param[in] Initpoly, The polynomials needed for linear transformations of multiplication
     @param[out] res, The output ciphertext encrypting a matrix (A * B)
     */
    void HErmatmul_preprocessing(Ciphertext& res, Ciphertext*& Actxts, Ciphertext& Bctxt, ZZX*& Initpoly);
};


#endif
