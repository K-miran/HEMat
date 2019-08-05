//!
//! @file       HEmatrix.cpp
//! @brief      implementing functions for matrix operations
//!
//! @author     Miran Kim
//! @date       Dec. 1, 2017
//! @copyright  GNU Pub License
//!

#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>

#include <cmath>
#include <chrono>
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

#include "../src/Scheme.h"
#include "../src/SecretKey.h"
#include "../src/SchemeAlgo.h"

#include "matrix.h"
#include "HEmatrix.h"


#define timing 0 // use for printing out the timing


//!@ Function: read the parameters
void readHEMatpar(HEMatpar& HEmatpar, long nrows, long ncols, long pBits, long cBits, long logQ, long subdim, long nbatching){
    HEmatpar.nrows = nrows;
    HEmatpar.ncols = ncols;
    
    long dim = nrows;
    if(dim < ncols) {
        dim = ncols;
    }
    
    HEmatpar.dim = (1<< (long)ceil(log2(dim)));  //! power of two
    HEmatpar.dim1 = HEmatpar.dim - 1;
    HEmatpar.logdim = (long) log2(HEmatpar.dim);
    HEmatpar.sqrdim = (long) ceil(sqrt(HEmatpar.dim));
    
    HEmatpar.pBits = pBits;
    HEmatpar.cBits = cBits;
    HEmatpar.logQ = logQ;
    
    HEmatpar.subdim = subdim;  //! used for non-squared matrix multiplication
    HEmatpar.nbatching = (1<< (long)ceil(log2(nbatching)));
    
    if(nbatching == 1){
        HEmatpar.nslots = HEmatpar.dim * HEmatpar.dim;  //! just encode a single matrix
    }
    else{
        HEmatpar.nslots = (HEmatpar.dim * HEmatpar.dim) * HEmatpar.nbatching;
    }
}


void HEmatrix::encryptRmat(Ciphertext& ctxt, mat_RR& mat, long logp){
    complex<double>* cmsg = new complex<double>[HEmatpar.nslots];
    
    NTL_EXEC_RANGE(HEmatpar.nrows, first, last);
    for(int i = first; i < last; ++i){
        for(long j = 0; j < HEmatpar.ncols; ++j){
            double dtemp;
            conv(dtemp, mat[i][j]);
            cmsg[i*HEmatpar.dim + j].real(dtemp);
        }
    }
    NTL_EXEC_RANGE_END;
    
    ctxt = scheme.encrypt(cmsg, HEmatpar.nslots, logp, HEmatpar.logQ);
    delete[] cmsg;
}

void HEmatrix::decryptRmat(mat_RR& mat, Ciphertext& ctxt){
    complex<double>* cmsg = scheme.decrypt(secretKey, ctxt);
    mat.SetDims(HEmatpar.dim, HEmatpar.dim);
    
    long k = 0;
    for(long i = 0; i < HEmatpar.dim; ++i){
        for(long j = 0; j < HEmatpar.dim; ++j){
            mat[i][j] = to_RR(cmsg[k].real());
            k++;
        }
    }
    delete[] cmsg;
}

void HEmatrix::encryptParallelRmat(Ciphertext& ctxt, mat_RR*& mat, long logp, long nbatching){
    complex<double>* cmsg = new complex<double>[HEmatpar.nslots];
    double* dtemp = new double[nbatching];
    
    NTL_EXEC_RANGE(nbatching, first, last);
    for(int k = first; k < last; ++k){
        // encode the kth matrix
        for(long i = 0; i < HEmatpar.nrows; ++i){
            for(long j = 0; j < HEmatpar.ncols; ++j){
                conv(dtemp[k], mat[k][i][j]);
                cmsg[(i*HEmatpar.dim + j)* HEmatpar.nbatching + k].real(dtemp[k]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    ctxt = scheme.encrypt(cmsg, HEmatpar.nslots, logp, HEmatpar.logQ);
    
    delete[] cmsg;
    delete[] dtemp;
}

void HEmatrix::decryptParallelRmat(mat_RR*& mat, Ciphertext& ctxt){
    complex<double>* cmsg = scheme.decrypt(secretKey, ctxt);
    mat = new mat_RR[HEmatpar.nbatching];
    
    for(long k = 0; k < HEmatpar.nbatching; ++k){
        mat[k].SetDims(HEmatpar.dim, HEmatpar.dim);
        long l = 0;
        for(long i = 0; i < HEmatpar.dim; ++i){
            for(long j = 0; j < HEmatpar.dim; ++j){
                mat[k][i][j] = to_RR(cmsg[l * HEmatpar.nbatching + k].real());
                l++;
            }
        }
    }
    delete[] cmsg;
}

void HEmatrix::msgleftRotate(complex<double>*& res, complex<double>* vals, long dim, long nrot){
    long nshift = (nrot)%HEmatpar.nslots;
    
    long k = dim - nshift;
    for(long j = 0; j < k; ++j){
        res[j] = vals[j + nshift];
    }
    for(long j = k; j < dim; ++j){
        res[j] = vals[j - k];
    }
}

void HEmatrix::msgrightRotate(complex<double>*& res, complex<double>* vals, long dim, long nrot){
    long nshift = (nrot)%HEmatpar.nslots;
    long k = dim - nshift;
    
    //! vals[0],vals[1],....,vals[d-nrot-1])
    for(long j = 0; j < k; ++j){
        res[nshift + j] = vals[j];
    }
    
    //! (vals[d-nrot],vals[d-nrot+1],...,vals[d-1])
    for(long j = k; j < dim; ++j){
        res[j - k] = vals[j];
    }
}

void HEmatrix::msgleftRotateAndEqual(complex<double>*& vals, long dim, long nrot){
    complex<double>* res = new complex<double>[dim];

    long nshift = (nrot)%HEmatpar.nslots;
    long k = dim - nshift;
    
    for(long j = 0; j < k; ++j){
        res[j] = vals[j + nshift];
    }
    for(long j = k; j < dim; ++j){
        res[j] = vals[j - k];
    }
    for(long j = 0; j < dim; ++j){
        vals[j] = res[j];
    }
}

void HEmatrix::msgrightRotateAndEqual(complex<double>*& vals, long dim, long nrot){
    complex<double>* res = new complex<double>[dim];
    
    long nshift = (nrot)%HEmatpar.nslots;
    long k = dim - nrot;
    
    for(long j = 0; j < k; ++j){
        res[nshift + j] = vals[j];
    }
    for(long j = k; j < dim; ++j){
        res[j - k] = vals[j];
    }
    for(long j = 0; j < dim; ++j){
        vals[j] = res[j];
    }
}


//------------------------------------------------
//! Transposition
//------------------------------------------------

//! Output: transpoly ("2*dim"), generate the polynomial for transpose
//! Originally we need to generate
//! rho(p[k]; -(d-1)*sqrt(r)i) where k = sqrt(r) * i + j, r = nrows, 0 <= i, j < sqrt(r)
// (0, 1,   ..., d-1),  (-, d+1, ..., 2d-1)

void HEmatrix::genTransPoly(ZZX*& transpoly){
    complex<double>** lvals   = new complex<double>*[HEmatpar.dim];
    complex<double>** rvals  = new complex<double>*[HEmatpar.dim];
    long dsquare = (HEmatpar.nslots) - 1;
    transpoly = new ZZX[2 * HEmatpar.dim];
    long dimsqrdim = HEmatpar.sqrdim * HEmatpar.dim1;   //! (d-1) * sqrt(d)
    
    // k = i * sqrt(d) + j < d
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;   //! all the terms have the same numbers
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(int i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){   //! last term
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            
            lvals[k] = new complex<double>[HEmatpar.nslots];
            rvals[k] = new complex<double>[HEmatpar.nslots];
            
            for(long l = 0; l < HEmatpar.dim - k; ++l){
                long dtemp = l * (HEmatpar.dim + 1) + k;
                lvals[k][dtemp].real(1.0);
                rvals[k][dsquare - dtemp].real(1.0);
            }
            
            msgrightRotateAndEqual(lvals[k], HEmatpar.nslots, i*dimsqrdim);  //! Lrho(P[k],-dsqrt(d)*i)
            msgleftRotateAndEqual(rvals[k], HEmatpar.nslots, i*dimsqrdim);   //! Rrho(P[k],-dsqrt(d)*i)
            
            transpoly[k] = scheme.context.encode(lvals[k], HEmatpar.nslots, HEmatpar.cBits);
            transpoly[k + HEmatpar.dim] = scheme.context.encode(rvals[k], HEmatpar.nslots, HEmatpar.cBits);
        }
    }
    NTL_EXEC_RANGE_END;
    delete[] lvals;
    delete[] rvals;
}

void HEmatrix::transpose(Ciphertext& res, Ciphertext& ctxt, ZZX*& transpoly){
    
    long dimsqrdim = HEmatpar.sqrdim * HEmatpar.dim1;   // (d-1) * sqrt(d)
    
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;   //! all the terms have the same numbers
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    Ciphertext* Babyctxt1 = new Ciphertext[HEmatpar.sqrdim];
    Ciphertext* Babyctxt2 = new Ciphertext[HEmatpar.sqrdim];
    
    Ciphertext** ltemp = new Ciphertext*[HEmatpar.sqrdim];
    Ciphertext** rtemp = new Ciphertext*[HEmatpar.sqrdim];
    
    for(long i = 0; i < HEmatpar.sqrdim; ++i){
        ltemp[i] = new Ciphertext[HEmatpar.sqrdim];
        rtemp[i] = new Ciphertext[HEmatpar.sqrdim];
    }
    
    //! diagonal
    res = scheme.multByPoly(ctxt, transpoly[0], HEmatpar.cBits);
    
    Babyctxt1[0] = ctxt;
    Babyctxt2[0] = ctxt;
    
    //! Babyctxt1[j] = rho(v; j * (d-1))
    //! res[j] = rho(-, (d-1)*srt(d) * i) * rho(v; j * (d-1))
    //HEmatpar.sqrdim
    
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);  // 1 <= j1 <= sqrdim - 1

        Babyctxt1[j1] = scheme.leftRotate(ctxt, j1 * (HEmatpar.dim1));   // prepare baby ctxt
        Babyctxt2[j1] = scheme.rightRotate(ctxt, j1 * (HEmatpar.dim1));

        ltemp[0][j1] = scheme.multByPoly(Babyctxt1[j1], transpoly[j1], HEmatpar.cBits);
        rtemp[0][j1] = scheme.multByPoly(Babyctxt2[j1], transpoly[j1 + HEmatpar.dim], HEmatpar.cBits);

        scheme.addAndEqual(ltemp[0][j1], rtemp[0][j1]);
    }
    NTL_EXEC_RANGE_END;

    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        scheme.addAndEqual(res, ltemp[0][j]);
    }

    //! 1 <= i < d
    //! [k]: rho(p[k];-dsqrt(d)*i1) * rho(v, j*(d-1)) -> rot by dsqrt(d)*i
    //! res[j] = rho( rho(-, (d-1)*srt(d) * i) * rho(v; j * (d-1)) ; i * (d-1)\sqrt(d) )
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if (btmp &&(i == (ibound - 2))){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        long i1 = i + 1;
        long k = (i1) * HEmatpar.sqrdim;

        //! 0 <= j < sqrdim
        NTL_EXEC_RANGE(jbound, first, last);
        for(long j = first; j < last; ++j){
            ltemp[i1][j] = scheme.multByPoly(Babyctxt1[j], transpoly[k + j], HEmatpar.cBits);
            rtemp[i1][j] = scheme.multByPoly(Babyctxt2[j], transpoly[k + j + HEmatpar.dim], HEmatpar.cBits);
        }
        NTL_EXEC_RANGE_END;

        for(long j = 1; j < jbound; ++j){
            scheme.addAndEqual(ltemp[i1][0], ltemp[i1][j]);
            scheme.addAndEqual(rtemp[i1][0], rtemp[i1][j]);
        }

        scheme.leftRotateAndEqual(ltemp[i1][0], i1 * dimsqrdim);
        scheme.rightRotateAndEqual(rtemp[i1][0], i1 * dimsqrdim);

        scheme.addAndEqual(ltemp[i1][0], rtemp[i1][0]);
    }
    NTL_EXEC_RANGE_END;
    
    for(long i = 1; i < ibound; ++i){
        scheme.addAndEqual(res, ltemp[i][0]);
    }
    
    scheme.reScaleByAndEqual(res, HEmatpar.cBits);   // Rescaling by scalarBits
    
    delete[] ltemp;
    delete[] rtemp;
}

void HEmatrix::genTransPoly_Parallel(ZZX*& transpoly){
    
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    complex<double>** lvals   = new complex<double>*[HEmatpar.dim];
    complex<double>** rvals  = new complex<double>*[HEmatpar.dim];
    
    transpoly = new ZZX[2 * HEmatpar.dim];
    
    long shiftunit = HEmatpar.sqrdim * HEmatpar.dim1 * HEmatpar.nbatching ;   //! (d-1) * sqrt(d) * l
    long dsquare = ((HEmatpar.dim * HEmatpar.dim) - 1) *HEmatpar.nbatching ;  //! (d^2-1) * l
    
    //! k = (i * sqrt(d) + j)  < d
    //! poly[k] = poly[i * sqrt(d) + j]  -> rho(poly[k]; i * (d-1)*sqrt(d) * l)
    for(int i = 0; i < ibound; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 0; j < jbound; ++j){
            long k = (i * HEmatpar.sqrdim + j);
        
            lvals[k] = new complex<double>[HEmatpar.nslots];
            rvals[k] = new complex<double>[HEmatpar.nslots];
            
            for(long l = 0; l < HEmatpar.dim - k; ++l){
                long dtemp= (l*(HEmatpar.dim+1) + k) * HEmatpar.nbatching;  //! starting index
        
                for(long n = 0; n < HEmatpar.nbatching; ++n){
                    lvals[k][dtemp + n].real(1.0);
                    rvals[k][dsquare - dtemp + n].real(1.0);
                }
            }
        
            //! Lrho(P[k], - i * (d-1) * sqrt(d) * (nbathcing))
            //! Rrho(P[k],- i * (d-1) * sqrt(d) * (nbathcing))
            msgrightRotateAndEqual(lvals[k], HEmatpar.nslots, i * shiftunit);
            msgleftRotateAndEqual(rvals[k], HEmatpar.nslots, i * shiftunit);
            
            transpoly[k] = scheme.context.encode(lvals[k], HEmatpar.nslots, HEmatpar.cBits);
            transpoly[k + HEmatpar.dim] = scheme.context.encode(rvals[k], HEmatpar.nslots, HEmatpar.cBits);
            
        }
    }
    delete[] lvals;
    delete[] rvals;
}

void HEmatrix::transpose_Parallel(Ciphertext& res, Ciphertext& ctxt, ZZX*& transpoly){
    long shiftunit = HEmatpar.sqrdim * HEmatpar.dim1 * HEmatpar.nbatching;  // (d-1) * sqrt(d) * l
    long unit  = (HEmatpar.dim1) * HEmatpar.nbatching;
    
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    Ciphertext* Babyctxt1 = new Ciphertext[HEmatpar.sqrdim];
    Ciphertext* Babyctxt2 = new Ciphertext[HEmatpar.sqrdim];
    
    Ciphertext** ltemp = new Ciphertext*[HEmatpar.sqrdim];
    Ciphertext** rtemp = new Ciphertext*[HEmatpar.sqrdim];
    
    for(long i = 0; i < HEmatpar.sqrdim; ++i){
        ltemp[i] = new Ciphertext[HEmatpar.sqrdim];
        rtemp[i] = new Ciphertext[HEmatpar.sqrdim];
    }
    
    //-------------------------
    // i = 0, sqrt(d) polynomials
    //-------------------------
    res = scheme.multByPoly(ctxt, transpoly[0], HEmatpar.cBits);

    Babyctxt1[0] = ctxt;
    Babyctxt2[0] = ctxt;
    
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);  //! 1 <= j1 <= sqrdim - 1
        
        // rho(v; j * (d-1) * l)
        Babyctxt1[j1] = scheme.leftRotate(ctxt, j1 * unit);
        Babyctxt2[j1] = scheme.rightRotate(ctxt, j1 * unit);
        
        ltemp[0][j1] = scheme.multByPoly(Babyctxt1[j1], transpoly[j1], HEmatpar.cBits);
        rtemp[0][j1] = scheme.multByPoly(Babyctxt2[j1], transpoly[j1 + HEmatpar.dim], HEmatpar.cBits);
        
        scheme.addAndEqual(ltemp[0][j1], rtemp[0][j1]);
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        scheme.addAndEqual(res, ltemp[0][j]);
    }
    
    //! 1 <= i < d
    //! [k]: rho(p[k];-dsqrt(d)*i1) * rho(v, j*(d-1)) -> rot by dsqrt(d)*i
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        long i1 = i + 1;
        long k = (i1) * HEmatpar.sqrdim;
        
        NTL_EXEC_RANGE(jbound, first, last);
        for(long j = first; j < last; ++j){
            ltemp[i1][j] = scheme.multByPoly(Babyctxt1[j], transpoly[k + j], HEmatpar.cBits);
            rtemp[i1][j] = scheme.multByPoly(Babyctxt2[j], transpoly[k + j + HEmatpar.dim], HEmatpar.cBits);
        }
        NTL_EXEC_RANGE_END;
        
        for(long j = 1; j < jbound; ++j){
            scheme.addAndEqual(ltemp[i1][0], ltemp[i1][j]);
            scheme.addAndEqual(rtemp[i1][0], rtemp[i1][j]);
        }
        
        scheme.leftRotateAndEqual(ltemp[i1][0], i1 * shiftunit);
        scheme.rightRotateAndEqual(rtemp[i1][0], i1 * shiftunit);
        
        scheme.addAndEqual(ltemp[i1][0], rtemp[i1][0]);
    }
    NTL_EXEC_RANGE_END;
    
    for(long i = 1; i < ibound; ++i){
        scheme.addAndEqual(res, ltemp[i][0]);
    }

    scheme.reScaleByAndEqual(res, HEmatpar.cBits);   //! Rescaling by scalarBits
    
    delete[] ltemp;
    delete[] rtemp;
}


//------------------------------------------------
//! Shift
//------------------------------------------------

// num = 0: dimension d-1
// subdim - 1: dimension (subdim-1)
void HEmatrix::genShiftPoly(ZZX*& shiftpoly, long num){
    long length = HEmatpar.dim1;
    if(num != 0){
        length = num;
    }
    
    complex<double>** vals = new complex<double>*[length];
    shiftpoly = new ZZX[length];
    
    // i: shifted by (i+1)
    NTL_EXEC_RANGE(length, first, last);
    for(int i = first; i < last; ++i){
        vals[i] = new complex<double>[HEmatpar.nslots];
        
        for(long j = 0; j < HEmatpar.dim; ++j){
            for(long k = 0; k < i + 1; ++k){
                vals[i][j * HEmatpar.dim + k].real(1.0);
            }
        }
        shiftpoly[i] = scheme.context.encode(vals[i], HEmatpar.nslots, HEmatpar.cBits);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] vals;
}

void HEmatrix::shiftBycols(Ciphertext& res, Ciphertext& ctxt,  long k , ZZX*& shiftpoly){
    //! ctemp = Enc(m[1],.., m[k],0,,,0 | ....)
    Ciphertext ctemp = scheme.multByPoly(ctxt, shiftpoly[k-1], HEmatpar.cBits);
    scheme.reScaleByAndEqual(ctemp, HEmatpar.cBits);
    
    //! ctemp = Enc(0,,,0, m[k+1],...,m[d] | ....)
    res = scheme.modDownTo(ctxt, ctemp.logq);
    scheme.subAndEqual(res, ctemp);
    
    scheme.rightRotateAndEqual(ctemp, HEmatpar.dim - k);
    scheme.leftRotateAndEqual(res, k);
    scheme.addAndEqual(res, ctemp);
}

void HEmatrix::genShiftPoly_Parallel(ZZX*& shiftpoly){
    complex<double>** vals   = new complex<double>*[HEmatpar.dim1];
    shiftpoly = new ZZX[HEmatpar.dim1];
    
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(int i = first; i < last; ++i){
        vals[i] = new complex<double>[HEmatpar.nslots];
        
        for(long j = 0; j < HEmatpar.dim; ++j){
            for(long k = 0; k < i + 1; ++k){
                long dtemp = (j * HEmatpar.dim + k) * HEmatpar.nbatching;
                
                for(long n = 0; n < HEmatpar.nbatching; ++n){
                    vals[i][dtemp + n].real(1.0);
                }
            }
        }
        shiftpoly[i] = scheme.context.encode(vals[i], HEmatpar.nslots, HEmatpar.cBits);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] vals;
}

void HEmatrix::shiftBycols_Parallel(Ciphertext& res, Ciphertext& ctxt, long k, ZZX*& shiftpoly){
    //! ctemp = Enc(m[1],.., m[k],0,,,0 | ....)
    Ciphertext ctemp = scheme.multByPoly(ctxt, shiftpoly[k-1], HEmatpar.cBits);
    scheme.reScaleByAndEqual(ctemp, HEmatpar.cBits);
    
    //! ctemp = Enc(0,,,0, m[k+1],...,m[d] | ....)
    res = scheme.modDownTo(ctxt, ctemp.logq);
    scheme.subAndEqual(res, ctemp);
    
    scheme.rightRotateAndEqual(ctemp, (HEmatpar.dim - k) * HEmatpar.nbatching);
    scheme.leftRotateAndEqual(res, k * HEmatpar.nbatching);
    scheme.addAndEqual(res, ctemp);
}

//------------------------------------------------
//! Matrix multiplication
//------------------------------------------------

//!@ Output: Initpoly[0] (constant left-polynomials for Amat): 0,1,....,d-1
//!@         Initpoly[1] (constant right-polynomials for Amat):  -1,...,-(d-1)
//!!         Initpoly[2] (constant polynomials for Bmat)

void HEmatrix::genMultPoly(ZZX**& Initpoly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double) HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    Initpoly = new ZZX*[3];
    
    Initpoly[0] =  new ZZX[HEmatpar.dim];
    Initpoly[1] =  new ZZX[HEmatpar.dim];
    Initpoly[2] =  new ZZX[HEmatpar.dim];
    
    complex<double>** fvals1   = new complex<double>*[HEmatpar.dim];
    complex<double>** fvals2   = new complex<double>*[HEmatpar.dim];
    complex<double>** bvals    = new complex<double>*[HEmatpar.dim];

    NTL_EXEC_RANGE(ibound, first, last);
    for(int i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            fvals1[k] = new complex<double>[HEmatpar.nslots];
            fvals2[k] = new complex<double>[HEmatpar.nslots];
            bvals[k] = new complex<double>[HEmatpar.nslots];
            
            for(long l = 0; l < HEmatpar.dim - k; ++l){
                fvals1[k][k * HEmatpar.dim + l].real(1.0);
            }
            msgleftRotate(fvals2[k], fvals1[k], HEmatpar.nslots, k*(2*HEmatpar.dim - 1));
            
            msgrightRotateAndEqual(fvals1[k], HEmatpar.nslots, i*HEmatpar.sqrdim);
            Initpoly[0][k] = scheme.context.encode(fvals1[k], HEmatpar.nslots, HEmatpar.cBits);
        
            msgleftRotateAndEqual(fvals2[k], HEmatpar.nslots, i*HEmatpar.sqrdim);
            Initpoly[1][k] = scheme.context.encode(fvals2[k], HEmatpar.nslots, HEmatpar.cBits);
            
            for(long l = 0; l < HEmatpar.dim; ++l){
                bvals[k][l * HEmatpar.dim + k].real(1.0);
            }
            msgrightRotateAndEqual(bvals[k], HEmatpar.nslots, i*HEmatpar.sqrdim*HEmatpar.dim);
            Initpoly[2][k] = scheme.context.encode(bvals[k], HEmatpar.nslots, HEmatpar.cBits);
        }
    }
    NTL_EXEC_RANGE_END;

    delete[] fvals1;
    delete[] fvals2;
    delete[] bvals;
}

void HEmatrix::genInitCtxt(Ciphertext& resA, Ciphertext& resB, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& poly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    Ciphertext** Actemp1 = new Ciphertext*[HEmatpar.sqrdim]; //! update right polynomial
    Ciphertext** Actemp2 = new Ciphertext*[HEmatpar.sqrdim]; //! update right polynomial
    Ciphertext** Bctemp = new Ciphertext*[HEmatpar.sqrdim];
    
    for(long i = 0; i < HEmatpar.sqrdim; ++i){
        Actemp1[i] = new Ciphertext[HEmatpar.sqrdim];
        Actemp2[i] = new Ciphertext[HEmatpar.sqrdim];
        Bctemp[i]  = new Ciphertext[HEmatpar.sqrdim];
    }
    
    //! 0. Store some ciphertexts (0,1,...,d-1), (,d+1,...2d-1)
    //! v, lrho(v;1), lrho(v;2), ..., lrho(v;d-1)
    //! -, rrho(v;1), rrho(v;2), ..., rrho(v;d-1)
    
    Ciphertext* BaByctxt1  = new Ciphertext[HEmatpar.dim];
    Ciphertext* BaByctxt2  = new Ciphertext[HEmatpar.dim];
    Ciphertext* BaByctxtB  = new Ciphertext[HEmatpar.dim];
    
    BaByctxt1[0] = Actxt;
    BaByctxt2[0] = Actxt;
    BaByctxtB[0] = Bctxt;
    
    //! i = 0:   Actxts[0] = v[0] + p1 * v[1] + ... + p[sqr(d)-1] *  v[sqr(d)-1]
    resA = scheme.multByPoly(BaByctxt1[0], poly[0][0], HEmatpar.cBits);
    resB = scheme.multByPoly(BaByctxtB[0], poly[2][0], HEmatpar.cBits);
    
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxt1[j1] = scheme.leftRotate(Actxt, j1);
        BaByctxt2[j1] = scheme.rightRotate(Actxt, j1);
        BaByctxtB[j1] = scheme.leftRotate(Bctxt, (j1)*HEmatpar.dim);
        
        Actemp1[0][j1] = scheme.multByPoly(BaByctxt1[j1], poly[0][j1], HEmatpar.cBits);
        Actemp2[0][j1] = scheme.multByPoly(BaByctxt2[j1], poly[1][j1], HEmatpar.cBits);
        Bctemp[0][j1]  = scheme.multByPoly(BaByctxtB[j1], poly[2][j1], HEmatpar.cBits);
        
        scheme.addAndEqual(Actemp1[0][j1], Actemp2[0][j1]);
    }
    NTL_EXEC_RANGE_END;

    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        scheme.addAndEqual(resA, Actemp1[0][j]);
        scheme.addAndEqual(resB, Bctemp[0][j]);
    }
    
    Ciphertext* Actxts1 = new Ciphertext[HEmatpar.dim];
    Ciphertext* Actxts2 = new Ciphertext[HEmatpar.dim];
    Ciphertext* Bctxts  = new Ciphertext[HEmatpar.dim];
    
    NTL_EXEC_RANGE(HEmatpar.dim - HEmatpar.sqrdim, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + HEmatpar.sqrdim);
        long i = (long)(k1 / HEmatpar.sqrdim);
        long j = (long)(k1 % HEmatpar.sqrdim);
        
        Actxts1[k1] = scheme.multByPoly(BaByctxt1[j], poly[0][k1], HEmatpar.cBits);
        Actxts2[k1] = scheme.multByPoly(BaByctxt2[j], poly[1][k1], HEmatpar.cBits);
        Bctxts[k1]  = scheme.multByPoly(BaByctxtB[j], poly[2][k1], HEmatpar.cBits);
        
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k+1) * HEmatpar.sqrdim;
        long jbound = HEmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 1; j < jbound; ++j){
            scheme.addAndEqual(Actxts1[k1], Actxts1[k1+j]);
            scheme.addAndEqual(Actxts2[k1], Actxts2[k1+j]);
            scheme.addAndEqual(Bctxts[k1], Bctxts[k1+j]);
        }
        
        long k2 = (k1 - (k1 % HEmatpar.sqrdim));
        scheme.leftRotateAndEqual(Actxts1[k1], k2);
        scheme.rightRotateAndEqual(Actxts2[k1], k2);
        scheme.leftRotateAndEqual(Bctxts[k1], k2 * HEmatpar.dim);
        
        scheme.addAndEqual(Actxts1[k1], Actxts2[k1]);
    }
    NTL_EXEC_RANGE_END;
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * HEmatpar.sqrdim;
        scheme.addAndEqual(resA, Actxts1[k1]);
        scheme.addAndEqual(resB, Bctxts[k1]);
    }
    
    scheme.reScaleByAndEqual(resA, HEmatpar.cBits);
    scheme.reScaleByAndEqual(resB, HEmatpar.cBits);

    delete[] Actxts1;
    delete[] Actxts2;
    delete[] Bctxts;
    
    delete[] Actemp1;
    delete[] Actemp2;
    delete[] Bctemp;
    
    delete[] BaByctxt1;
    delete[] BaByctxt2;
    delete[] BaByctxtB;
}

void HEmatrix::HEmatmul_Hadamard(Ciphertext& res, Ciphertext* Actxts, Ciphertext* Bctxts, long num){
    //! logq: Actxt > Actxts[0] = Actxts[i + 1] + 2*cBits
    //! logq: Bctxt > Bctxts[0] = Bctxt[i + 1]
    
    long num1 = num - 1;
    //! 1) Bctxts[0] < Actxts[1] < Actxts[0]
    if(Actxts[1].logq > Bctxts[0].logq){
        res = Actxts[0];
        scheme.modDownToAndEqual(res, Bctxts[0].logq);
        scheme.multAndEqual(res, Bctxts[0]);   //! log(res) = log Bctxt[0]
        
        NTL_EXEC_RANGE(num1, first, last);
        for(int i = first; i < last; ++i){
            long i1 = (i + 1);
            scheme.modDownToAndEqual(Actxts[i1], Bctxts[0].logq);
            scheme.multAndEqual(Actxts[i1], Bctxts[i1]);   //! log(Actxt[i+1]) = log Bctxt[0]
        }
        NTL_EXEC_RANGE_END;
    }
    
    //! 2)  Actxts[1] < Actxts[0] < Bctxts[0]
    else if(Actxts[0].logq < Bctxts[0].logq){
        res = Bctxts[0];
        scheme.modDownToAndEqual(res, Actxts[0].logq);
        scheme.multAndEqual(res, Actxts[0]);   //! log(res) = log Actxt[0]
        
        NTL_EXEC_RANGE(num1, first, last);
        for(int i = first; i < last; ++i){
            long i1 = (i + 1);
            scheme.modDownToAndEqual(Bctxts[i1], Actxts[i1].logq);
            scheme.multAndEqual(Actxts[i1], Bctxts[i1]);
        }
        NTL_EXEC_RANGE_END;
        
        scheme.modDownToAndEqual(res, Actxts[1].logq);
    }
    
    //! 3)  Actxts[1] < Bctxts[i] < Actxts[0]: the common case
    else{
        res = Actxts[0];
        scheme.modDownToAndEqual(res, Bctxts[0].logq);
        scheme.multAndEqual(res, Bctxts[0]);   //! log(res) = log Bctxt[0]
        
        NTL_EXEC_RANGE(num1, first, last);
        for(int i = first; i < last; ++i){
            long i1 = (i + 1);
            scheme.modDownToAndEqual(Bctxts[i1], Actxts[i1].logq);
            scheme.multAndEqual(Actxts[i1], Bctxts[i1]);
        }
        NTL_EXEC_RANGE_END;
        
        scheme.modDownToAndEqual(res, Actxts[1].logq);
    }
    
    //! aggregate the results and rescaling
    for(int i = 1; i < num; ++i){
        scheme.addAndEqual(res, Actxts[i]);
    }
    scheme.reScaleByAndEqual(res, res.logp);
  
}

void HEmatrix::HEmatmul(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly){
    Ciphertext* Actxts = new Ciphertext[HEmatpar.dim];
    Ciphertext* Bctxts = new Ciphertext[HEmatpar.dim];
   
    //! 1. Generate the initial ciphertexts
    genInitCtxt(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);

    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(int i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = scheme.leftRotate(Bctxts[0], HEmatpar.dim * (i1));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard multiplication
    HEmatmul_Hadamard(res, Actxts, Bctxts, HEmatpar.dim);
    
    delete[] Actxts;
    delete[] Bctxts;
}

void HEmatrix::genMultPoly_Parallel(ZZX**& Initpoly){
    long btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"

    
    Initpoly = new ZZX*[3];
    
    Initpoly[0] =  new ZZX[HEmatpar.dim];
    Initpoly[1] =  new ZZX[HEmatpar.dim];
    Initpoly[2] =  new ZZX[HEmatpar.dim];
    
    complex<double>** fvals1   = new complex<double>*[HEmatpar.dim];
    complex<double>** fvals2   = new complex<double>*[HEmatpar.dim];
    complex<double>** bvals    = new complex<double>*[HEmatpar.dim];
    
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(int i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp) && (i == (ibound - 1))){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            fvals1[k] = new complex<double>[HEmatpar.nslots];
            fvals2[k] = new complex<double>[HEmatpar.nslots];
            bvals[k] = new complex<double>[HEmatpar.nslots];
            
            //! original: k * d <= index <= k * d + d - k -1
            //! parallel: (k * d) * l + 0 <= ... <= (k * d) * l + (l-1)
            //!           (k * d + d - k -1) * l + 0 <= ... <= (k * d + d - k -1) * l + (l-1) <  (k * d + d - k -1) * l + l
            //!           from k * d * l <= index <  (k * d + d - k) * l
            long start = k * HEmatpar.dim * HEmatpar.nbatching;
            long end = ( k * HEmatpar.dim + HEmatpar.dim - k) * HEmatpar.nbatching;
            for(long l = start; l < end; ++l){
                fvals1[k][l].real(1.0);
            }
            
            msgleftRotate(fvals2[k], fvals1[k], HEmatpar.nslots, k*(2*HEmatpar.dim - 1) * HEmatpar.nbatching);
            
            msgrightRotateAndEqual(fvals1[k], HEmatpar.nslots, i*HEmatpar.sqrdim * HEmatpar.nbatching);
            Initpoly[0][k] = scheme.context.encode(fvals1[k], HEmatpar.nslots, HEmatpar.cBits);
            
            msgleftRotateAndEqual(fvals2[k], HEmatpar.nslots, i*HEmatpar.sqrdim * HEmatpar.nbatching);
            Initpoly[1][k] = scheme.context.encode(fvals2[k], HEmatpar.nslots, HEmatpar.cBits);
            
            for(long l = 0; l < HEmatpar.dim; ++l){
                long dtemp = (l * HEmatpar.dim + k) * HEmatpar.nbatching;
                
                for(long n = 0; n < HEmatpar.nbatching; ++n){
                    bvals[k][dtemp + n].real(1.0);
                }
            }
            msgrightRotateAndEqual(bvals[k], HEmatpar.nslots, i*HEmatpar.sqrdim*HEmatpar.dim * HEmatpar.nbatching);
            Initpoly[2][k] = scheme.context.encode(bvals[k], HEmatpar.nslots, HEmatpar.cBits);
        }
    }
    NTL_EXEC_RANGE_END;
    
    delete[] fvals1;
    delete[] fvals2;
    delete[] bvals;
}

void HEmatrix::genInitCtxt_Parallel(Ciphertext& resA, Ciphertext& resB, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& poly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"

    
    Ciphertext** Actemp1 = new Ciphertext*[HEmatpar.sqrdim]; //! update right polynomial
    Ciphertext** Actemp2 = new Ciphertext*[HEmatpar.sqrdim]; //! update right polynomial
    Ciphertext** Bctemp = new Ciphertext*[HEmatpar.sqrdim];
    
    for(long i = 0; i < HEmatpar.sqrdim; ++i){
        Actemp1[i] = new Ciphertext[HEmatpar.sqrdim];
        Actemp2[i] = new Ciphertext[HEmatpar.sqrdim];
        Bctemp[i]  = new Ciphertext[HEmatpar.sqrdim];
    }

    //---------------------------------------------------------------
    // 0. Store some ciphertexts (0,1,...,d-1), ( ,d+1,...2d-1)
    // v, lrho(v;1), lrho(v;2), ..., lrho(v;d-1)
    // -, rrho(v;1), rrho(v;2), ..., rrho(v;d-1)
    
    Ciphertext* BaByctxt1  = new Ciphertext[HEmatpar.dim];
    Ciphertext* BaByctxt2  = new Ciphertext[HEmatpar.dim];
    Ciphertext* BaByctxtB  = new Ciphertext[HEmatpar.dim];
    
    BaByctxt1[0] = Actxt;
    BaByctxt2[0] = Actxt;
    BaByctxtB[0] = Bctxt;
    
    resA = scheme.multByPoly(BaByctxt1[0], poly[0][0], HEmatpar.cBits);
    resB = scheme.multByPoly(BaByctxtB[0], poly[2][0], HEmatpar.cBits);
    
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxt1[j1] = scheme.leftRotate(Actxt, j1 * HEmatpar.nbatching);
        BaByctxt2[j1] = scheme.rightRotate(Actxt, j1 * HEmatpar.nbatching);
        BaByctxtB[j1] = scheme.leftRotate(Bctxt, (j1)*HEmatpar.dim * HEmatpar.nbatching);
        
        Actemp1[0][j1] = scheme.multByPoly(BaByctxt1[j1], poly[0][j1], HEmatpar.cBits);
        Actemp2[0][j1] = scheme.multByPoly(BaByctxt2[j1], poly[1][j1], HEmatpar.cBits);
        Bctemp[0][j1]  = scheme.multByPoly(BaByctxtB[j1], poly[2][j1], HEmatpar.cBits);
        
        scheme.addAndEqual(Actemp1[0][j1], Actemp2[0][j1]);
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        scheme.addAndEqual(resA, Actemp1[0][j]);
        scheme.addAndEqual(resB, Bctemp[0][j]);
    }
 
    //---------------------------
    Ciphertext* Actxts1 = new Ciphertext[HEmatpar.dim];
    Ciphertext* Actxts2 = new Ciphertext[HEmatpar.dim];
    Ciphertext* Bctxts  = new Ciphertext[HEmatpar.dim];
    
    NTL_EXEC_RANGE(HEmatpar.dim - ibound, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + HEmatpar.sqrdim);
        long i = (long)(k1 / HEmatpar.sqrdim);
        long j = (long)(k1 % HEmatpar.sqrdim);
        
        Actxts1[k1] = scheme.multByPoly(BaByctxt1[j], poly[0][k1], HEmatpar.cBits);
        Actxts2[k1] = scheme.multByPoly(BaByctxt2[j], poly[1][k1], HEmatpar.cBits);
        Bctxts[k1]  = scheme.multByPoly(BaByctxtB[j], poly[2][k1], HEmatpar.cBits);
        
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k+1) * HEmatpar.sqrdim;
        
        long jbound = HEmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 1; j < jbound; ++j){
            scheme.addAndEqual(Actxts1[k1], Actxts1[k1+j]);
            scheme.addAndEqual(Actxts2[k1], Actxts2[k1+j]);
            scheme.addAndEqual(Bctxts[k1], Bctxts[k1+j]);
        }
        
        long k2 = (k1 - (k1 % HEmatpar.sqrdim)) * HEmatpar.nbatching;
        scheme.leftRotateAndEqual(Actxts1[k1], k2);
        scheme.rightRotateAndEqual(Actxts2[k1], k2);
        scheme.leftRotateAndEqual(Bctxts[k1], k2 * HEmatpar.dim);
        
        scheme.addAndEqual(Actxts1[k1], Actxts2[k1]);
    }
    NTL_EXEC_RANGE_END;
    
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * HEmatpar.sqrdim;
        scheme.addAndEqual(resA, Actxts1[k1]);
        scheme.addAndEqual(resB, Bctxts[k1]);
    }
    
    scheme.reScaleByAndEqual(resA, HEmatpar.cBits);
    scheme.reScaleByAndEqual(resB, HEmatpar.cBits);
    
    delete[] Actxts1;
    delete[] Actxts2;
    delete[] Bctxts;
    
    delete[] Actemp1;
    delete[] Actemp2;
    delete[] Bctemp;
    
    delete[] BaByctxt1;
    delete[] BaByctxt2;
    delete[] BaByctxtB;
}

void HEmatrix::HEmatmul_Parallel(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly){
    Ciphertext* Actxts = new Ciphertext[HEmatpar.dim];
    Ciphertext* Bctxts = new Ciphertext[HEmatpar.dim];
    
    //! 1. Generate the initial ciphertexts
    genInitCtxt_Parallel(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);

    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    long unit = HEmatpar.dim  * HEmatpar.nbatching;
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(int i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols_Parallel(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = scheme.leftRotate(Bctxts[0], unit * (i1));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    HEmatmul_Hadamard(res, Actxts, Bctxts, HEmatpar.dim);

    delete[] Actxts;
    delete[] Bctxts;
}

void HEmatrix::HErmatmul(Ciphertext& res, Ciphertext& Actxt, Ciphertext& Bctxt, ZZX**& Initpoly, ZZX*& shiftpoly){
    Ciphertext* Actxts = new Ciphertext[HEmatpar.subdim];
    Ciphertext* Bctxts = new Ciphertext[HEmatpar.subdim];
    
    //! 1. Generate the initial ciphertexts
    genInitCtxt(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);
    
    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.subdim - 1, first, last);
    for(int i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = scheme.leftRotate(Bctxts[0], HEmatpar.dim * (i1));
    }
    NTL_EXEC_RANGE_END;
 
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    HEmatmul_Hadamard(res, Actxts, Bctxts, HEmatpar.subdim);
    
    //! 4. shift and aggregate the results
    long index = (long) log2(HEmatpar.dim/HEmatpar.subdim);
    
    for(long i = 0; i < index; ++i){
        Ciphertext ctemp = scheme.leftRotate(res, HEmatpar.dim*HEmatpar.subdim * (1<<i));
        scheme.addAndEqual(res, ctemp);
    }

    delete[] Actxts;
    delete[] Bctxts;
}

void HEmatrix::genMultBPoly(ZZX*& Initpoly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    Initpoly =  new ZZX[HEmatpar.dim];
    complex<double>** bvals  = new complex<double>*[HEmatpar.dim];
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(int i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            bvals[k] = new complex<double>[HEmatpar.nslots];
            
            for(long l = 0; l < HEmatpar.dim; ++l){
                bvals[k][l * HEmatpar.dim + k].real(1.0);
            }
            msgrightRotateAndEqual(bvals[k], HEmatpar.nslots, i*HEmatpar.sqrdim*HEmatpar.dim);
            Initpoly[k] = scheme.context.encode(bvals[k], HEmatpar.nslots, HEmatpar.cBits);
        }
    }
    NTL_EXEC_RANGE_END;
    
    delete[] bvals;
}

void HEmatrix::genInitActxt(Ciphertext*& Actxts, Mat<RR>& mat){
    Mat<RR>* Amat = new Mat<RR>[HEmatpar.dim];
    Actxts = new Ciphertext[HEmatpar.dim];
    complex<double>** cmsg = new complex<double>*[HEmatpar.dim];

    NTL_EXEC_RANGE(HEmatpar.dim, first, last);
    for(long k = first; k < last; ++k){
        Amat[k].SetDims(HEmatpar.dim, HEmatpar.dim);
        long dimk = HEmatpar.dim - k;
    
        //! 0 <= i < d - k: shift by (k+i)-positions from mat[i]
        for(long i = 0; i < dimk; ++i){
            long nshift = k + i;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = mat[i][j + nshift];
            }
            
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = mat[i][j - nshift2];
            }
        }
        
        //! i = d - k : Amat[k][i] <- mat[i]
        if(k!=0){
            for(long j = 0; j < HEmatpar.dim; ++j){
                Amat[k][dimk][j] = mat[dimk][j];
            }
        }
        
        //! d - k + 1 <= i < d: shift by (k+i-d)-positions from mat[i]
        for(long i = dimk + 1; i < HEmatpar.dim; ++i){
            long nshift =  i - dimk;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = mat[i][j + nshift];
            }
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = mat[i][j - nshift2];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encryption of d Rmat
    NTL_EXEC_RANGE(HEmatpar.dim, first, last);
    for(long k = first; k < last; ++k){
        cmsg[k] = new complex<double>[HEmatpar.nslots];
    
        for(long i = 0; i < HEmatpar.nrows; ++i){
            for(long j = 0; j < HEmatpar.ncols; ++j){
                double dtemp;
                conv(dtemp, Amat[k][i][j]);
                cmsg[k][i*HEmatpar.dim + j].real(dtemp);
            }
        }
        Actxts[k] = scheme.encrypt(cmsg[k], HEmatpar.nslots, HEmatpar.pBits, HEmatpar.logQ);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] cmsg;
}

void HEmatrix::genInitBctxt(Ciphertext& resB, Ciphertext& Bctxt, ZZX*& poly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim);
    
    
    Ciphertext** Bctemp = new Ciphertext*[HEmatpar.sqrdim];
    for(long i = 0; i < HEmatpar.sqrdim; ++i){
        Bctemp[i]  = new Ciphertext[HEmatpar.sqrdim];
    }
    
    //! 0. Store some ciphertexts (0,1,...,d-1), ( ,d+1,...2d-1)
    Ciphertext* BaByctxtB  = new Ciphertext[HEmatpar.dim];
    
    BaByctxtB[0] = Bctxt;
    
    //! i = 0:   Actxts[0] = v[0] + p1 * v[1] + ... + p[sqr(d)-1] *  v[sqr(d)-1]
    resB = scheme.multByPoly(BaByctxtB[0], poly[0], HEmatpar.cBits);

    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxtB[j1] = scheme.leftRotate(Bctxt, (j1)*HEmatpar.dim);
        Bctemp[0][j1]  = scheme.multByPoly(BaByctxtB[j1], poly[j1], HEmatpar.cBits);
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        scheme.addAndEqual(resB, Bctemp[0][j]);
    }
    
    Ciphertext* Bctxts  = new Ciphertext[HEmatpar.dim];
    
    NTL_EXEC_RANGE(HEmatpar.dim - HEmatpar.sqrdim, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + HEmatpar.sqrdim);
        long i = (long)(k1 / HEmatpar.sqrdim);
        long j = (long)(k1 % HEmatpar.sqrdim);
        
        Bctxts[k1]  = scheme.multByPoly(BaByctxtB[j], poly[k1], HEmatpar.cBits);
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k+1) * HEmatpar.sqrdim;
        long jbound = HEmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 1; j < jbound; ++j){
            scheme.addAndEqual(Bctxts[k1], Bctxts[k1+j]);
        }
        
        long k2 = (k1 - (k1 % HEmatpar.sqrdim));
        scheme.leftRotateAndEqual(Bctxts[k1], k2 * HEmatpar.dim);
    }
    NTL_EXEC_RANGE_END;
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * HEmatpar.sqrdim;
        scheme.addAndEqual(resB, Bctxts[k1]);
    }

    scheme.reScaleByAndEqual(resB, HEmatpar.cBits);
    
    delete[] Bctxts;
    delete[] Bctemp;
    delete[] BaByctxtB;
}

void HEmatrix::HEmatmul_preprocessing(Ciphertext& res, Ciphertext*& Actxts, Ciphertext& Bctxt, ZZX*& Initpoly){

    //! 1. Generate the initial ciphertexts
    Ciphertext* Bctxts = new Ciphertext[HEmatpar.dim];
    genInitBctxt(Bctxts[0], Bctxt, Initpoly);
    
    //! 2. Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(int i = first; i < last; ++i){
        long i1 = (i + 1);
        Bctxts[i1] = scheme.leftRotate(Bctxts[0], HEmatpar.dim * (i1));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    NTL_EXEC_RANGE(HEmatpar.dim, first, last);
    for(int i = first; i < last; ++i){
        scheme.modDownToAndEqual(Actxts[i], Bctxts[0].logq);
        scheme.multAndEqual(Actxts[i], Bctxts[i]);   // log(Actxt[i+1]) = log Bctxt[0]
    }
    NTL_EXEC_RANGE_END;
    
    res = Actxts[0];
    for(int i = 1; i < HEmatpar.dim; ++i){
        scheme.addAndEqual(res, Actxts[i]);
    }
    
    scheme.reScaleByAndEqual(res, res.logp);

    delete[] Bctxts;
}

void HEmatrix::genInitRecActxt(Ciphertext*& Actxts, Mat<RR>& mat){
    // rep_mat: (mat; mat; ... ;mat) square mat
    Mat<RR> replicate_mat;
    replicate_mat.SetDims(HEmatpar.dim, HEmatpar.dim);
    
    long index_rows = HEmatpar.dim/HEmatpar.subdim;
    
    NTL_EXEC_RANGE(index_rows, first, last);
    for(long k = first; k < last; ++k){
        for(long i = 0; i < HEmatpar.subdim; ++i){
            for(long j = 0; j < HEmatpar.dim; ++j){
                replicate_mat[k*HEmatpar.subdim + i][j] = mat[i][j];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! generate the (linear transformed) matrices
    Mat<RR>* Amat = new Mat<RR>[HEmatpar.subdim];
    Actxts = new Ciphertext[HEmatpar.subdim];
    
    NTL_EXEC_RANGE(HEmatpar.subdim, first, last);
    for(long k = first; k < last; ++k){
        Amat[k].SetDims(HEmatpar.dim, HEmatpar.dim);
        long dimk= HEmatpar.dim - k;
        
        //! 0 <= i < d - k: shift by (k+i)-positions from mat[i]
        for(long i = 0; i < dimk; ++i){
            long nshift = k + i;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = replicate_mat[i][j + nshift];
            }
            
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = replicate_mat[i][j - nshift2];
            }
        }
        
        //! i = d - k : Amat[k][i] <- mat[i]
        if(k!=0){
            for(long j = 0; j < HEmatpar.dim; ++j){
                Amat[k][dimk][j] = replicate_mat[dimk][j];
            }
        }
        
        //! d - k + 1 <= i < d: shift by (k+i-d)-positions from mat[i]
        for(long i = dimk + 1; i < HEmatpar.dim; ++i){
            long nshift =  i - dimk;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = replicate_mat[i][j + nshift];
            }
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = replicate_mat[i][j - nshift2];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encryption of d Rmat
    complex<double>** cmsg = new complex<double>*[HEmatpar.subdim];
    NTL_EXEC_RANGE(HEmatpar.subdim, first, last);
    for(long k = first; k < last; ++k){
        cmsg[k] = new complex<double>[HEmatpar.nslots];
        for(int i = 0; i < HEmatpar.nrows; ++i){
            for(long j = 0; j < HEmatpar.ncols; ++j){
                double dtemp;
                conv(dtemp, Amat[k][i][j]);
                cmsg[k][i*HEmatpar.dim + j].real(dtemp);
            }
        }
        Actxts[k] = scheme.encrypt(cmsg[k], HEmatpar.nslots, HEmatpar.pBits, HEmatpar.logQ);
    }
    NTL_EXEC_RANGE_END;
    
    delete[] cmsg;
}

void HEmatrix::HErmatmul_preprocessing(Ciphertext& res, Ciphertext*& Actxts, Ciphertext& Bctxt, ZZX*& Initpoly){
    
    //! 1. Generate the initial ciphertexts
    Ciphertext* Bctxts = new Ciphertext[HEmatpar.subdim];
    genInitBctxt(Bctxts[0], Bctxt, Initpoly);
    
    //! 2. Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.subdim - 1, first, last);
    for(int i = first; i < last; ++i){
        long i1 = (i + 1);
        Bctxts[i1] = scheme.leftRotate(Bctxts[0], HEmatpar.dim * (i1));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    NTL_EXEC_RANGE(HEmatpar.subdim, first, last);
    for(int i = first; i < last; ++i){
        scheme.modDownToAndEqual(Actxts[i], Bctxts[0].logq);
        scheme.multAndEqual(Actxts[i], Bctxts[i]);   //! log(Actxt[i+1]) = log Bctxt[0]
    }
    NTL_EXEC_RANGE_END;
    
    //! aggregate the results
    res = Actxts[0];
    for(int i = 1; i < HEmatpar.subdim; ++i){
        scheme.addAndEqual(res, Actxts[i]);
    }
    scheme.reScaleByAndEqual(res, res.logp);
    
    //! 4. shift and aggregate the results
    long index = (long) log2(HEmatpar.dim/HEmatpar.subdim);
    for(long i = 0; i < index; ++i){
        Ciphertext ctemp = scheme.leftRotate(res, HEmatpar.dim*HEmatpar.subdim * (1<<i));
        scheme.addAndEqual(res, ctemp);
    }
    
    delete[] Bctxts;
}
