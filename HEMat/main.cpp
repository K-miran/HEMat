//!
//!@ detals     Test the main.cpp
//!
//!@ author     Miran Kim
//!@ date       Dec.1, 2017

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>
#include <chrono>

#include <NTL/RR.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include "NTL/RR.h"
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#include <NTL/BasicThreadPool.h>

#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>

#include "../src/Context.h"
#include "../src/Scheme.h"
#include "../src/SecretKey.h"
#include "../src/TimeUtils.h"

#include "matrix.h"
#include "HEmatrix.h"
#include "TestHEmatrix.h"


using namespace NTL;
using namespace std;

int main() {
    
    SetNumThreads(4);
    
    long dim, dim1, dim2;        //! dimension of (square) matrix
    long nbatching;   //! number of multiple matrices in a single ciphertext
    long niter;       //! number of iterations
    
    while (true)
    {
        cout << "+------------------------------------------+" << endl;
        cout << "| Basic Examples                           |" << endl;
        cout << "+------------------------------------------+" << endl;
        cout << "| 1. testHEAAN                             |" << endl;
        cout << "| 2. testEnc                               |" << endl;
        cout << "| 3. Matrix addition                       |" << endl;
        cout << "| 4. Matrix transposition                  |" << endl;
        cout << "| 5. Matrix multiplication                 |" << endl;
        cout << "| 6. Rectangular matrix multiplication     |" << endl;
        cout << "+------------------------------------------+" << endl;
        cout << "+------------------------------------------+" << endl;
        cout << "| Pre-processing Examples                  |" << endl;
        cout << "+------------------------------------------+" << endl;
        cout << "| 7. Matrix multiplication                 |" << endl;
        cout << "| 8. Rectangular matrix multiplication     |" << endl;
        cout << "+------------------------------------------+" << endl;
        cout << "+------------------------------------------+" << endl;
        cout << "| Parallel Computation Examples            |" << endl;
        cout << "+------------------------------------------+" << endl;
        cout << "| 9. Parallel Matrix addition              |" << endl;
        cout << "| 10. Parallel Matrix transposition        |" << endl;
        cout << "| 11. Parallel Matrix multiplication       |" << endl;
        cout << "+------------------------------------------+" << endl;
        
        int selection = 0;
        bool invalid = true;
        do
        {
            cout << endl << "> Run example (1 ~ 11) or exit (0): ";
            if (!(cin >> selection))
            {
                invalid = false;
            }
            else if (selection < 0 || selection > 11)
            {
                invalid = false;
            }
            else
            {
                invalid = true;
            }
            if (!invalid)
            {
                cout << "  Invalid option: type 0 ~ 11" << endl;
                cin.clear();
            }
        } while (!invalid);
        
        switch (selection)
        {
            case 1:
                cout << endl << "> Enter log(N) and loq(q): ";
                long logn, logq;
                cin >> logn;
                cin >> logq;
                TestHEmatrix::testHEAAN(logn, logq);
                break;
                
            case 2:
                cout << endl << "> Enter a dim of matrix (as power-of-two): ";
                cin >> dim;
                TestHEmatrix::testEnc(dim);
                break;
                
            case 3:
                cout << endl << "> Enter a dim of matrix (as power-of-two): ";
                cin >> dim;
                TestHEmatrix::testAdd(dim);
                break;
                
            case 4:
                cout << endl << "> Enter a dim of matrix (as power-of-two): ";
                cin >> dim;
                TestHEmatrix::testTrans(dim);
                break;
                
            case 5:
                cout << endl << "> Enter a dim of matrix (as power-of-two): ";
                cin >> dim;
                TestHEmatrix::testMult(dim);
                break;
                
            case 6:
                cout << endl << "> Enter dim-1 & dim-2 (dim-1 >= dim-2, e.g., 32 8): ";
                cin >> dim1;
                cin >> dim2;
                cout << "(" << dim1 << "," << dim2 << ") * (" << dim2 << "," << dim2 << ")" << endl;
                TestHEmatrix::testRMult(dim1, dim2); //! (d2 * d1) * (d1 * d1)
                break;
                
            case 7:
                cout << endl << "> Enter a dim of matrix (as power-of-two): ";
                cin >> dim;
                TestHEmatrix::testMult_preprocessing(dim);
                break;
                
            case 8:
                cout << endl << "> Enter dim-1 & dim-2 (dim-1 >= dim-2, e.g., 32 8): ";
                cin >> dim1 >> dim2;
                cout << "(" << dim1 << "," << dim2 << ") * (" << dim2 << "," << dim2 << ")" << endl;
                TestHEmatrix::testRMult_preprocessing(dim1, dim2); //! (d1 * d2) * (d2 * d2)
                break;
                
            case 9:
                cout << endl << "> Enter dim & #(parallel matrices) & #(iterations): ";
                cin >> dim >> nbatching >> niter;
                TestHEmatrix::testSIMDAdd(dim, nbatching, niter);
                break;
                
            case 10:
                cout << endl << "> Enter dim & #(parallel matrices) & #(iterations): ";
                cin >> dim >> nbatching >> niter;
                TestHEmatrix::testSIMDTrans(dim, nbatching, niter);
                break;
                
            case 11:
                cout << endl << "> Enter dim & #(parallel matrices) & #(iterations): ";
                cin >> dim >> nbatching >> niter;
                TestHEmatrix::testSIMDMult(dim, nbatching, niter);  //! square matrix
                break;
                
            case 0:
                return 0;
        }
    }
    
    return 0;
}

