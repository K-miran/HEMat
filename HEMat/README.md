# HEMat

This is an implementation of secure outsource matrix computation based on approximate homomorphic encryption. 
This also contains an implementation of an evaluation of encrypted neural networks model on encrypted data.

Our code requires a c++ compiler and the following libraries:

a. the GMP (GNU Multi-Precision library), which is available at https://gmplib.org,

b. the NTL library, which is available at http://www.shoup.net/ntl/,  (with pThread)

c. the approximate HE library, which is an implementation of the paper "Homomorphic Encryption for Arithmetic of Approximate Numbers" (https://eprint.iacr.org/2016/421.pdf). We refered to the underlying HE library in the "src" folder. You can build the libarary “libheaan.a" by typing "$make all" in the "/src" directory.

## How to build and test the HEMat library?

You can build our libarary "libHEMat.a" by typing "$make new" in the "/HEMat" directory. 
You can run a test program in the the directory by typing "$make test".

## How to build and test the EDM?

You can build our libarary "libEDM.a" by typing "$make new" in the "/EDM" directory. 
You can run a test program in the the directory by typing "$make test".
