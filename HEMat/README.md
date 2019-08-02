# HEMat

HEMat is a software package for performing a secure outsourced matrix computation using homomorphic encryption (https://dl.acm.org/citation.cfm?id=3243837). 

## Setting up HEMat library 

### Dependencies

Our library requires a c++ compiler and the following libraries:

a. the GMP (GNU Multi-Precision library, >= 6.1.2), which is available at https://gmplib.org,

b. the NTL library (>=11.0.0), which is available at http://www.shoup.net/ntl/,  (with pThread)

c. the approximate HE library, which is an implementation of the paper "Homomorphic Encryption for Arithmetic of Approximate Numbers" (https://eprint.iacr.org/2016/421.pdf). We refered to the underlying HE library in the "src" folder. You can build the libarary “libheaan.a" by typing "$make all" in the "/src" directory.


### Installing HEMat library

HEMat is easy to configure and build in Linux and macOS. 

First, clone this Github repository into a directory of your choice:

```sh
git clone https://github.com/K-miran/HEMat
```

You can then install our library by the following work with command line tools in the "/HEMat" directory. :

```sh
cd HEMat
make new
```

This will build our library "libHEMat.a".  

## Running a source code

The program will run with several secure matrix algorithms (e.g. matrix addition, multiplication, transposition, and rectangular matrix multiplication).  
For example, you run a test program in the the directory by doing as follows: 

```sh
make test
./foo 
```

