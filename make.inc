
# specify compilers
#CXX=g++
#CC=gcc
CXX=icpc
CC=icc

#customize the name of the generated static kblas library
LIB_KBLAS_NAME=kblas-cpu

#which library provides the CBLAS/BLAS implementation
BLAS_LIB = -mkl=parallel

#ifeq (${CC}, gcc)
#  COPTS+=-std=c99
#endif

