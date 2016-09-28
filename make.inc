
# specify compilers
#CXX=g++
#CC=gcc
CXX=icpc
CC=icc

#customize the name of the generated kblas library
LIB_KBLAS_NAME=libkblas.a

CBLAS_LIB = -mkl=parallel

#ifeq (${CC}, gcc)
#  COPTS+=-std=c99
#endif

