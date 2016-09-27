#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "kblas.h"
#include "testing_utils.h"

#define cblas_Xtrsm cblas_ctrsm
#define kblas_Xtrsm kblas_ctrsm
#define Xrand_matrix crand_matrix
#define kblas_Xlange kblas_Clange
#define kblas_Xaxpy kblas_Caxpy
#define cblas_Xtrmm cblas_ctrmm
#define kblas_Xmake_hpd kblas_Cmake_hpd
#define T void
#define TT float
#define COMPSIZE 2
#define Xget_max_error_matrix cget_max_error_matrix

//==============================================================================================
int main(int argc, char** argv)
{

  kblas_opts opts;
  if(!parse_opts( argc, argv, &opts )){
    USAGE;
    return -1;
  }
  
  TT mone[2] = {-1.f,0.f};
  TT alpha[2] = {0.29f,0.54f};
  TT im = ((TT*)alpha)[0]*((TT*)alpha)[0]+((TT*)alpha)[1]*((TT*)alpha)[1];
  TT invAlpha[2] = {((TT*)alpha)[0] / im , -((TT*)alpha)[1] / im};
  
#include "test_Xtrsm.ch"
  
}


