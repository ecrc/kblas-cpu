#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "kblas.h"
#include "testing_utils.h"

#define cblas_Xtrsm cblas_ztrsm
#define kblas_Xtrsm kblas_ztrsm
#define Xrand_matrix zrand_matrix
#define kblas_Xlange kblas_Clange
#define kblas_Xaxpy kblas_Caxpy
#define cblas_Xtrmm cblas_ztrmm
#define kblas_Xmake_hpd kblas_Cmake_hpd
#define T void
#define TT double
#define COMPSIZE 2
#define Xget_max_error_matrix zget_max_error_matrix

//==============================================================================================
int main(int argc, char** argv)
{

  kblas_opts opts;
  if(!parse_opts( argc, argv, &opts )){
    USAGE;
    return -1;
  }

  TT mone[2] = {-1.0,0.};
  TT alpha[2] = {0.29,0.54};
  TT im = ((TT*)alpha)[0]*((TT*)alpha)[0]+((TT*)alpha)[1]*((TT*)alpha)[1];
  TT invAlpha[2] = {((TT*)alpha)[0] / im , -((TT*)alpha)[1] / im};
  
#include "test_Xtrsm.ch"
  
}


