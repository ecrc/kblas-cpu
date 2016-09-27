#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "kblas.h"
#include "testing_utils.h"

#define cblas_Xtrsm cblas_dtrsm
#define kblas_Xtrsm kblas_dtrsm
#define Xrand_matrix drand_matrix
#define kblas_Xlange kblas_Rlange
#define kblas_Xaxpy kblas_Raxpy
#define cblas_Xtrmm cblas_dtrmm
#define kblas_Xmake_hpd kblas_Rmake_hpd
#define T double
#define TT double
#define COMPSIZE 1
#define Xget_max_error_matrix dget_max_error_matrix

//==============================================================================================
int main(int argc, char** argv)
{

  kblas_opts opts;
  if(!parse_opts( argc, argv, &opts )){
    USAGE;
    return -1;
  }
  
  TT mone = -1.0;
  TT alpha = 0.29;
  TT invAlpha = 1. / alpha;
  
#include "test_Xtrsm.ch"

  
}


