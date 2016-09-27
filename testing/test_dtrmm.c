#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "kblas.h"
#include "testing_utils.h"

#define cblas_Xtrmm cblas_dtrmm
#define kblas_Xtrmm kblas_dtrmm
#define Xrand_matrix drand_matrix
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

  double alpha = 0.29;
  
#include "test_Xtrmm.ch"

  
}


