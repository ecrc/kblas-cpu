#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "kblas.h"
#include "testing_utils.h"

#define cblas_Xtrmm cblas_ctrmm
#define kblas_Xtrmm kblas_ctrmm
#define Xrand_matrix crand_matrix
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

  float alpha[2] = {0.29f,0.54f};
  
#include "test_Xtrmm.ch"
  
}


