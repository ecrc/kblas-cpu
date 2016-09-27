

#include <set>
#include <stdio.h>
#include <string>
#include "cblas.h"
#include "Xtr_common.ch"
//#include "cholinv_util.h"


//==============================================================================================

#define A(r_,c_) A_+(order == CblasColMajor ? (r_)+(c_)*lda : (r_)*lda+(c_) )
#define B(r_,c_) B_+(order == CblasColMajor ? (r_)+(c_)*ldb : (r_)*ldb+(c_) )
#define kblas_Xtrsm kblas_strsm
#define cblas_Xtrsm cblas_strsm
#define cblas_Xgemm cblas_sgemm

void kblas_strsm(const enum CBLAS_ORDER order,
                 const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const enum CBLAS_DIAG diag,
                 const int m, const int n,
                 const float alpha, const float *A_, const int lda,
                 float *B_, const int ldb){


  float one = 1.0f;
  float mone = -1.0f;
  float mInvAlpha = mone / alpha;

#include "Xtrsm.ch"
  
}