

#include <set>
#include <stdio.h>
#include <string>
// #include "cblas.h"
#include "Xtr_common.ch"
//#include "cholinv_util.h"

//==============================================================================================

#define A(r_,c_) A_+(order == CblasColMajor ? (r_)*2*sizeof(float)+(c_)*2*sizeof(float)*lda : (r_)*2*sizeof(float)*lda+(c_)*2*sizeof(float))
#define B(r_,c_) B_+(order == CblasColMajor ? (r_)*2*sizeof(float)+(c_)*2*sizeof(float)*ldb : (r_)*2*sizeof(float)*ldb+(c_)*2*sizeof(float))
#define kblas_Xtrmm kblas_ctrmm
#define cblas_Xtrmm cblas_ctrmm
#define cblas_Xgemm cblas_cgemm

void kblas_ctrmm(const CBLAS_ORDER order,
                 const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                 const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag,
                 const int m, const int n,
                 const void* alpha, const void *A_, const int lda,
                 void *B_, const int ldb){

  float one[2] = {1.0f,0.f};

#include "Xtrmm.ch"

}