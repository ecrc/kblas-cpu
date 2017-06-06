

#include <set>
#include <stdio.h>
#include <string>
// #include "cblas.h"
#include "Xtr_common.ch"
//#include "cholinv_util.h"

//==============================================================================================

#define A(r_,c_) A_+(order == CblasColMajor ? (r_)+(c_)*lda : (r_)*lda+(c_))
#define B(r_,c_) B_+(order == CblasColMajor ? (r_)+(c_)*ldb : (r_)*ldb+(c_))
#define kblas_Xtrmm kblas_dtrmm
#define cblas_Xtrmm cblas_dtrmm
#define cblas_Xgemm cblas_dgemm

void kblas_dtrmm(const CBLAS_ORDER order,
                const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag,
                const int m, const int n,
                const double alpha, const double *A_, const int lda,
                double *B_, const int ldb){

  double one = 1.0f;

#include "Xtrmm.ch"

}