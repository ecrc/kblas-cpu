

#include <set>
#include <stdio.h>
#include <string>
#include "cblas.h"
#include "Xtr_common.ch"
//#include "cholinv_util.h"

//==============================================================================================

#define A(r_,c_) A_+(order == CblasColMajor ? (r_)*2*sizeof(double)+(c_)*2*sizeof(double)*lda : (r_)*2*sizeof(double)*lda+(c_)*2*sizeof(double) )
#define B(r_,c_) B_+(order == CblasColMajor ? (r_)*2*sizeof(double)+(c_)*2*sizeof(double)*ldb : (r_)*2*sizeof(double)*ldb+(c_)*2*sizeof(double) )
#define kblas_Xtrsm kblas_ztrsm
#define cblas_Xtrsm cblas_ztrsm
#define cblas_Xgemm cblas_zgemm

void kblas_ztrsm(const enum CBLAS_ORDER order,
                 const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const enum CBLAS_DIAG diag,
                 const int m, const int n,
                 const void* alpha, const void *A_, const int lda,
                                          void *B_, const int ldb){
   
  double one[2] = {1.0f,0.f};
  double mone[2] = {-1.0f,0.f};
  double im = ((double*)alpha)[0]*((double*)alpha)[0]+((double*)alpha)[1]*((double*)alpha)[1];
  double mInvAlpha[2] = {-((double*)alpha)[0] / im , ((double*)alpha)[1] / im};

#include "Xtrsm.ch"
  
}