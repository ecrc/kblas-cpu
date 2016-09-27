/**
  - -* (C) Copyright 2013 King Abdullah University of Science and Technology
  Authors:
  Ali Charara (ali.charara@kaust.edu.sa)
  David Keyes (david.keyes@kaust.edu.sa)
  Hatem Ltaief (hatem.ltaief@kaust.edu.sa)
  
  Redistribution  and  use  in  source and binary forms, with or without
  modification,  are  permitted  provided  that the following conditions
  are met:
  
  * Redistributions  of  source  code  must  retain  the above copyright
  * notice,  this  list  of  conditions  and  the  following  disclaimer.
  * Redistributions  in  binary  form must reproduce the above copyright
  * notice,  this list of conditions and the following disclaimer in the
  * documentation  and/or other materials provided with the distribution.
  * Neither  the  name of the King Abdullah University of Science and
  * Technology nor the names of its contributors may be used to endorse
  * or promote products derived from this software without specific prior
  * written permission.
  * 
  THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
  LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/
#include <stdlib.h>
#include <stdio.h>
//#include "operators.h"
#include "Xtr_common.ch"


int kblas_back_door[KBLAS_BACKDOORS] = {-1};
//==============================================================================================
extern "C"{
const char* kblasGetErrorString( int error )
{
  switch( error ) {
    case KBLAS_UnknownError:
      return "unknown KBLAS error code";
    case KBLAS_NotSupported:
      return "Operation not supported";
    case KBLAS_NotImplemented:
      return "Operation not implemented yet";
    default:
      return "unknown KBLAS error code";
  }
}
int kblas_init(){
  for(int i = 0; i < KBLAS_BACKDOORS; i++)
    kblas_back_door[i] = -1;
  return 1;
}
}

// --------------------
//inline
int _kblas_error( int err, const char* func, const char* file, int line )
{
  if ( err != KBLAS_Success ) {
    fprintf( stderr, "KBLAS error: %s (%d) in %s at %s:%d\n",
             kblasGetErrorString( err ), err, func, file, line );
    return 0;
  }
  return 1;
}

//#define check_error( err ) \
//{if(!_kblas_error( (err), __func__, __FILE__, __LINE__ )) return 0;}

//==============================================================================================
bool REG_SIZE(int n){
  return ((n > 0) && !(n & (n - 1)));
}
int CLOSEST_REG_SIZE(int n){
  //TODO validate input
  if(n > 0){
    int res = 1;
    while (res < n){
      res = res << 1;
    }
    return res >> 1;
  }else{
    return 0;    
  }
}


int kblas_trmm_ib = 128;
int kblas_trsm_ib = 128;

/*/==============================================================================================
void cblas_Xgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc,
                 bool ISFLOAT){
  cblas_sgemm(Order, TransA, TransB,
              M, N, K,
              alpha, A, lda,
                     B, ldb,
              beta,  C, ldc);
}

void cblas_Xgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc,
                 bool ISFLOAT){
  cblas_dgemm(Order, TransA, TransB,
              M, N, K,
              alpha, A, lda,
                     B, ldb,
              beta,  C, ldc);
}
void cblas_Xgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc,
                 bool ISFLOAT){
  if(ISFLOAT)
    cblas_cgemm(Order, TransA, TransB,
                M, N, K,
                alpha, A, lda,
                       B, ldb,
                beta,  C, ldc);
  else
    cblas_zgemm(Order, TransA, TransB,
                M, N, K,
                alpha, A, lda,
                       B, ldb,
                beta,  C, ldc);
}
*/
