 /**
 -- (C) Copyright 2013 King Abdullah University of Science and Technology
  Authors:
  Ahmad Abdelfattah (ahmad.ahmad@kaust.edu.sa)
  David Keyes (david.keyes@kaust.edu.sa)
  Hatem Ltaief (hatem.ltaief@kaust.edu.sa)

  Redistribution  and  use  in  source and binary forms, with or without
  modification,  are  permitted  provided  that the following conditions
  are met:

  * Redistributions  of  source  code  must  retain  the above copyright
    notice,  this  list  of  conditions  and  the  following  disclaimer.
  * Redistributions  in  binary  form must reproduce the above copyright
    notice,  this list of conditions and the following disclaimer in the
    documentation  and/or other materials provided with the distribution.
  * Neither  the  name of the King Abdullah University of Science and
    Technology nor the names of its contributors may be used to endorse 
    or promote products derived from this software without specific prior 
    written permission.

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
#ifndef _KBLAS_H_
#define _KBLAS_H_

#include "defs.h"
#ifdef USE_MKL
#include <mkl.h>
#else
#include "cblas.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
  

//CPU API
void kblas_strmm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const float alpha, const float *A, const int lda,
                                          float *B, const int ldb);
void kblas_dtrmm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const double alpha, const double *A, const int lda,
                                           double *B, const int ldb);
void kblas_ctrmm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                                          void *B, const int ldb);
void kblas_ztrmm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                                          void *B, const int ldb);


void kblas_strsm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);
void kblas_dtrsm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const double alpha, const double *A, const int lda,
                                            double *B, const int ldb);
void kblas_ctrsm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                                          void *B, const int ldb);
void kblas_ztrsm(const CBLAS_ORDER Order,
                 const CBLAS_SIDE Side, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                                          void *B, const int ldb);
#ifdef __cplusplus
}
#endif

#endif // _KBLAS_H_
