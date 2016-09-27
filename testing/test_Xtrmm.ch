#ifndef _TEST_TRMM_
#define _TEST_TRMM_
  
  int nruns = opts.nruns;
  int M, N;
  int Am, An, Bm, Bn;
  int sizeA, sizeB;
  int lda, ldb;
  int ione     = 1;
  int ISEED[4] = {0,0,0,1};
  
  T *h_A, *h_B, *h_R;
  
  
  USING
 
  printf("    M     N     kblasTRMM_REC GF/s (ms)  cblasTRMM GF/s (ms)  SP_REC   Error\n");
  printf("====================================================================\n");
  for( int i = 0; i < opts.ntest; ++i ) {
    for( int iter = 0; iter < opts.niter; ++iter ) {
      double  gflops, perf,
              ref_avg_perf = 0.0, ref_sdev_perf = 0.0, ref_avg_time = 0.0,
              kblas_avg_perf = 0.0, kblas_sdev_perf = 0.0, kblas_avg_time = 0.0,
              ref_error = 0.0;
      M = opts.msize[i];
      N = opts.nsize[i];
      
      gflops = FLOPS_TRMM(COMPSIZE, opts.side, M, N ) / 1e9;
      
      printf("%5d %5d   ",
             (int) M, (int) N);
      fflush( stdout );
      
      if ( opts.side == CblasLeft ) {
        lda = Am = M;
        An = M;
      } else {
        lda = Am = N;
        An = N;
      }
      ldb = Bm = M;
      Bn = N;
            
      sizeA = lda*An*COMPSIZE;
      sizeB = ldb*Bn*COMPSIZE;
      
      TESTING_MALLOC_CPU( h_A, TT, sizeA);
      TESTING_MALLOC_CPU( h_B, TT, sizeB);
      TESTING_MALLOC_CPU( h_R, TT, sizeB);
      
      if(opts.check)
      {
        nruns = 1;
      }
      // Initialize matrix and vector
      //printf("Initializing on cpu .. \n");
      Xrand_matrix(Am, An, (TT*)h_A, lda);
      Xrand_matrix(Bm, Bn, (TT*)h_B, ldb);
      kblas_Xmake_hpd( Am, (TT*)h_A, lda );
      memcpy(h_R, h_B, sizeB * sizeof(TT));
      
      if(opts.warmup){
        cblas_Xtrmm(opts.order, opts.side, opts.uplo, opts.transA, opts.diag,
                    M, N,
                    alpha, h_A, lda,
                           h_R, ldb);
      }
      double time = 0;
      
      for(int r = 0; r < nruns; r++)
      {
        memcpy(h_R, h_B, sizeB * sizeof(TT));
        
        time = -gettime();
        kblas_Xtrmm(opts.order, opts.side, opts.uplo, opts.transA, opts.diag,
                    M, N,
                    alpha, h_A, lda,
                           h_R, ldb);
        time += gettime();
        perf = gflops / time;
        kblas_avg_perf += perf;
        kblas_sdev_perf += perf * perf;
        kblas_avg_time += time;
      }
      
      if(opts.check || opts.time){
        if(opts.time)
          memcpy(h_R, h_B, sizeB * sizeof(TT));
        
        for(int r = 0; r < nruns; r++)
        {
          if(opts.time)
            memcpy(h_B, h_R, sizeB * sizeof(TT));
          
          time = -gettime();
          cblas_Xtrmm(opts.order, opts.side, opts.uplo, opts.transA, opts.diag,
                      M, N,
                      alpha, h_A, lda,
                             h_B, ldb);
          time += gettime();
          perf = gflops / time;
          ref_avg_perf += perf;
          ref_sdev_perf += perf * perf;
          ref_avg_time += time;
        }

        
        if(opts.check)
          ref_error = Xget_max_error_matrix((TT*)h_B, (TT*)h_R, Bm, Bn, ldb);
      }


      free( h_A );
      free( h_B );
      free( h_R );
        ref_sdev_perf = sqrt((ref_sdev_perf - (ref_avg_perf * ref_avg_perf / nruns))/nruns);
      kblas_sdev_perf = sqrt((kblas_sdev_perf - (kblas_avg_perf * kblas_avg_perf / nruns))/nruns);

      printf(" %7.2f %7.2f %7.2f    %7.2f %7.2f %7.2f    %2.2f   %8.2e\n",
             gflops * nruns / kblas_avg_time, kblas_avg_time * 1000 / nruns, kblas_sdev_perf,
             gflops * nruns / ref_avg_time, ref_avg_time * 1000 / nruns, ref_sdev_perf,
             ref_avg_time / kblas_avg_time,
             ref_error );
    }
    if ( opts.niter > 1 ) {
      printf( "\n" );
    }
  }


#endif