
#ifndef	_TESTING_UTILS_
#define _TESTING_UTILS_

#include <stdarg.h>

#define NRUNS	(10)


const char *usage =
"  --range start:stop:step\n"
"                   Adds test cases with range for sizes m,n. Can be repeated.\n"
"  -N m[:n]         Adds one test case with sizes m,n. Can be repeated.\n"
"                   If only m given then n=m.\n"
"  -m m             Sets m for all tests, overriding -N and --range.\n"
"  -n n             Sets n for all tests, overriding -N and --range.\n"
"\n"
"  -c  --[no]check  Whether to check results against CUBLAS, default is on.\n"
"  --dev x          GPU device to use, default 0.\n"
"\n"
"  --niter x        Number of iterations to repeat each test, default 1.\n"
"  -L -U            uplo   = Lower*, Upper.\n"
"  -[NTC][NTC]      transA = NoTrans*, Trans, or ConjTrans (first letter) and\n"
"                   transB = NoTrans*, Trans, or ConjTrans (second letter).\n"
"  -[TC]            transA = Trans or ConjTrans. Default is NoTrans. Doesn't change transB.\n"
"  -S[LR]           side   = Left*, Right.\n"
"  -D[NU]           diag   = NonUnit*, Unit.\n"
"                   * default values\n"
"\n"
"examples: \n"
"to test trmm with matrix A[512,512], B[2000,512] do\n"
"   test_dtrmm -N 2000:512 \n"
"to test trmm for range of sizes starting at 1024, stoping at 4096, steping at 1024, sizes will be for both A and B, with A upper traingular and transposed, do\n"
"   test_dtrmm --range 1024:4096:1024 -U -T";


#define USAGE printf("usage: -N m[:n] --range m-start:m-end:m-step -m INT -n INT -L|U -SL|SR -DN|DU -[NTC][NTC] -c --niter INT --dev devID\n\n"); \
printf("%s\n", usage);

#define USING printf("side %c, uplo %c, trans %c, diag %c, db %d\n", opts.side, opts.uplo, opts.transA, opts.diag, opts.db);

// random number gen --------------
float kblas_srand()
{
    return ( (float)rand() / (float)RAND_MAX );
}

double kblas_drand()
{
    return ( (double)rand() / (double)RAND_MAX );
}

/*cuFloatComplex kblas_crand()
{
    float x =  ( (float)rand() / (float)RAND_MAX );
    float y =  ( (float)rand() / (float)RAND_MAX );
    return make_cuFloatComplex(x, y);
}

cuDoubleComplex kblas_zrand()
{
    double x =  ( (double)rand() / (double)RAND_MAX );
    double y =  ( (double)rand() / (double)RAND_MAX );
    return make_cuDoubleComplex(x, y);
}

// --------------------------------
*/
float cget_magnitude(float* a){ return sqrt(a[0] * a[0] + a[1] * a[1]); }

double zget_magnitude(double* a) { return sqrt(a[0] * a[0] + a[1] * a[1]); }

float sget_max_error(float* ref, float *res, int n, int inc)
{
  int i;
  float max_err = -1.0;
  float err = -1.0;
  inc = abs(inc);
  for(i = 0; i < n; i++)
  {
    err = fabs(res[i * inc] - ref[i * inc]);
    if(ref[i * inc] != 0.0)err /= fabs(ref[i * inc]);
    if(err > max_err)max_err = err;
    //printf("[%2d]   %-.2f   %-.2f   %-.2e \n", i, ref[i], res[i], err);
  }
  return max_err;
}

double dget_max_error(double* ref, double *res, int n, int inc)
{
  int i;
  double max_err = -1.0;
  double err = -1.0;
  inc = abs(inc);
  for(i = 0; i < n; i++)
  {
    err = fabs(res[i * inc] - ref[i * inc]);
    if(ref[i * inc] != 0.0)err /= fabs(ref[i * inc]);
    if(err > max_err)max_err = err;
    //printf("[%2d]   %-.2f   %-.2f   %-.2e \n", i, ref[i], res[i], err);
  }
  return max_err;
}

float cget_max_error(float* ref, float *res, int n, int inc)
{
  int i;
  float max_err = 0.f;
  float err, m, diff[2];
  inc = abs(inc);
  for(i = 0; i < n*2; i+=2)
  {
    diff[0] = res[i * inc] - ref[i * inc];
    diff[1] = res[i * inc + 1] - ref[i * inc + 1];
    err = cget_magnitude( diff );
    diff[0] = ref[i * inc];
    diff[1] = ref[i * inc + 1];
    m = cget_magnitude(diff); 
    if( m > 0)
      err /= m;
    if(err > max_err)
      max_err = err;
  }
  return max_err;
}

double zget_max_error(double* ref, double* res, int n, int inc)
{
  int i;
  double max_err = 0.f;
  double err, m, diff[2];
  inc = abs(inc);
  for(i = 0; i < n*2; i+=2)
  {
    diff[0] = res[i * inc] - ref[i * inc];
    diff[1] = res[i * inc + 1] - ref[i * inc + 1];
    err = zget_magnitude( diff );
    diff[0] = ref[i * inc];
    diff[1] = ref[i * inc + 1];
    m = zget_magnitude(diff);
    if( m > 0)
      err /= m;
    if(err > max_err)
      max_err = err;
  }
  return max_err;
}

//======================================================================
float sget_max_error_matrix(float* ref, float *res, long m, long n, long lda)
{
  long i, j;
  float max_err = -1.0;
  float err = -1.0;
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++)
    {
      float ref_ = ref[j * lda + i];
      float res_ = res[j * lda + i];
      err = fabs(res_ - ref_);
      if(ref_ != 0.0)err /= fabs(ref_);
      if(err > max_err)max_err = err;
      //printf("\n[%2d]   %-.2f   %-.2f   %-.2e \n", i, ref_, res_, err);
    }
  }
  return max_err;
}

double dget_max_error_matrix(double* ref, double *res, long m, long n, long lda)
{
  long i, j;
  double max_err = -1.0;
  double err = -1.0;
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++)
    {
      double ref_ = ref[j * lda + i];
      double res_ = res[j * lda + i];
      err = fabs(res_ - ref_);
      if(ref_ != 0.0)err /= fabs(ref_);
      if(err > max_err)max_err = err;
      //printf("\n[%2d]   %-.2f   %-.2f   %-.2e \n", i, ref_, res_, err);
    }
  }
  return max_err;
}

float cget_max_error_matrix(float* ref, float* res, long m, long n, long lda)
{
  long i, j;
  float max_err = -1.0f;
  float err = -1.0, r, diff[2];
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++)
    {
      diff[0] = res[2*(j * lda + i)] - ref[2*(j * lda + i)];
      diff[1] = res[2*(j * lda + i)+1] - ref[2*(j * lda + i)+1];
      err = cget_magnitude( diff );
      diff[0] = ref[2*(j * lda + i)];
      diff[1] = ref[2*(j * lda + i)+1];
      r = cget_magnitude( diff );
      
      if(r > 0)
        err /= r;
      if(err > max_err)
        max_err = err;
    }
  }
  return max_err;
}

double zget_max_error_matrix(double* ref, double* res, long m, long n, long lda)
{
  long i, j;
  double max_err = -1.0;
  double err = -1.0, r, diff[2];
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++)
    {
      diff[0] = res[2*(j * lda + i)] - ref[2*(j * lda + i)];
      diff[1] = res[2*(j * lda + i)+1] - ref[2*(j * lda + i)+1];
      err = zget_magnitude( diff );
      diff[0] = ref[2*(j * lda + i)];
      diff[1] = ref[2*(j * lda + i)+1];
      r = zget_magnitude( diff );
      
      if(r > 0)
        err /= r;
      if(err > max_err)
        max_err = err;
    }
  }
  return max_err;
}

//=====================================================================================
void srand_matrix(long rows, long cols, float* A, long LDA)
{
    // fill in the entire matrix with random values
    long i;
    long size_a = cols * LDA;
    for(i = 0; i < size_a; i++) A[i] = ( (float)rand() ) / (float)RAND_MAX;
}

void drand_matrix(long rows, long cols, double* A, long LDA)
{
    // fill in the entire matrix with random values
    long i;
    long size_a = cols * LDA;
    for(i = 0; i < size_a; i++) A[i] = ( (double)rand() ) / (double)RAND_MAX;
}

void crand_matrix(long rows, long cols, float* A, long LDA)
{
  // fill in the entire matrix with random values
  long i;
  long size_a = cols * LDA * 2;
  for(i = 0; i < size_a; i++)
  {
    A[i]   = ( (float)rand() ) / (float)RAND_MAX;
  }
}

void zrand_matrix(long rows, long cols, double* A, long LDA)
{
  // fill in the entire matrix with random values
  long i;
  long size_a = cols * LDA * 2;
  for(i = 0; i < size_a; i++)
  {
    A[i] = ( (double)rand() ) / (double)RAND_MAX;
  }
}
//==============================================================================================

#if defined(__cplusplus)
extern "C"{
#endif
  
  void kblas_assert( int condition, const char* msg, ... )
  {
    if ( ! condition ) {
      printf( "Assert failed: " );
      va_list va;
      va_start( va, msg );
      vprintf( msg, va );
      printf( "\n" );
      exit(1);
    }
  }
  //==============================================================================================
#define MAX_NTEST 1000
  typedef struct kblas_opts
  {
    // matrix size
    int ntest;
    int msize[ MAX_NTEST ];
    int nsize[ MAX_NTEST ];
    int ksize[ MAX_NTEST ];
    int qsize[ MAX_NTEST ];
    
    // scalars
    int niter;
    int nruns;
    double      tolerance;
    int check;
    int verbose;
    int nb;
    int db;
    int custom;
    int warmup;
    int time;
    int bd;

    // lapack flags
    CBLAS_UPLO uplo;
    CBLAS_TRANSPOSE transA;
    CBLAS_TRANSPOSE transB;
    CBLAS_SIDE side;
    CBLAS_DIAG diag;
    CBLAS_ORDER order;
    
  } kblas_opts;
  
  //adapted from magma/testings.cpp
  int parse_opts( int argc, char** argv, kblas_opts *opts )
  {
    // negative flag indicating -m, -n, -k not given
    int m = -1;
    int n = -1;
    int k = -1;
    int n_start = 0, n_stop = 0, n_step = 0;
    int m_start = 0, m_stop = 0, m_step = 0;
    
    // fill in default values
    
    opts->niter    = 1;
    opts->nruns    = 4;
    opts->nb       = 64;  // ??
    opts->db       = 512;  // ??
    opts->tolerance = 30.;
    opts->check     = 0;
    opts->verbose   = 0;
    opts->custom   = 0;
    opts->warmup    = 0;
    opts->time    = 0;
    opts->bd    = -1;

    opts->uplo      = CblasLower;      
    opts->transA    = CblasNoTrans;    
    opts->transB    = CblasNoTrans;    
    opts->side      = CblasLeft;       
    opts->diag      = CblasNonUnit;    
    opts->order     = CblasColMajor;

    if(argc < 2){
      USAGE
      exit(0);
    }
    
    
    int info;
    int ntest = 0;
    for( int i = 1; i < argc; ++i ) {
      // ----- matrix size
      // each -N fills in next entry of msize, nsize, ksize and increments ntest
      if ( strcmp("-N", argv[i]) == 0 && i+1 < argc ) {
        kblas_assert( ntest < MAX_NTEST, "error: -N %s, max number of tests exceeded, ntest=%d.\n",
                     argv[i], ntest );
        i++;
        int m2, n2, k2, q2;
        info = sscanf( argv[i], "%d:%d:%d:%d", &m2, &n2, &k2, &q2 );
        if ( info == 4 && m2 >= 0 && n2 >= 0 && k2 >= 0 && q2 >= 0 ) {
          opts->msize[ ntest ] = m2;
          opts->nsize[ ntest ] = n2;
          opts->ksize[ ntest ] = k2;
          opts->qsize[ ntest ] = q2;
        }
        else if ( info == 3 && m2 >= 0 && n2 >= 0 && k2 >= 0 ) {
          opts->msize[ ntest ] = m2;
          opts->nsize[ ntest ] = n2;
          opts->ksize[ ntest ] = k2;
          opts->qsize[ ntest ] = k2;  // implicitly
        }
        else if ( info == 2 && m2 >= 0 && n2 >= 0 ) {
          opts->msize[ ntest ] = m2;
          opts->nsize[ ntest ] = n2;
          opts->ksize[ ntest ] = n2;  // implicitly
          opts->qsize[ ntest ] = n2;  // implicitly
        }
        else if ( info == 1 && m2 >= 0 ) {
          opts->msize[ ntest ] = m2;
          opts->nsize[ ntest ] = m2;  // implicitly
          opts->ksize[ ntest ] = m2;  // implicitly
          opts->qsize[ ntest ] = m2;  // implicitly
        }
        else {
          fprintf( stderr, "error: -N %s is invalid; ensure m >= 0, n >= 0, k >= 0, info=%d, m2=%d, n2=%d, k2=%d, q2=%d.\n",
                  argv[i],info,m2,n2,k2,q2 );
          exit(1);
        }
        ntest++;
      }
      // --range start:stop:step fills in msize[ntest:], nsize[ntest:], ksize[ntest:]
      // with given range and updates ntest
      else if ( strcmp("--range", argv[i]) == 0 && i+1 < argc ) {
        i++;
        int start, stop, step;
        info = sscanf( argv[i], "%d:%d:%d", &start, &stop, &step );
        if ( info == 3 && start >= 0 && stop >= 0 && step != 0 ) {
          for( int n = start; (step > 0 ? n <= stop : n >= stop); n += step ) {
            if ( ntest >= MAX_NTEST ) {
              printf( "warning: --range %s, max number of tests reached, ntest=%d.\n",
                     argv[i], ntest );
              break;
            }
            opts->msize[ ntest ] = n;
            opts->nsize[ ntest ] = n;
            opts->ksize[ ntest ] = n;
            opts->qsize[ ntest ] = n;
            ntest++;
          }
        }
        else {
          fprintf( stderr, "error: --range %s is invalid; ensure start >= 0, stop >= start, step > 0.\n",
                  argv[i] );
          exit(1);
        }
      }
      else if ( strcmp("--nrange", argv[i]) == 0 && i+1 < argc ) {
        i++;
        int start, stop, step;
        info = sscanf( argv[i], "%d:%d:%d", &start, &stop, &step );
        if ( info == 3 && start >= 0 && stop >= 0 && step != 0 ) {
          n_start = start;
          n_stop = stop;
          n_step = step;
        }
        else {
          fprintf( stderr, "error: --nrange %s is invalid; ensure start >= 0, stop >= 0, step != 0.\n",
                  argv[i] );
          exit(1);
        }
      }
      else if ( strcmp("--mrange", argv[i]) == 0 && i+1 < argc ) {
        i++;
        int start, stop, step;
        info = sscanf( argv[i], "%d:%d:%d", &start, &stop, &step );
        if ( info == 3 && start >= 0 && stop >= 0 && step != 0 ) {
          m_start = start;
          m_stop = stop;
          m_step = step;
        }
        else {
          fprintf( stderr, "error: --mrange %s is invalid; ensure start >= 0, stop >= 0, step != 0.\n",
                  argv[i] );
          exit(1);
        }
      }
      // save m, n, k if -m, -n, -k is given; applied after loop
      else if ( strcmp("-m", argv[i]) == 0 && i+1 < argc ) {
        m = atoi( argv[++i] );
        kblas_assert( m >= 0, "error: -m %s is invalid; ensure m >= 0.\n", argv[i] );
      }
      else if ( strcmp("-n", argv[i]) == 0 && i+1 < argc ) {
        n = atoi( argv[++i] );
        kblas_assert( n >= 0, "error: -n %s is invalid; ensure n >= 0.\n", argv[i] );
      }
      else if ( strcmp("-k", argv[i]) == 0 && i+1 < argc ) {
        k = atoi( argv[++i] );
        kblas_assert( k >= 0, "error: -k %s is invalid; ensure k >= 0.\n", argv[i] );
      }
      
      // ----- scalar arguments
      else if ( strcmp("--niter",   argv[i]) == 0 && i+1 < argc ) {
        opts->niter = atoi( argv[++i] );
        kblas_assert( opts->niter > 0,
                     "error: --niter %s is invalid; ensure niter > 0.\n", argv[i] );
      }
      else if ( strcmp("--nruns",   argv[i]) == 0 && i+1 < argc ) {
        opts->nruns = atoi( argv[++i] );
        kblas_assert( opts->nruns > 0,
                      "error: --nruns %s is invalid; ensure nruns > 0.\n", argv[i] );
      }
      else if ( strcmp("--tolerance", argv[i]) == 0 && i+1 < argc ) {
        opts->tolerance = atof( argv[++i] );
        kblas_assert( opts->tolerance >= 0 && opts->tolerance <= 1000,
                     "error: --tolerance %s is invalid; ensure tolerance in [0,1000].\n", argv[i] );
      }
      else if ( strcmp("--nb",      argv[i]) == 0 && i+1 < argc ) {
        opts->nb = atoi( argv[++i] );
        kblas_assert( opts->nb > 0,
                     "error: --nb %s is invalid; ensure nb > 0.\n", argv[i] );
      }
      else if ( strcmp("--db",      argv[i]) == 0 && i+1 < argc ) {
        opts->db = atoi( argv[++i] );
        kblas_assert( opts->db > 0,
                     "error: --db %s is invalid; ensure db > 0.\n", argv[i] );
      }
      else if ( strcmp("--bd",      argv[i]) == 0 && i+1 < argc ) {
        opts->bd = atoi( argv[++i] );
        kblas_assert( opts->bd >= 0,
                     "error: --bd %s is invalid; ensure db >= 0.\n", argv[i] );
      }
      // ----- boolean arguments
      // check results
      else if ( strcmp("-c",         argv[i]) == 0 ) { opts->check  = 1; }
      else if ( strcmp("-t",         argv[i]) == 0 ) { opts->time  = 1; }
      else if ( strcmp("-v",  argv[i]) == 0 ) { opts->verbose= 1;  }
      else if ( strcmp("-cu",         argv[i]) == 0 ) { opts->custom  = 1; }
      else if ( strcmp("-w",  argv[i]) == 0 ) { opts->warmup = 1;  }
      
      // ----- lapack flag arguments
      else if ( strcmp("-L",  argv[i]) == 0 ) { opts->uplo = CblasLower; }
      else if ( strcmp("-U",  argv[i]) == 0 ) { opts->uplo = CblasUpper; }
      else if ( strcmp("-NN", argv[i]) == 0 ) { opts->transA = CblasNoTrans;   opts->transB = CblasNoTrans;   }
      else if ( strcmp("-NT", argv[i]) == 0 ) { opts->transA = CblasNoTrans;   opts->transB = CblasTrans;     }
      else if ( strcmp("-TN", argv[i]) == 0 ) { opts->transA = CblasTrans;     opts->transB = CblasNoTrans;   }
      else if ( strcmp("-TT", argv[i]) == 0 ) { opts->transA = CblasTrans;     opts->transB = CblasTrans;     }
      else if ( strcmp("-T",  argv[i]) == 0 ) { opts->transA = CblasTrans;     }
      
      else if ( strcmp("-SL", argv[i]) == 0 ) { opts->side  = CblasLeft;  }
      else if ( strcmp("-SR", argv[i]) == 0 ) { opts->side  = CblasRight; }
      
      else if ( strcmp("-DN", argv[i]) == 0 ) { opts->diag  = CblasNonUnit; }
      else if ( strcmp("-DU", argv[i]) == 0 ) { opts->diag  = CblasUnit;    }
      
      // ----- usage
      else if ( strcmp("-h",     argv[i]) == 0 || strcmp("--help", argv[i]) == 0 ) {
        USAGE
        exit(0);
      }
      else {
        fprintf( stderr, "error: unrecognized option %s\n", argv[i] );
        exit(1);
      }
    }
    
    // fill in msize[:], nsize[:], ksize[:] if -m, -n, -k were given
    if(m_step != 0 && n_step != 0){
      for( int m = m_start; (m_step > 0 ? m <= m_stop : m >= m_stop); m += m_step ) {
        for( int n = n_start; (n_step > 0 ? n <= n_stop : n >= n_stop); n += n_step ) {
          if ( ntest >= MAX_NTEST ) {
            printf( "warning: --m/n_range, max number of tests reached, ntest=%d.\n",
                    ntest );
            break;
          }
          opts->msize[ ntest ] = m;
          opts->nsize[ ntest++ ] = n;
        }
      }
    }else
    if(m_step != 0 && n >= 0){
      for( int m = m_start; (m_step > 0 ? m <= m_stop : m >= m_stop); m += m_step ) {
        if ( ntest >= MAX_NTEST ) {
          printf( "warning: --m/n_range, max number of tests reached, ntest=%d.\n",
                  ntest );
          break;
        }
        opts->msize[ ntest ] = m;
        opts->nsize[ ntest++ ] = n;
        
      }
    }else
    if(n_step != 0 && m >= 0){
      for( int n = n_start; (n_step > 0 ? n <= n_stop : n >= n_stop); n += n_step ) {
        if ( ntest >= MAX_NTEST ) {
          printf( "warning: --m/n_range, max number of tests reached, ntest=%d.\n",
                  ntest );
          break;
        }
        opts->msize[ ntest ] = m;
        opts->nsize[ ntest++ ] = n;
      }
    }else{
      if ( m >= 0 ) {
        for( int j = 0; j < MAX_NTEST; ++j ) {
          opts->msize[j] = m;
        }
      }
      if ( n >= 0 ) {
        for( int j = 0; j < MAX_NTEST; ++j ) {
          opts->nsize[j] = n;
        }
      }
      if ( k >= 0 ) {
        for( int j = 0; j < MAX_NTEST; ++j ) {
          opts->ksize[j] = k;
        }
      }
      if ( m > 0 && n > 0) {
        ntest = 1;
      }
    }
    // if size not specified
    if ( ntest == 0 ) {
      fprintf(stderr, "please specify matrix size\n\n");
      //USAGE
      exit(0);
    }
    kblas_assert( ntest <= MAX_NTEST, "error: tests exceeded max allowed tests!\n" );
    opts->ntest = ntest;
    
    return 1;
  }// end parse_opts
  // Make a matrix symmetric/symmetric positive definite.
  // Increases diagonal by N.
  // Sets Aji = ( Aij ) for j < i, that is, copy lower triangle to upper triangle.
#if defined(__cplusplus)
}
#endif


//==============================================================================================
#if defined(__cplusplus)
template<class T>
void kblas_Rmake_hpd( int N, T* A, int lda )
{
  int i, j;
  for( i=0; i<N; ++i ) {
    A[i*(1+lda)] = A[i*(1+lda)] + N;
    for( j=0; j<i; ++j ) {
      A[j + i*lda] = A[i + j*lda];
    }
  }
}
template<class T>
void kblas_Cmake_hpd( int N, T* A, int lda )
{
  int i, j;
  for( i=0; i<N; ++i ) {
    A[2*i*(1+lda)] = A[2*i*(1+lda)] + N;
    for( j=0; j<i; ++j ) {
      A[2*(j + i*lda)] = A[2*(i + j*lda)];
      A[2*(j + i*lda)+1] = A[2*(i + j*lda)+1];
    }
  }
}
#endif//__cplusplus

/*double gettime(void)
{
  struct timeval tp;
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec;
}*/

//==============================================================================================
#define FMULS_TRMM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)+1))
#define FADDS_TRMM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)-1))
#define FMULS_TRMM(side_, m_, n_) ( ( (side_) == CblasLeft ) ? FMULS_TRMM_2((m_), (n_)) : FMULS_TRMM_2((n_), (m_)) )
#define FADDS_TRMM(side_, m_, n_) ( ( (side_) == CblasLeft ) ? FADDS_TRMM_2((m_), (n_)) : FADDS_TRMM_2((n_), (m_)) )


double FLOPS_TRMM(int cs, CBLAS_SIDE side, int m, int n){
  if(cs == 1)
    return FMULS_TRMM(side, (double)(m), (double)(n)) + FADDS_TRMM(side, (double)(m), (double)(n));
  else
    return 6. * FMULS_TRMM(side, (double)(m), (double)(n)) + 2. * FADDS_TRMM(side, (double)(m), (double)(n));
}

//==============================================================================================
#define FMULS_TRSM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)+1))
#define FADDS_TRSM_2(m_, n_) (0.5 * (n_) * (m_) * ((m_)-1))
#define FMULS_TRSM(side_, m_, n_) ( ( (side_) == CblasLeft ) ? FMULS_TRSM_2((m_), (n_)) : FMULS_TRSM_2((n_), (m_)) )
#define FADDS_TRSM(side_, m_, n_) ( ( (side_) == CblasLeft ) ? FADDS_TRSM_2((m_), (n_)) : FADDS_TRSM_2((n_), (m_)) )


double FLOPS_TRSM(float p, CBLAS_SIDE side, int m, int n){
  return FMULS_TRSM(side, (double)(m), (double)(n)) + FADDS_TRSM(side, (double)(m), (double)(n));
}
double FLOPS_TRSM(double p, CBLAS_SIDE side, int m, int n){
  return FMULS_TRSM(side, (double)(m), (double)(n)) + FADDS_TRSM(side, (double)(m), (double)(n));
}/*
double FLOPS_TRSM(void p, CBLAS_SIDE side, int m, int n){
  return 6. * FMULS_TRSM(side, (double)(m), (double)(n)) + 2. * FADDS_TRSM(side, (double)(m), (double)(n));
}*/


template<typename T>
void kblas_Raxpy (int n, T alpha, const T *x, int incx, T *y, int incy){
  int ix = 0, iy = 0;
  if(incx < 0) ix = 1 - n * incx;
  if(incy < 0) iy = 1 - n * incy;
  for(int i = 0; i < n; i++, ix+=incx, iy+=incy){
    y[iy] += alpha * x[ix];
  }
}
template<typename T>
void kblas_Caxpy (int n, T* alpha, const T *x, int incx, T *y, int incy){
  int ix = 0, iy = 0;
  if(incx < 0) ix = 1 - n * incx;
  if(incy < 0) iy = 1 - n * incy;
  for(int i = 0; i < n; i++, ix+=incx, iy+=incy){
    y[2*iy] += alpha[0] * x[2*ix] - alpha[1] * x[2*ix+1];
    y[2*iy+1] += alpha[1] * x[2*ix] + alpha[0] * x[2*ix+1];
  }
}
//==============================================================================================

template<typename T>
bool kblas_laisnan(T val1, T val2){
  return val1 != val2;
}

template<typename T>
bool kblas_isnan(T val){
  return kblas_laisnan(val,val);
}

float Xabs(float a){return fabs(a);}
double Xabs(double a){return fabs(a);}
//float Xabs(void* a){return Xget_magnitude(a);}
//double Xabs(cuDoubleComplex a){return zget_magnitude(a);}

template<typename T>
T kblas_Rlange(char type, int M, int N, T* arr, int lda){
  T value = 0.;
  T temp;
  for(int j = 0; j < N; j++){
    for(int i = 0; i < N; i++){
      temp = Xabs(arr[i + j * lda]);
      if( kblas_isnan(temp) || value < temp)
        value = temp;
    }
  }
  return value;
}

template<typename T>
T kblas_Clange(char type, int M, int N, T* arr, int lda){
  T value = 0.;
  T temp;
  for(int j = 0; j < N; j++){
    for(int i = 0; i < N; i++){
      temp = sqrt(arr[2*(i + j * lda)]*arr[2*(i + j * lda)] + arr[2*(i + j * lda)+1]*arr[2*(i + j * lda)+1]);
      if( kblas_isnan(temp) || value < temp)
        value = temp;
    }
  }
  return value;
}

double gettime(void)
{
  struct timeval tp;
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec;
}


int _kblas_error( int err, const char* func, const char* file, int line );

#ifdef check_error
#undef check_error
#endif



#define check_error( err ) \
{ \
  if(!_kblas_error( (err), __func__, __FILE__, __LINE__ )) void(); \
}
 //   return 0;\


#define TESTING_MALLOC_CPU( ptr, T, size)                       \
  if ( (ptr = (T*) malloc( (size)*sizeof( T ) ) ) == NULL) {    \
    fprintf( stderr, "!!!! malloc_cpu failed for: %s\n", #ptr ); \
    exit(-1);                                                   \
  }
#endif	//_TESTING_UTILS_