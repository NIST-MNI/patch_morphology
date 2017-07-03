#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef USE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_wtime() ((double)(clock())/CLOCKS_PER_SEC)
#define omp_get_max_threads() 1
#endif //USE_OPENMP

#include "nnls.h"

int main(int argc,char **argv)
{
  // generate a plausable system
  // solve it using two methods
  const int n=128;
  const int m=4096;
  int i,j;
  double *A=malloc(sizeof(double)*m*n);
  double *b=malloc(sizeof(double)*m);
  double *x1=malloc(sizeof(double)*n);
  double *x2=malloc(sizeof(double)*n);
  double *x3=malloc(sizeof(double)*n);
  double *x4=malloc(sizeof(double)*n);
  NNLS_HANDLE work;
  
  double startT, endT;

  double noise_level=1e-4;
  int    errors_detected=0;
  
  int negative_vars=n/4;

  printf("Using OpenMP, max number of threads=%d\n",omp_get_max_threads());
  
  srand(1024);
  
  for(i=0;i<m;i++)
  {
    double noise=noise_level* (rand()-RAND_MAX/2)/(double)RAND_MAX;
    for(j=0;j<n;j++)
    {
      //column-major mode
      A[i+j*m]=( (i%n)==j?1.0:0);
    }

    b[i]=((i%n)<negative_vars?-(i%n+1):i%n+1)+noise;
  }
  
  
  printf("Solution1:\n");
  startT = omp_get_wtime();
  nnls(A, b, x1, 0, (m+n)*4, (m+n)*4, 1,  m, n, 1e-20);
  endT = omp_get_wtime();
  for(j=0;j<n;j++)
  {
    printf("%g,",x1[j]);
  }
  printf("\n");
  printf("Elapsed time %f\n", endT - startT);
  
  printf("Solution2:\n");
  startT = omp_get_wtime();
  nnls_updates(A, b, x2, 0, (m+n)*4, (m+n)*4, 1,  m, n, 1e-20);
  endT = omp_get_wtime();
  for(j=0;j<n;j++)
  {
    printf("%g,",x2[j]);
  }
  printf("\n");
  printf("Elapsed time %f\n", endT - startT);

  printf("Solution3 (single thread):\n");
  startT = omp_get_wtime();
  nnls_updates_single(A, b, x3, 0, (m+n)*4, (m+n)*4, 1,  m, n, 1e-20);
  endT = omp_get_wtime();
  for(j=0;j<n;j++)
  {
    printf("%g,",x3[j]);
  }
  printf("\n");
  printf("Elapsed time %f\n", endT - startT);
  
  work=allocate_nnls(1, m, n);
  printf("Solution4 (single thread, with pre-allocated work space):\n");
  startT = omp_get_wtime();
  nnls2_updates_single(work,A, b, x4, 0, (m+n)*4, (m+n)*4, 1,  m, n, 1e-20);
  endT = omp_get_wtime();
  for(j=0;j<n;j++)
  {
    printf("%g,",x4[j]);
  }
  printf("\n");
  printf("Elapsed time %f\n", endT - startT);
  free_nnls(work);
  
  //check if first negative_vars are zero, and the rest is approximately equal to i
  for(j=0;j<n;j++)
  {
    if(j<negative_vars)
    {
      if(x1[j]!=0.0)
      {
        fprintf(stderr,"Error: solution 1 element %d expected 0, got %g\n",j,x1[j]);
        errors_detected++;
      }
      if(x2[j]!=0.0)
      {
        fprintf(stderr,"Error: solution 2 element %d expected 0, got %g\n",j,x1[j]);
        errors_detected++;
      }
      if(x3[j]!=0.0)
      {
        fprintf(stderr,"Error: solution 3 element %d expected 0, got %g\n",j,x1[j]);
        errors_detected++;
      }
      if(x4[j]!=0.0)
      {
        fprintf(stderr,"Error: solution 4 element %d expected 0, got %g\n",j,x1[j]);
        errors_detected++;
      }
    } else {
      if(fabs(x1[j]-j-1)>noise_level)
      {
        fprintf(stderr,"Error: solution 1 element %d expected %g -+ %g, got %g\n",j,(double)j+1,noise_level,x1[j]);
        errors_detected++;
      }
      if(fabs(x2[j]-j-1)>noise_level)
      {
        fprintf(stderr,"Error: solution 2 element %d expected %g -+ %g, got %g\n",j,(double)j+1,noise_level,x2[j]);
        errors_detected++;
      }
      if(fabs(x3[j]-j-1)>noise_level)
      {
        fprintf(stderr,"Error: solution 3 element %d expected %g -+ %g, got %g\n",j,(double)j+1,noise_level,x3[j]);
        errors_detected++;
      }
      if(fabs(x4[j]-j-1)>noise_level)
      {
        fprintf(stderr,"Error: solution 4 element %d expected %g -+ %g, got %g\n",j,(double)j+1,noise_level,x4[j]);
        errors_detected++;
      }
    }
    return errors_detected;
  }
  
  return 0;
}
