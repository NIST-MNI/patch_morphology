/////////////////////////////////////
//lsqnoneg, solve for Ax=b, where x >= 0
//OpenMP + MKL
//Yuancheng Luo, 2/2011
////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif //USE_OPENMP
#include <cblas.h>
#include <lapacke.h>
#include "nnls_internal.h"
#include "nnls.h"

//OMP/MKL  systems, no update/downdate
void nnls(REAL *A, REAL *b, REAL *x, int isTransposed, int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION){

  CBLAS_TRANSPOSE transaN = CblasNoTrans;
  CBLAS_TRANSPOSE transaT = CblasTrans;
  REAL alpha = 1.0f, negAlpha = -1.0f,  beta = 0.0f;
  int idx = 1, idy = 1;
  int OMPT=omp_get_num_threads();
  
  REAL *f = (REAL*)malloc(sizeof(REAL) * OMPT * m);
  REAL *g = (REAL*)malloc(sizeof(REAL) * OMPT * n);  //g[OMPT][n]
  int *z = (int*)malloc(sizeof(int) * OMPT * n);         //z[OMPT][n]
  int *zRemoved =  (int*)malloc(sizeof(int) * OMPT* n); 

  int *kIdx = (int*)malloc(sizeof(int) * OMPT * n);

  int maxmn = MAX(m, n);
  REAL *Apt = (REAL*)malloc(sizeof(REAL) * OMPT * n * m);
  REAL *xp = (REAL*)malloc(sizeof(REAL) * OMPT * maxmn );
  REAL *tx = (REAL*) malloc(sizeof(REAL) * OMPT * n);

  int NB = 1;
  /*int lwork =  (MAX(m, n) * (2 + NB));
  REAL *work = (REAL*) malloc(sizeof(REAL) *OMPT * lwork);*/
  int info;

  int s;
  int i, j, k;

  //Set number of threads
  //mkl_set_num_threads(MKLT); 
  //omp_set_num_threads(OMPT); 

  //double startT, endT;
  //startT = omp_get_wtime();
  
  //Set initial solution to 0
  memset(x, 0, sizeof(REAL) * nSys * n);

#pragma omp parallel private(s, i, j, k, info)
  {
    int tid = omp_get_thread_num();
    int nnlsIter, lsIter;
    
#pragma omp for
    for(s = 0; s < nSys; ++s){
      //Initials
      nnlsIter = lsIter = 0;
      int kCols = 0;
      //Set active-set to all variables
      for(i = 0; i < n; ++i)
	z[M2(tid, i, n)] = 1;
      
      do{
	//printf("nnlsIter %d, lsIter %d\n", nnlsIter, lsIter);
	////Compute negative gradient of f=1/2|Ax-b|^2 w.r.t. x and store into g=A'(b-Ax)
	
#ifndef USE_DOUBLE
	if(isTransposed){
	  //Compute f=-Ax, treat as A transpose for fortran call
	  cblas_sgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax
	  cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f, treat as A transpose for fortran call
	  cblas_sgemv(CblasColMajor,transaN, n, m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}else{
	  //Compute f=-Ax
	  cblas_sgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax
	  cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f
	  cblas_sgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}
#else
	if(isTransposed){
	  //Compute f=-Ax, treat as A transpose for fortran call
	  cblas_dgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax
	  cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f, treat as A transpose for fortran call
	  cblas_dgemv(CblasColMajor,transaN, n, m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}else{
	  //Compute f=-Ax
	  cblas_dgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax
	  cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f
	  cblas_dgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}
#endif
	//Check for termination condition by finding a variable in z, s.t. max(g)>0 and remove from z
	int vmax = -1;
	REAL maxg = TOL_TERMINATION;
	for(i = 0; i < n; ++i){
	  if(z[M2(tid, i, n)] && g[M2(tid, i, n)] > maxg){
	    maxg = g[M2(tid, i, n)];
	    vmax = i;
	  }
	}
	//All gradients are non-positive so terminate
	if(vmax == -1)
	  break;
	
	//Remove vmax from z
	z[M2(tid, vmax, n)] = 0;
	int colAdd = 1;
	
	//Solve uncontrained linear least squares subproblem
	do{
	  ++lsIter;
	  
	  if(colAdd){
	    //Build unconstrained system
	    for(i = 0; i < kCols; ++i){
	      int kdx = kIdx[M2(tid, i, n)];
	      if(isTransposed){
		for(j = 0; j < m; ++j)
		  Apt[M3(tid, i, j, n, m)] = A[M3(s, j, kdx, m, n)];
	      }else{
		memcpy(&Apt[M3(tid, i, 0, n, m)], &A[M3(s, kdx, 0, n, m)], sizeof(REAL) * m);
	      }
	    }
	    //Add new column
	    if(isTransposed){
	      for(j = 0; j < m; ++j)
		Apt[M3(tid, kCols, j, n, m)] = A[M3(s, j, vmax, m, n)];
	    }else{
	      memcpy(&Apt[M3(tid, kCols, 0, n, m)], &A[M3(s, vmax, 0, n, m)], sizeof(REAL) * m);
	    }
	    kIdx[M2(tid, kCols, n)] = vmax;
	    kCols++;
	    
	  }else{
	    //Shift sK (column list in set P) list left by one
	    for(i = kCols - 1; i >= 0; --i){
	      int kdx =  kIdx[M2(tid, i, n)];
	      //Deleted column
	      if(zRemoved[M2(tid, kdx, n)]){
		for(j = i + 1; j < kCols; ++j)
		  kIdx[M2(tid, j-1, n)] = kIdx[M2(tid, j, n)];
		--kCols;
	      }
	    }
	    //Build unconstrained system
	    for(i = 0; i < kCols; ++i){
	      int kdx = kIdx[M2(tid, i, n)];
	      if(isTransposed){
		for(j = 0; j < m; ++j)
		  Apt[M3(tid, i, j, n, m)] = A[M3(s, j, kdx, m, n)];
	      }else{
		memcpy(&Apt[M3(tid, i, 0, n, m)], &A[M3(s, kdx, 0, n, m)], sizeof(REAL) * m);
	      }
	    }
	  }
	  //Solve unconstrained system
	  int one = 1;
	  memcpy(&xp[M2(tid, 0, maxmn)], &b[M2(s, 0, m)], sizeof(REAL) * m);
#ifndef USE_DOUBLE
	  LAPACKE_sgels(LAPACK_COL_MAJOR,'N', m, kCols, one, &Apt[M3(tid, 0, 0, n, m)], m, &xp[M2(tid, 0, maxmn)], m );
#else
	  LAPACKE_dgels(LAPACK_COL_MAJOR,'N', m, kCols, one, &Apt[M3(tid, 0, 0, n, m)], m, &xp[M2(tid, 0, maxmn)], m );
#endif
	  //Load solution xp int tx
	  REAL minTx = HUGE_VAL;
	  memset(&tx[M2(tid, 0, n)], 0, sizeof(REAL) * n);
	  for(i = kCols - 1; i >= 0; --i){
	    int kdx = kIdx[M2(tid, i, n)];
	    tx[M2(tid, kdx, n)] = xp[M2(tid, i, maxmn)];
	    minTx = MIN(tx[M2(tid, kdx, n)], minTx);
	  }

	  if(minTx > 0){
	    //Accept solution, update x
	    //printf("Accept\n");
	    memcpy(&x[M2(s, 0, n)], &tx[M2(tid, 0, n)], sizeof(REAL) * n);
	    break;
	  }else{
	    //Reject solution, update subproblem
	    //Find index q in set P for negative z such that x/(x-z) is minimized
	    //printf("Reject\n");
	    REAL minAlpha = HUGE_VAL;
	    for(i = 0; i < kCols; ++i){
	      int kdx = kIdx[M2(tid, i, n)];
	      if(tx[M2(tid, kdx, n)] <= 0)
		minAlpha = MIN(minAlpha, x[M2(s, kdx, n)] / (x[M2(s, kdx, n)] - tx[M2(tid, kdx, n)]) );
	    }	  
	    for(i = 0; i < n; ++i)
	      tx[M2(tid, i, n)] -=  x[M2(s, i, n)];
	    
#ifndef USE_DOUBLE
	    cblas_saxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#else
	    cblas_daxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#endif	    
	    memset(&zRemoved[M2(tid, 0, n)], 0, sizeof(int) * n);
	    //Move from set P to set Z all elements whose corresponding X is 0 (guaranteed to have one element)
	    for(i = 0; i < kCols; ++i){
	      int kdx = kIdx[M2(tid, i, n)];
	      zRemoved[M2(tid, kdx, n)] = z[M2(tid, kdx, n)] = fabs(x[M2(s, kdx, n)]) <= TOL_TERMINATION;
	    }
	    
	    colAdd = 0;
	  }
	  
	}while(lsIter < maxLSIters);
	
	if(lsIter >= maxLSIters)
	  break;
	
	++nnlsIter;
      }while(nnlsIter < maxNNLSIters);
      
      //printf("nSys %d, nnlsIter %d, lsIter %d\n", s, nnlsIter, lsIter);  
    }
  }

  //endT = omp_get_wtime();
  //printf("Elapsed time %f\n", endT - startT);
  
  /*free(work); */
  free(Apt); 
  free(xp); 
  free(tx); 
  free(kIdx);
  free(zRemoved); 
  free(z); 
  free(f);
  free(g);
}//End of nnlsOMPSysMKL

//OMP systems, MKL updates
void nnls_updates(REAL *A, REAL *b, REAL *x, int isTransposed,  int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION){

  CBLAS_TRANSPOSE transaN = CblasNoTrans;
  CBLAS_TRANSPOSE transaT = CblasTrans;

  REAL alpha = 1.0f, negAlpha = -1.0f,  beta = 0.0f;
  int idx = 1, idy = 1;

  int OMPT=omp_get_num_threads();

  
  REAL *f = (REAL*)malloc(sizeof(REAL) * OMPT * m);
  REAL *g = (REAL*)malloc(sizeof(REAL) * OMPT * n);  //g[OMPT][n]
  int *z = (int*)malloc(sizeof(int) * OMPT * n);         //z[OMPT][n]
  int *zRemoved =  (int*)malloc(sizeof(int) * OMPT* n); 

  REAL *newCol = (REAL*)malloc(sizeof(REAL) * OMPT * m);
  REAL *oldCol = (REAL*)malloc(sizeof(REAL) * OMPT * m);

  int *kIdx = (int*)malloc(sizeof(int) * OMPT * n);

  REAL *Qt = (REAL*) malloc(sizeof(REAL) * OMPT * n * m);
  REAL *R =  (REAL*) malloc(sizeof(REAL) * OMPT * n * n );
  REAL *Qtb = (REAL*) malloc(sizeof(REAL) * OMPT * n);

  REAL *tx = (REAL*) malloc(sizeof(REAL) * OMPT * n);

  int s;
  int i, j, k;

  //Set number of threads
  //mkl_set_num_threads(MKLT);
  //omp_set_num_threads(OMPT);

  //double startT, endT;
  //startT = omp_get_wtime();
  
  //Set initial solution to 0
  memset(x, 0, sizeof(REAL) * nSys * n);

#pragma omp parallel private(s, i, j, k)
  {
    int tid = omp_get_thread_num();
    int nnlsIter, lsIter;
    
#pragma omp for
    for(s = 0; s < nSys; ++s){
      //Initials
      nnlsIter = lsIter = 0;
      int kCols = 0;
      //Set active-set to all variables
      for(i = 0; i < n; ++i)
	z[M2(tid, i, n)] = 1;
      
      do{
	//printf("nnlsIter %d, lsIter %d\n", nnlsIter, lsIter);
	////Compute negative gradient of f=1/2|Ax-b|^2 w.r.t. x and store into g=A'(b-Ax)
#ifndef USE_DOUBLE
	if(isTransposed){
	  //Compute f=-Ax, treat as A transpose for fortran call
	  cblas_sgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax
	  cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f, treat as A transpose for fortran call
	  cblas_sgemv(CblasColMajor,transaN, n, m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}else{
	  //Compute f=-Ax
	  cblas_sgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax 
	  cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f
	  cblas_sgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}
#else	
	if(isTransposed){
	  //Compute f=-Ax, treat as A transpose for fortran call
	  cblas_dgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax
	  cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f, treat as A transpose for fortran call
	  cblas_dgemv(CblasColMajor,transaN, n,m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}else{
	  //Compute f=-Ax
	  cblas_dgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
	  //Compute f=b-Ax 
	  cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
	  //Compute g=A'f
	  cblas_dgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
	}
#endif
	//Check for termination condition by finding a variable in z, s.t. max(g)>0 and remove from z
	int vmax = -1;
	REAL maxg = TOL_TERMINATION;
	for(i = 0; i < n; ++i){
	  if(z[M2(tid, i, n)] && g[M2(tid, i, n)] > maxg){
	    maxg = g[M2(tid, i, n)];
	    vmax = i;
	  }
	}
	//All gradients are non-positive so terminate
	if(vmax == -1)
	  break;
	
	//Remove vmax from z
	z[M2(tid, vmax, n)] = 0;
	int colAdd = 1;
	
	//Solve uncontrained linear least squares subproblem
	do{
	  ++lsIter;
	  
	  if(colAdd){
	    //printf("Column update\n");
	    ////Modified Gram-Schmidt update
	    //load old col
	    if(isTransposed){
	      for(i = 0 ; i < m; ++i) 
		newCol[M2(tid, i, m)] = A[M3(s, i, vmax, m, n)];
	    }else{
	      memcpy(&newCol[M2(tid, 0, m)], &A[M3(s, vmax, 0, n, m)], sizeof(REAL) * m);
	    }
	    memcpy(&oldCol[M2(tid, 0, m)], &newCol[M2(tid, 0, m)], sizeof(REAL) * m);
	 
	    //Orthogonalize newCol with all previous kCols columns in Q
	    for(i = 0; i < kCols ;++i){
	      int kdx = kIdx[M2(tid, i, n)];
#ifndef USE_DOUBLE
	      REAL sc = -cblas_sdot(m, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
	      cblas_saxpy(m, sc, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
	      //Update R
	      R[M3(tid, kdx, vmax, n, n)] = cblas_sdot(m, &Qt[M3(tid, kdx, 0, n, m)], idx, &oldCol[M2(tid, 0, m)], idy);
#else
	      REAL sc = -cblas_ddot(m, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
	      cblas_daxpy(m, sc, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
	      //Update R
	      R[M3(tid, kdx, vmax, n, n)] = cblas_ddot(m, &Qt[M3(tid, kdx, 0, n, m)], idx, &oldCol[M2(tid, 0, m)], idy);
#endif    
	    }
	   
	    //REAL norm2 = sqrtf( cblas_sdot(m, newCol, idx, newCol, idy));
#ifndef USE_DOUBLE
	    REAL norm2 = cblas_snrm2(m, &newCol[M2(tid, 0, m)], idx);
#else
	    REAL norm2 = cblas_dnrm2(m, &newCol[M2(tid, 0, m)], idx);
#endif
	    R[M3(tid, vmax, vmax, n, n)] = norm2;
	    
	    norm2 = 1.0f / norm2;
#ifndef USE_DOUBLE
	    cblas_sscal(m, norm2, &newCol[M2(tid, 0, m)], idx);
#else   
	    cblas_dscal(m, norm2, &newCol[M2(tid, 0, m)], idx);
#endif
	    memcpy(&Qt[M3(tid, vmax, 0, n, m)], &newCol[M2(tid, 0, m)], sizeof(REAL) * m);
#ifndef USE_DOUBLE
	    Qtb[M2(tid, vmax, n)] = cblas_sdot(m, &newCol[M2(tid, 0, m)], idx, &b[M2(s, 0, m)], idy);
#else
	    Qtb[M2(tid, vmax, n)] = cblas_ddot(m, &newCol[M2(tid, 0, m)], idx, &b[M2(s, 0, m)], idy);
#endif
	    //Append vmax to kIdx
	    kIdx[M2(tid, kCols, n)] = vmax;
	    kCols++;
	    
	  }else{
	    //printf("Column Downdate\n");
	    //Given's Rotations downdate
	    
	    //Rebuild kIdx list
	    for(i = kCols - 1; i >= 0; --i){
	      int kdx = kIdx[M2(tid, i, n)];
	      
	      //Deleted column
	      if(zRemoved[M2(tid, kdx, n)]){
		for(j = i + 1; j < kCols; ++j){
		  int jdx  = kIdx[M2(tid, j, n)];
		  
		  REAL gx = R[M3(tid, kdx, jdx, n, n)];
		  REAL ax = R[M3(tid, jdx, jdx, n, n)];
		  
		  //Compute c,s,-s,c coefficients for row/col jdx and row/col kdx
		  REAL givenc, givens;
		  //cblas_srotG(&ax, &gx, &givenc, &givens); //Not accurate, use below
		  if(ax == 0 && gx == 0){
		    givenc = 0;
		    givens = 0;
		  }else if(ax == 0){
		    givenc = 0;
		    givens = 1;
		  }else if(gx == 0){
		    givenc = 1;
		    givens = 0;
		  }else if(fabs(gx) > fabs(ax)){
		    REAL r = -ax / gx;
		    givens = 1 / sqrt(1 + r * r);
		    givenc = -givens * r;
		  }else{
		    REAL r = -gx / ax;
		    givenc = 1 / sqrt(1 + r * r);
		    givens = -givenc * r;
		  }
		  
#ifndef USE_DOUBLE		  
		  //Update R
		  cblas_srot(n, &R[M3(tid, jdx, 0, n, n)], idx, &R[M3(tid, kdx, 0, n, n)], idy, givenc, givens);
		  
		  //Update Qt
		  cblas_srot(m, &Qt[M3(tid, jdx, 0, n, m)], idx, &Qt[M3(tid, kdx, 0, n, m)], idy, givenc, givens);
		  
		  //Update Qtb
		  int one = 1;
		  cblas_srot(one, &Qtb[M2(tid, jdx, n)], idx, &Qtb[M2(tid, kdx, n)], idy, givenc, givens);
#else
		  //Update R
		  cblas_drot(n, &R[M3(tid, jdx, 0, n, n)], idx, &R[M3(tid, kdx, 0, n, n)], idy, givenc, givens);
		  
		  //Update Qt
		  cblas_drot(m, &Qt[M3(tid, jdx, 0, n, m)], idx, &Qt[M3(tid, kdx, 0, n, m)], idy, givenc, givens);
		  
		  //Update Qtb
		  int one = 1;
		  cblas_drot(one, &Qtb[M2(tid, jdx, n)], idx, &Qtb[M2(tid, kdx, n)], idy, givenc, givens);
#endif
		}
		
		//Shift sK (column list in set P) list left by one
		for(j = i + 1; j < kCols; ++j)
		  kIdx[M2(tid, j-1, n)] = kIdx[M2(tid, j, n)];
		
		--kCols;
	      }
	    }
	  }
	  
	  //Compute solution tx	
	  REAL minTx = HUGE_VAL;
	  memset(&tx[M2(tid, 0, n)], 0, sizeof(REAL) * n);
	  for(i = kCols - 1; i >= 0; --i){
	    int kdx = kIdx[M2(tid, i, n)];
	    //Compute tx via backsubstitution
	    REAL coeff = R[M3(tid, kdx, kdx, n, n)];
	    REAL pSum = 0;
	    for(j = 0; j < kCols; ++j){
	      int jdx = kIdx[M2(tid, j, n)];
	      pSum += (kdx != jdx) * R[M3(tid, kdx, jdx, n, n)] * tx[M2(tid, jdx, n)];
	    }
	    tx[M2(tid, kdx, n)] = (Qtb[M2(tid, kdx, n)] - pSum) / coeff;
	    minTx = MIN(tx[M2(tid, kdx, n)], minTx);
	  }
	  
	  if(minTx > 0){
	    //Accept solution, update x
	    //printf("Accept\n");
	    memcpy(&x[M2(s, 0, n)], &tx[M2(tid, 0, n)], sizeof(REAL) * n);
	    break;
	  }else{
	    //Reject solution, update subproblem
	    //Find index q in set P for negative z such that x/(x-z) is minimized
	    //printf("Reject\n");
	    REAL minAlpha = HUGE_VAL;
	    for(i = 0; i < kCols; ++i){
	      int kdx = kIdx[M2(tid, i, n)];
	      if(tx[M2(tid, kdx, n)] <= 0)
		minAlpha = MIN(minAlpha, x[M2(s, kdx, n)] / (x[M2(s, kdx, n)] - tx[M2(tid, kdx, n)]) );
	    }	  
	    for(i = 0; i < n; ++i)
	      tx[M2(tid, i, n)] -=  x[M2(s, i, n)];
	    
#ifndef USE_DOUBLE
	    cblas_saxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#else
	    cblas_daxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#endif	    
	    memset(&zRemoved[M2(tid, 0, n)], 0, sizeof(int) * n);
	    //Move from set P to set Z all elements whose corresponding X is 0 (guaranteed to have one element)
	    for(i = 0; i < kCols; ++i){
	      int kdx = kIdx[M2(tid, i, n)];
	      zRemoved[M2(tid, kdx, n)] = z[M2(tid, kdx, n)] = fabs(x[M2(s, kdx, n)]) <= TOL_TERMINATION;
	    }
	    
	    colAdd = 0;
	  }
	  
	}while(lsIter < maxLSIters);
	
	if(lsIter >= maxLSIters)
	  break;
	
	++nnlsIter;
      }while(nnlsIter < maxNNLSIters);
      
      //printf("nSys %d, nnlsIter %d, lsIter %d\n", s, nnlsIter, lsIter);  
    }
  }

  //endT = omp_get_wtime();
  //printf("Elapsed time %f\n", endT - startT);

  free(tx);
  free(Qtb);
  free(R);
  free(Qt);
  free(kIdx);
  free(newCol);
  free(oldCol);
  free(zRemoved);
  free(z);
  free(f);
  free(g);
}//End of nnlsOMPSysMKLUpdates


//single thread version
void nnls_updates_single(REAL *A, REAL *b, REAL *x, int isTransposed,  int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION){

  CBLAS_TRANSPOSE transaN = CblasNoTrans;
  CBLAS_TRANSPOSE transaT = CblasTrans;

  REAL alpha = 1.0f, negAlpha = -1.0f,  beta = 0.0f;
  int idx = 1, idy = 1;

  int OMPT=1;

  
  REAL *f = (REAL*)malloc(sizeof(REAL) * OMPT * m);
  REAL *g = (REAL*)malloc(sizeof(REAL) * OMPT * n);  //g[OMPT][n]
  int *z = (int*)malloc(sizeof(int) * OMPT * n);         //z[OMPT][n]
  int *zRemoved =  (int*)malloc(sizeof(int) * OMPT* n); 

  REAL *newCol = (REAL*)malloc(sizeof(REAL) * OMPT * m);
  REAL *oldCol = (REAL*)malloc(sizeof(REAL) * OMPT * m);

  int *kIdx = (int*)malloc(sizeof(int) * OMPT * n);

  REAL *Qt = (REAL*) malloc(sizeof(REAL) * OMPT * n * m);
  REAL *R =  (REAL*) malloc(sizeof(REAL) * OMPT * n * n );
  REAL *Qtb = (REAL*) malloc(sizeof(REAL) * OMPT * n);

  REAL *tx = (REAL*) malloc(sizeof(REAL) * OMPT * n);

  int s;
  int i, j, k;

  //Set number of threads
  //mkl_set_num_threads(MKLT);
  //omp_set_num_threads(OMPT);

  //double startT, endT;
  //startT = omp_get_wtime();
  
  //Set initial solution to 0
  memset(x, 0, sizeof(REAL) * nSys * n);

  {
    int tid = 0;
    int nnlsIter, lsIter;
    
    for(s = 0; s < nSys; ++s){
      //Initials
      nnlsIter = lsIter = 0;
      int kCols = 0;
      //Set active-set to all variables
      for(i = 0; i < n; ++i)
        z[M2(tid, i, n)] = 1;
      
      do{
        //printf("nnlsIter %d, lsIter %d\n", nnlsIter, lsIter);
        ////Compute negative gradient of f=1/2|Ax-b|^2 w.r.t. x and store into g=A'(b-Ax)
#ifndef USE_DOUBLE
        if(isTransposed){
          //Compute f=-Ax, treat as A transpose for fortran call
          cblas_sgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax
          cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f, treat as A transpose for fortran call
          cblas_sgemv(CblasColMajor,transaN, n, m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }else{
          //Compute f=-Ax
          cblas_sgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax 
          cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f
          cblas_sgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }
#else   
        if(isTransposed){
          //Compute f=-Ax, treat as A transpose for fortran call
          cblas_dgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax
          cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f, treat as A transpose for fortran call
          cblas_dgemv(CblasColMajor,transaN, n,m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }else{
          //Compute f=-Ax
          cblas_dgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax 
          cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f
          cblas_dgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }
#endif
        //Check for termination condition by finding a variable in z, s.t. max(g)>0 and remove from z
        int vmax = -1;
        REAL maxg = TOL_TERMINATION;
        for(i = 0; i < n; ++i){
          if(z[M2(tid, i, n)] && g[M2(tid, i, n)] > maxg){
            maxg = g[M2(tid, i, n)];
            vmax = i;
          }
        }
        //All gradients are non-positive so terminate
        if(vmax == -1)
          break;
        
        //Remove vmax from z
        z[M2(tid, vmax, n)] = 0;
        int colAdd = 1;
        
        //Solve uncontrained linear least squares subproblem
        do{
          ++lsIter;
          
          if(colAdd){
            //printf("Column update\n");
            ////Modified Gram-Schmidt update
            //load old col
            if(isTransposed){
              for(i = 0 ; i < m; ++i) 
                newCol[M2(tid, i, m)] = A[M3(s, i, vmax, m, n)];
            }else{
              memcpy(&newCol[M2(tid, 0, m)], &A[M3(s, vmax, 0, n, m)], sizeof(REAL) * m);
            }
            memcpy(&oldCol[M2(tid, 0, m)], &newCol[M2(tid, 0, m)], sizeof(REAL) * m);
         
            //Orthogonalize newCol with all previous kCols columns in Q
            for(i = 0; i < kCols ;++i){
              int kdx = kIdx[M2(tid, i, n)];
#ifndef USE_DOUBLE
              REAL sc = -cblas_sdot(m, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              cblas_saxpy(m, sc, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              //Update R
              R[M3(tid, kdx, vmax, n, n)] = cblas_sdot(m, &Qt[M3(tid, kdx, 0, n, m)], idx, &oldCol[M2(tid, 0, m)], idy);
#else
              REAL sc = -cblas_ddot(m, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              cblas_daxpy(m, sc, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              //Update R
              R[M3(tid, kdx, vmax, n, n)] = cblas_ddot(m, &Qt[M3(tid, kdx, 0, n, m)], idx, &oldCol[M2(tid, 0, m)], idy);
#endif    
            }
           
            //REAL norm2 = sqrtf( cblas_sdot(m, newCol, idx, newCol, idy));
#ifndef USE_DOUBLE
            REAL norm2 = cblas_snrm2(m, &newCol[M2(tid, 0, m)], idx);
#else
            REAL norm2 = cblas_dnrm2(m, &newCol[M2(tid, 0, m)], idx);
#endif
            R[M3(tid, vmax, vmax, n, n)] = norm2;
            
            norm2 = 1.0f / norm2;
#ifndef USE_DOUBLE
            cblas_sscal(m, norm2, &newCol[M2(tid, 0, m)], idx);
#else   
            cblas_dscal(m, norm2, &newCol[M2(tid, 0, m)], idx);
#endif
            memcpy(&Qt[M3(tid, vmax, 0, n, m)], &newCol[M2(tid, 0, m)], sizeof(REAL) * m);
#ifndef USE_DOUBLE
            Qtb[M2(tid, vmax, n)] = cblas_sdot(m, &newCol[M2(tid, 0, m)], idx, &b[M2(s, 0, m)], idy);
#else
            Qtb[M2(tid, vmax, n)] = cblas_ddot(m, &newCol[M2(tid, 0, m)], idx, &b[M2(s, 0, m)], idy);
#endif
            //Append vmax to kIdx
            kIdx[M2(tid, kCols, n)] = vmax;
            kCols++;
            
          }else{
            //printf("Column Downdate\n");
            //Given's Rotations downdate
            
            //Rebuild kIdx list
            for(i = kCols - 1; i >= 0; --i){
              int kdx = kIdx[M2(tid, i, n)];
              
              //Deleted column
              if(zRemoved[M2(tid, kdx, n)]){
                for(j = i + 1; j < kCols; ++j){
                  int jdx  = kIdx[M2(tid, j, n)];
                  
                  REAL gx = R[M3(tid, kdx, jdx, n, n)];
                  REAL ax = R[M3(tid, jdx, jdx, n, n)];
                  
                  //Compute c,s,-s,c coefficients for row/col jdx and row/col kdx
                  REAL givenc, givens;
                  //cblas_srotG(&ax, &gx, &givenc, &givens); //Not accurate, use below
                  if(ax == 0 && gx == 0){
                    givenc = 0;
                    givens = 0;
                  }else if(ax == 0){
                    givenc = 0;
                    givens = 1;
                  }else if(gx == 0){
                    givenc = 1;
                    givens = 0;
                  }else if(fabs(gx) > fabs(ax)){
                    REAL r = -ax / gx;
                    givens = 1 / sqrt(1 + r * r);
                    givenc = -givens * r;
                  }else{
                    REAL r = -gx / ax;
                    givenc = 1 / sqrt(1 + r * r);
                    givens = -givenc * r;
                  }
                  
#ifndef USE_DOUBLE                
                  //Update R
                  cblas_srot(n, &R[M3(tid, jdx, 0, n, n)], idx, &R[M3(tid, kdx, 0, n, n)], idy, givenc, givens);
                  
                  //Update Qt
                  cblas_srot(m, &Qt[M3(tid, jdx, 0, n, m)], idx, &Qt[M3(tid, kdx, 0, n, m)], idy, givenc, givens);
                  
                  //Update Qtb
                  int one = 1;
                  cblas_srot(one, &Qtb[M2(tid, jdx, n)], idx, &Qtb[M2(tid, kdx, n)], idy, givenc, givens);
#else
                  //Update R
                  cblas_drot(n, &R[M3(tid, jdx, 0, n, n)], idx, &R[M3(tid, kdx, 0, n, n)], idy, givenc, givens);
                  
                  //Update Qt
                  cblas_drot(m, &Qt[M3(tid, jdx, 0, n, m)], idx, &Qt[M3(tid, kdx, 0, n, m)], idy, givenc, givens);
                  
                  //Update Qtb
                  int one = 1;
                  cblas_drot(one, &Qtb[M2(tid, jdx, n)], idx, &Qtb[M2(tid, kdx, n)], idy, givenc, givens);
#endif
                }
                
                //Shift sK (column list in set P) list left by one
                for(j = i + 1; j < kCols; ++j)
                  kIdx[M2(tid, j-1, n)] = kIdx[M2(tid, j, n)];
                
                --kCols;
              }
            }
          }
          
          //Compute solution tx 
          REAL minTx = HUGE_VAL;
          memset(&tx[M2(tid, 0, n)], 0, sizeof(REAL) * n);
          for(i = kCols - 1; i >= 0; --i){
            int kdx = kIdx[M2(tid, i, n)];
            //Compute tx via backsubstitution
            REAL coeff = R[M3(tid, kdx, kdx, n, n)];
            REAL pSum = 0;
            for(j = 0; j < kCols; ++j){
              int jdx = kIdx[M2(tid, j, n)];
              pSum += (kdx != jdx) * R[M3(tid, kdx, jdx, n, n)] * tx[M2(tid, jdx, n)];
            }
            tx[M2(tid, kdx, n)] = (Qtb[M2(tid, kdx, n)] - pSum) / coeff;
            minTx = MIN(tx[M2(tid, kdx, n)], minTx);
          }
          
          if(minTx > 0){
            //Accept solution, update x
            //printf("Accept\n");
            memcpy(&x[M2(s, 0, n)], &tx[M2(tid, 0, n)], sizeof(REAL) * n);
            break;
          }else{
            //Reject solution, update subproblem
            //Find index q in set P for negative z such that x/(x-z) is minimized
            //printf("Reject\n");
            REAL minAlpha = HUGE_VAL;
            for(i = 0; i < kCols; ++i){
              int kdx = kIdx[M2(tid, i, n)];
              if(tx[M2(tid, kdx, n)] <= 0)
                minAlpha = MIN(minAlpha, x[M2(s, kdx, n)] / (x[M2(s, kdx, n)] - tx[M2(tid, kdx, n)]) );
            }     
            for(i = 0; i < n; ++i)
              tx[M2(tid, i, n)] -=  x[M2(s, i, n)];
            
#ifndef USE_DOUBLE
            cblas_saxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#else
            cblas_daxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#endif      
            memset(&zRemoved[M2(tid, 0, n)], 0, sizeof(int) * n);
            //Move from set P to set Z all elements whose corresponding X is 0 (guaranteed to have one element)
            for(i = 0; i < kCols; ++i){
              int kdx = kIdx[M2(tid, i, n)];
              zRemoved[M2(tid, kdx, n)] = z[M2(tid, kdx, n)] = fabs(x[M2(s, kdx, n)]) <= TOL_TERMINATION;
            }
            
            colAdd = 0;
          }
          
        }while(lsIter < maxLSIters);
        
        if(lsIter >= maxLSIters)
          break;
        
        ++nnlsIter;
      }while(nnlsIter < maxNNLSIters);
      
      //printf("nSys %d, nnlsIter %d, lsIter %d\n", s, nnlsIter, lsIter);  
    }
  }

  //endT = omp_get_wtime();
  //printf("Elapsed time %f\n", endT - startT);

  free(tx);
  free(Qtb);
  free(R);
  free(Qt);
  free(kIdx);
  free(newCol);
  free(oldCol);
  free(zRemoved);
  free(z);
  free(f);
  free(g);
}//End of nnlsOMPSysMKLUpdates


typedef struct
{
  REAL *f;
  REAL *g;  //g[OMPT][n]
  int *z;         //z[OMPT][n]
  int *zRemoved ; 

  REAL *newCol;
  REAL *oldCol;

  int *kIdx;

  REAL *Qt;
  REAL *R;
  REAL *Qtb;

  REAL *tx;
  int nSys;
  int m;
  int n;
} NNLS_DATA;

NNLS_HANDLE allocate_nnls( int nSys, int m, int n)
{
  const int OMPT=1;
  NNLS_DATA * data=(NNLS_DATA*)malloc(sizeof(NNLS_DATA));
  
  data->nSys = nSys;
  data->m = m;
  data->n = n;
  
  data->f = (REAL*)malloc(sizeof(REAL) * OMPT * m);
  data->g = (REAL*)malloc(sizeof(REAL) * OMPT * n);  //g[OMPT][n]
  data->z = (int*)malloc(sizeof(int) * OMPT * n);         //z[OMPT][n]
  data->zRemoved =  (int*)malloc(sizeof(int) * OMPT* n); 

  data->newCol = (REAL*)malloc(sizeof(REAL) * OMPT * m);
  data->oldCol = (REAL*)malloc(sizeof(REAL) * OMPT * m);

  data->kIdx = (int*)malloc(sizeof(int) * OMPT * n);

  data->Qt = (REAL*) malloc(sizeof(REAL) * OMPT * n * m);
  data->R =  (REAL*) malloc(sizeof(REAL) * OMPT * n * n );
  data->Qtb = (REAL*) malloc(sizeof(REAL) * OMPT * n);

  data->tx = (REAL*) malloc(sizeof(REAL) * OMPT * n);
  
  return data;
}

void free_nnls(NNLS_HANDLE _data)
{
  NNLS_DATA * data=(NNLS_DATA*)_data;
  
  free(data->tx);
  free(data->Qtb);
  free(data->R);
  free(data->Qt);
  free(data->kIdx);
  free(data->newCol);
  free(data->oldCol);
  free(data->zRemoved);
  free(data->z);
  free(data->f);
  free(data->g);
  
  free(_data);
}

void nnls2_updates_single(NNLS_HANDLE _data, REAL *A, REAL *b, REAL *x, int isTransposed,  int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION)
{
  const int OMPT=1;
  NNLS_DATA * data=(NNLS_DATA*)_data;

  CBLAS_TRANSPOSE transaN = CblasNoTrans;
  CBLAS_TRANSPOSE transaT = CblasTrans;

  REAL alpha = 1.0f, negAlpha = -1.0f,  beta = 0.0f;
  int idx = 1, idy = 1;
  
  REAL *f = data->f;
  REAL *g = data->g;  //g[OMPT][n]
  int *z = data->z;         //z[OMPT][n]
  int *zRemoved =  data->zRemoved; 

  REAL *newCol = data->newCol;
  REAL *oldCol = data->oldCol;

  int *kIdx = data->kIdx;

  REAL *Qt = data->Qt;
  REAL *R =  data->R;
  REAL *Qtb = data->Qtb;

  REAL *tx = data->tx;

  int s;
  int i, j, k;

  //Set number of threads
  //mkl_set_num_threads(MKLT);
  //omp_set_num_threads(OMPT);

  //double startT, endT;
  //startT = omp_get_wtime();
  
  //Set initial solution to 0
  memset(x, 0, sizeof(REAL) * nSys * n);

  {
    int tid = 0;
    int nnlsIter, lsIter;
    
    for(s = 0; s < nSys; ++s){
      //Initials
      nnlsIter = lsIter = 0;
      int kCols = 0;
      //Set active-set to all variables
      for(i = 0; i < n; ++i)
        z[M2(tid, i, n)] = 1;
      
      do{
        //printf("nnlsIter %d, lsIter %d\n", nnlsIter, lsIter);
        ////Compute negative gradient of f=1/2|Ax-b|^2 w.r.t. x and store into g=A'(b-Ax)
#ifndef USE_DOUBLE
        if(isTransposed){
          //Compute f=-Ax, treat as A transpose for fortran call
          cblas_sgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax
          cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f, treat as A transpose for fortran call
          cblas_sgemv(CblasColMajor,transaN, n, m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }else{
          //Compute f=-Ax
          cblas_sgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax 
          cblas_saxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f
          cblas_sgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }
#else   
        if(isTransposed){
          //Compute f=-Ax, treat as A transpose for fortran call
          cblas_dgemv(CblasColMajor,transaT, n,  m,  negAlpha, &A[M3(s, 0, 0, m, n)], n, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax
          cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f, treat as A transpose for fortran call
          cblas_dgemv(CblasColMajor,transaN, n,m, alpha, &A[M3(s, 0, 0, m, n)], n, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }else{
          //Compute f=-Ax
          cblas_dgemv(CblasColMajor,transaN, m,  n,  negAlpha, &A[M3(s, 0, 0, n, m)], m, &x[M2(s, 0, n)], idx, beta, &f[M2(tid, 0, m)], idy);
          //Compute f=b-Ax 
          cblas_daxpy(m, alpha, &b[M2(s, 0, m)], idx, &f[M2(tid, 0, m)], idy);
          //Compute g=A'f
          cblas_dgemv(CblasColMajor,transaT, m, n, alpha, &A[M3(s, 0, 0, n, m)], m, &f[M2(tid, 0, m)], idx, beta, &g[M2(tid, 0, n)], idy);
        }
#endif
        //Check for termination condition by finding a variable in z, s.t. max(g)>0 and remove from z
        int vmax = -1;
        REAL maxg = TOL_TERMINATION;
        for(i = 0; i < n; ++i){
          if(z[M2(tid, i, n)] && g[M2(tid, i, n)] > maxg){
            maxg = g[M2(tid, i, n)];
            vmax = i;
          }
        }
        //All gradients are non-positive so terminate
        if(vmax == -1)
          break;
        
        //Remove vmax from z
        z[M2(tid, vmax, n)] = 0;
        int colAdd = 1;
        
        //Solve uncontrained linear least squares subproblem
        do{
          ++lsIter;
          
          if(colAdd){
            //printf("Column update\n");
            ////Modified Gram-Schmidt update
            //load old col
            if(isTransposed){
              for(i = 0 ; i < m; ++i) 
                newCol[M2(tid, i, m)] = A[M3(s, i, vmax, m, n)];
            }else{
              memcpy(&newCol[M2(tid, 0, m)], &A[M3(s, vmax, 0, n, m)], sizeof(REAL) * m);
            }
            memcpy(&oldCol[M2(tid, 0, m)], &newCol[M2(tid, 0, m)], sizeof(REAL) * m);
         
            //Orthogonalize newCol with all previous kCols columns in Q
            for(i = 0; i < kCols ;++i){
              int kdx = kIdx[M2(tid, i, n)];
#ifndef USE_DOUBLE
              REAL sc = -cblas_sdot(m, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              cblas_saxpy(m, sc, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              //Update R
              R[M3(tid, kdx, vmax, n, n)] = cblas_sdot(m, &Qt[M3(tid, kdx, 0, n, m)], idx, &oldCol[M2(tid, 0, m)], idy);
#else
              REAL sc = -cblas_ddot(m, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              cblas_daxpy(m, sc, &Qt[M3(tid, kdx, 0, n,  m)], idx, &newCol[M2(tid, 0, m)], idy);
              //Update R
              R[M3(tid, kdx, vmax, n, n)] = cblas_ddot(m, &Qt[M3(tid, kdx, 0, n, m)], idx, &oldCol[M2(tid, 0, m)], idy);
#endif    
            }
           
            //REAL norm2 = sqrtf( cblas_sdot(m, newCol, idx, newCol, idy));
#ifndef USE_DOUBLE
            REAL norm2 = cblas_snrm2(m, &newCol[M2(tid, 0, m)], idx);
#else
            REAL norm2 = cblas_dnrm2(m, &newCol[M2(tid, 0, m)], idx);
#endif
            R[M3(tid, vmax, vmax, n, n)] = norm2;
            
            norm2 = 1.0f / norm2;
#ifndef USE_DOUBLE
            cblas_sscal(m, norm2, &newCol[M2(tid, 0, m)], idx);
#else   
            cblas_dscal(m, norm2, &newCol[M2(tid, 0, m)], idx);
#endif
            memcpy(&Qt[M3(tid, vmax, 0, n, m)], &newCol[M2(tid, 0, m)], sizeof(REAL) * m);
#ifndef USE_DOUBLE
            Qtb[M2(tid, vmax, n)] = cblas_sdot(m, &newCol[M2(tid, 0, m)], idx, &b[M2(s, 0, m)], idy);
#else
            Qtb[M2(tid, vmax, n)] = cblas_ddot(m, &newCol[M2(tid, 0, m)], idx, &b[M2(s, 0, m)], idy);
#endif
            //Append vmax to kIdx
            kIdx[M2(tid, kCols, n)] = vmax;
            kCols++;
            
          }else{
            //printf("Column Downdate\n");
            //Given's Rotations downdate
            
            //Rebuild kIdx list
            for(i = kCols - 1; i >= 0; --i){
              int kdx = kIdx[M2(tid, i, n)];
              
              //Deleted column
              if(zRemoved[M2(tid, kdx, n)]){
                for(j = i + 1; j < kCols; ++j){
                  int jdx  = kIdx[M2(tid, j, n)];
                  
                  REAL gx = R[M3(tid, kdx, jdx, n, n)];
                  REAL ax = R[M3(tid, jdx, jdx, n, n)];
                  
                  //Compute c,s,-s,c coefficients for row/col jdx and row/col kdx
                  REAL givenc, givens;
                  //cblas_srotG(&ax, &gx, &givenc, &givens); //Not accurate, use below
                  if(ax == 0 && gx == 0){
                    givenc = 0;
                    givens = 0;
                  }else if(ax == 0){
                    givenc = 0;
                    givens = 1;
                  }else if(gx == 0){
                    givenc = 1;
                    givens = 0;
                  }else if(fabs(gx) > fabs(ax)){
                    REAL r = -ax / gx;
                    givens = 1 / sqrt(1 + r * r);
                    givenc = -givens * r;
                  }else{
                    REAL r = -gx / ax;
                    givenc = 1 / sqrt(1 + r * r);
                    givens = -givenc * r;
                  }
                  
#ifndef USE_DOUBLE                
                  //Update R
                  cblas_srot(n, &R[M3(tid, jdx, 0, n, n)], idx, &R[M3(tid, kdx, 0, n, n)], idy, givenc, givens);
                  
                  //Update Qt
                  cblas_srot(m, &Qt[M3(tid, jdx, 0, n, m)], idx, &Qt[M3(tid, kdx, 0, n, m)], idy, givenc, givens);
                  
                  //Update Qtb
                  int one = 1;
                  cblas_srot(one, &Qtb[M2(tid, jdx, n)], idx, &Qtb[M2(tid, kdx, n)], idy, givenc, givens);
#else
                  //Update R
                  cblas_drot(n, &R[M3(tid, jdx, 0, n, n)], idx, &R[M3(tid, kdx, 0, n, n)], idy, givenc, givens);
                  
                  //Update Qt
                  cblas_drot(m, &Qt[M3(tid, jdx, 0, n, m)], idx, &Qt[M3(tid, kdx, 0, n, m)], idy, givenc, givens);
                  
                  //Update Qtb
                  int one = 1;
                  cblas_drot(one, &Qtb[M2(tid, jdx, n)], idx, &Qtb[M2(tid, kdx, n)], idy, givenc, givens);
#endif
                }
                
                //Shift sK (column list in set P) list left by one
                for(j = i + 1; j < kCols; ++j)
                  kIdx[M2(tid, j-1, n)] = kIdx[M2(tid, j, n)];
                
                --kCols;
              }
            }
          }
          
          //Compute solution tx 
          REAL minTx = HUGE_VAL;
          memset(&tx[M2(tid, 0, n)], 0, sizeof(REAL) * n);
          for(i = kCols - 1; i >= 0; --i){
            int kdx = kIdx[M2(tid, i, n)];
            //Compute tx via backsubstitution
            REAL coeff = R[M3(tid, kdx, kdx, n, n)];
            REAL pSum = 0;
            for(j = 0; j < kCols; ++j){
              int jdx = kIdx[M2(tid, j, n)];
              pSum += (kdx != jdx) * R[M3(tid, kdx, jdx, n, n)] * tx[M2(tid, jdx, n)];
            }
            tx[M2(tid, kdx, n)] = (Qtb[M2(tid, kdx, n)] - pSum) / coeff;
            minTx = MIN(tx[M2(tid, kdx, n)], minTx);
          }
          
          if(minTx > 0){
            //Accept solution, update x
            //printf("Accept\n");
            memcpy(&x[M2(s, 0, n)], &tx[M2(tid, 0, n)], sizeof(REAL) * n);
            break;
          }else{
            //Reject solution, update subproblem
            //Find index q in set P for negative z such that x/(x-z) is minimized
            //printf("Reject\n");
            REAL minAlpha = HUGE_VAL;
            for(i = 0; i < kCols; ++i){
              int kdx = kIdx[M2(tid, i, n)];
              if(tx[M2(tid, kdx, n)] <= 0)
                minAlpha = MIN(minAlpha, x[M2(s, kdx, n)] / (x[M2(s, kdx, n)] - tx[M2(tid, kdx, n)]) );
            }     
            for(i = 0; i < n; ++i)
              tx[M2(tid, i, n)] -=  x[M2(s, i, n)];
            
#ifndef USE_DOUBLE
            cblas_saxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#else
            cblas_daxpy(n, minAlpha, &tx[M2(tid, 0, n)], idx, &x[M2(s, 0, n)], idy);
#endif      
            memset(&zRemoved[M2(tid, 0, n)], 0, sizeof(int) * n);
            //Move from set P to set Z all elements whose corresponding X is 0 (guaranteed to have one element)
            for(i = 0; i < kCols; ++i){
              int kdx = kIdx[M2(tid, i, n)];
              zRemoved[M2(tid, kdx, n)] = z[M2(tid, kdx, n)] = fabs(x[M2(s, kdx, n)]) <= TOL_TERMINATION;
            }
            
            colAdd = 0;
          }
          
        }while(lsIter < maxLSIters);
        
        if(lsIter >= maxLSIters)
          break;
        
        ++nnlsIter;
      }while(nnlsIter < maxNNLSIters);
      
      //printf("nSys %d, nnlsIter %d, lsIter %d\n", s, nnlsIter, lsIter);  
    }
  }
  
}



/*
void writeOutput(REAL *A, REAL *b, REAL *x,  int nSys, int m, int n, char* fileName){
  //Write x in binary
  FILE *f_cpu = fopen(fileName, "wb");
  fwrite(x, 1,  sizeof(REAL) * nSys *  n, f_cpu);
  fclose(f_cpu);
  
  //Write x, Ax, and b in txt format
  //Compute error ||Ax-b||2

  f_cpu = fopen("cpu_result.txt","w");
  int s, i, j, k;

  for(s = 0; s < nSys; ++s){
    fprintf(f_cpu, "\nnSys %d Ax\n", s);
    REAL norm2 = 0;
    for(i = 0; i < m; ++i){
      REAL temp = 0;
      for(j = 0; j < n; ++j)
	temp += A[M3(s, i, j, m, n)] * x[M2(s, j, n) ]; 
      fprintf(f_cpu, "%f ", temp);
      temp -= b[M2(s, i, m)];
      norm2 += temp * temp;
    }
    norm2 = sqrt(norm2);
    printf("nSys %d, norm %f\n", s, norm2);
    fprintf(f_cpu, "\nnorm %f", norm2);
    
    fprintf(f_cpu, "\n\nb\n");
    for(i = 0; i < m; ++i) 
      fprintf(f_cpu, "%f ", b[M2(s, i, m) ]);
    fprintf(f_cpu, "\n");

    fprintf(f_cpu, "\n\nx\n");
    for(i = 0; i < n; ++i) 
      fprintf(f_cpu, "%f ", x[M2(s, i, n) ]);
    fprintf(f_cpu, "\n");
    
  }
  fclose(f_cpu);
} //end of writeOutput
*/
/*
int main(int argc, char *argv[]){
  int method = 0;
  int  M = 512; //Equations
  int  N = 512; //Unknowns
  int  NSYS = 256; //1-256
  int isTransposed = 0;
  REAL TOL_TERMINATION = 1e-6;
  int MKLThreads = 1, OMPThreads = 1;
  REAL *A = NULL, *b = NULL, *x = NULL;
  if(argc < 9){
    printf("Expect: method, nSys, isTransposed, mEquations, nUnknowns, tol, MKLThreads, OMPThreads\n");
  }else{
    method =  atoi(argv[1]);
    NSYS = atoi(argv[2]);
    isTransposed = atoi(argv[3]);
    M = atoi(argv[4]);
    N = atoi(argv[5]);
    TOL_TERMINATION = atof(argv[6]);
    MKLThreads = atoi(argv[7]);
    OMPThreads = atoi(argv[8]);

    srand(2011);
    loadData(&A, &b, &x, NSYS, M, N, "sysA.bin", "sysB.bin");
    
    if(method == 0)
      nnlsOMPSysMKL(A, b, x, 0, MAX_ITER_NNLS(M, N), MAX_ITER_LS(M, N), NSYS, M, N, MKLThreads, OMPThreads, TOL_TERMINATION);
    else
      nnlsOMPSysMKLUpdates(A, b, x, isTransposed,  MAX_ITER_NNLS(M, N), MAX_ITER_LS(M, N), NSYS, M, N, MKLThreads, OMPThreads, TOL_TERMINATION);  
    
    writeOutput(A, b, x, NSYS, M, N, "sysX.bin");
    
    if(A) free(A);
    if(b) free(b);
    if(x) free(x);
  }
  return 0;
}//End of main
*/

#if 0
//For mex calls
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // printf("in mex function\n");
  if(nrhs != 10 || nlhs != 1)
    mexErrMsgTxt("Expect: x = NNLS(method, A, b, nSys, isTransposed, mEquations, nUnknowns, tol, MKLThreads, OMPThreads)\n");
#ifndef USE_DOUBLE
  else if(!mxIsSingle(prhs[1]) || !mxIsSingle(prhs[2]) )
    mexErrMsgTxt("A, b must be single precision\n");
#else
  else if(!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) )
    mexErrMsgTxt("A, b must be double precision\n");
#endif
  else{
    int method = mxGetScalar(prhs[0]);       
    REAL *A = (REAL*)mxGetData(prhs[1]);
    REAL *b = (REAL*)mxGetData(prhs[2]);
    int NSYS =  mxGetScalar(prhs[3]);
    int isTransposed =  mxGetScalar(prhs[4]); //!0 => A transposed in memory, 0 => A not transposed in memory
    int M = mxGetScalar(prhs[5]);           //Number of equations
    int N = mxGetScalar(prhs[6]);           //Number of unknowns
    REAL TOL_TERMINATION = mxGetScalar(prhs[7]);           //Tolerance for termination (0)
    int MKLT = mxGetScalar(prhs[8]);        //Number of MKL threads
    int OMPT = mxGetScalar(prhs[9]);        //Number of OMP threads

#ifndef USE_DOUBLE
    plhs[0] = mxCreateNumericMatrix(N, NSYS, mxSINGLE_CLASS, mxREAL);
#else
    plhs[0] = mxCreateNumericMatrix(N, NSYS, mxDOUBLE_CLASS, mxREAL);
#endif
    REAL *x = (REAL*)mxGetPr(plhs[0]);   
   
    //Naive (no update/downdates)
    if(method == 0)
      nnlsOMPSysMKL(A, b, x, isTransposed, MAX_ITER_NNLS(M, N), MAX_ITER_LS(M, N), NSYS, M, N, MKLT, OMPT, TOL_TERMINATION);           //Naive implementation, no update/downdate,  MKL, OMP per system
    else 
      nnlsOMPSysMKLUpdates(A, b, x, isTransposed,  MAX_ITER_NNLS(M, N), MAX_ITER_LS(M, N), NSYS, M, N, MKLT, OMPT, TOL_TERMINATION);   //OMP per system, MKL update/downdate 
  }
  
}//end of mexFunction
#endif

// kate: space-indent on; hl C;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8 
