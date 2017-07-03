#ifndef __NNLS_H__
#define __NNLS_H__



#ifdef __cplusplus
extern "C" {
#endif
  
#include "nnls_internal.h"


// all data is stored in Column major mode (fortran)

// A  - matrix (fortran) m rows, n columns
// b  - right side (vector of m)
// x  - output (vector of m)
// isTransposed 0 => A transposed in memory, 0 => A not transposed in memory
// maxNNLSIters - maximum NNLS iterations, use (m+n)*2
// maxNNLSIters - maximum LS iterations, use (m+n)*2
// nSys - number of systems
// m  - number of equations
// n  - number of unknowns
// TOL_TERMINATION Tolerance for termination (0)

//NNLS  systems, no update/downdate
void nnls(REAL *A, REAL *b, REAL *x, int isTransposed, int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION);
  
//NNLS systems, update/downdate
void nnls_updates(REAL *A, REAL *b, REAL *x, int isTransposed,  int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION);

//NNLS systems, update/downdate single thread
void nnls_updates_single(REAL *A, REAL *b, REAL *x, int isTransposed,  int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION);

typedef void* NNLS_HANDLE;
//Allocate work memory for NNLS2 method to avoid repeated 
//nnls2_updates_single can run with n <= allocated n and/or m<= allocated m
NNLS_HANDLE allocate_nnls( int nSys, int m, int n);
void free_nnls(NNLS_HANDLE data);
//NNLS systems, update/downdate using pre-allocated storage
void nnls2_updates_single(NNLS_HANDLE data, REAL *A, REAL *b, REAL *x, int isTransposed,  int maxNNLSIters, int maxLSIters, int nSys, int m, int n, REAL TOL_TERMINATION);



#ifdef __cplusplus
}
#endif

#endif //__NNLS_H__