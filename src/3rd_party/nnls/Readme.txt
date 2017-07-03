Author: 	Yuancheng Luo
E-mail: 	yluo1@mail.umd.edu
Date: 	2011-02-21
Modified: 	2011-02-21

NNLS OMP/MKL v1.1
v1.0 to 1.1
-added double precision support
-fixed a bug in zRemoved
================================================================
Introduction:
Non-negative least squares with updating/downdating on multi-core processors using OpenMP/MKL:

The package contains code for solving a set of Ax=b systems on multi-core machines. The work is based on the "active-set" method originally proposed by Lawson and Hanson (1974). A mex interface is also provided for use on Mathworks Matlab.

Two implementations are given: The first is a port of Matlab's lsqnonneg without the matrix updating and downdating steps for solving unconstrained sub-problems. The second uses a column-reordered update (Modified Gram-Schmidt) and downdate (Given's Rotations) step for modifying Q and R matrices. Both implementations make use of OpenMP for solving multiple linear systems and the Intel Math kernel Library (MKL) for BLAS / LAPACK routines.

For a link to the home-page and paper references, visit
http://www.cs.umd.edu/~yluo1/Projects/NNLS.html

================================================================
Compilation/Testing:
Compiled under intel icc 11.1 or gcc 4.5.1 via make
Tested under Dual Quad-Core Intel(R) Xeon(R) X5560 CPU @ 2.80GHz, Mathwork Matlab 2010b, 2009b, 2008b

For matlab, NNLS.mexa64 file is generated in the matlab folder.
For a standalone executable, runNNLS is generated in the current folder.

To generate alternate between single and double precision builds, comment out the definition USE_DOUBLE in the file headers.h.

================================================================
Matlab Interface:
x = NNLS(method, A, b, nSys, isTransposed, mEquations, nUnknowns, tol, MKLThreads, OMPThreads)

method: 0, 1 for without/with updating/downdating
A, b: Single or double precision matrix/array in column major order
nSys: Number of systems to solve
isTransposed: A is transposed in memory
mEquations: Number of equations (rows in system A)
nUnknowns: Number of variables (cols in system A)
tol: Tolerance used in termination condition and 0
MKLThreads: Number of MLK threads to use
OMPThreads: Number of OpenMP threads to use

x: output array

For a simple test between lsqnoneg and NNLS.mexa64, see matlab/randomNNLS.m

================================================================
Executable interface:
Command line:  method, nSys, isTransposed, mEquations, nUnknowns, tol, MKLThreads, OMPThreads

Input Files:
sysA.bin: A binary input file containing nSys matrices A stored in column major order
sysB.bin: A binary input file containing nSys arrays b

Output Files:
sysX.bin: A binary output file containing the non-negative solution Ax=b
cpu_results.txt: An ascii output of the system and equations

================================================================