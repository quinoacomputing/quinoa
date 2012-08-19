// -----------------------------------------------------------------------------
// \file    src/LinearAlgebra/CG.C
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Conjugate gradient solver
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

//  Functions that implement the solution of a linear system with the method
//  of conjugate gradients, using a matrix stored in symmetric compressed
//  sparse row (CSR) format, with only the upper triangle stored including the
//  main diagonal.

#include <string.h>
#include <cstdio>
#include <cmath>
#include "Macros.h"
#include "Const.h"
#include "SparseMatrix.h"
#include "CG.h"

static void pc_SGS(sparsemat *A)
// -----------------------------------------------------------------------------
// Routine: pc_SGS - Prepare for symmetric Gauss-Seidel preconditioning
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i;

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( i = 0; i < A->rsize; i++ )
    A->d[i] = getma(A,i,i);		// extract diagonal
}

static void pc_Jacobi(sparsemat *A)
// -----------------------------------------------------------------------------
// Routine: pc_Jacobi - Extract main diagonal of sparsematrix A into its vector
//                      d that will be the Jacobi preconditioner
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i;

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( i = 0; i < A->rsize; i++ )
    A->d[i] = getma(A,i,i);
}

static void pc_none(sparsemat *A)
// -----------------------------------------------------------------------------
// Routine: pc_none - Fill vector d of sparsematrix A with 1.0,
//                    this means no preconditioner during conjugate gradients
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i;

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( i = 0; i < A->rsize; i++ )
    A->d[i] = 1.0;
}

static int cg_pcnone(sparsemat *A, double *x, double normb)
// -----------------------------------------------------------------------------
// Routine: cg_pcnone - CG solver with no preconditioner
// Author : J. Bakosi
// -----------------------------------------------------------------------------
// Solves linear system A * x = b with the conjugate gradient method without
// preconditioning.
//
// Input parameters:
//   A     - a symmetric positive definite sparse matrix in compressed sparse row
//           (CSR) format,
//   normb - norm of right hand side vector b
//
// Output parameters:
//   x     - the solution,
//   returns the number of iterations performed
// -----------------------------------------------------------------------------
{
  int i, j, it;
  double rho, rhoo, alpha, dtmp;

  // perform conjugate gradients iteration without preconditioner
  it = 0;
  do
  {
    rho = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:rho)
    #endif
    for ( i = 0; i < A->rsize; i++ )
       rho += A->r[i] * A->r[i];			// rho = (r,r)

    if ( it == 0 ) alpha = 0.0; else alpha = rho/rhoo;
    rhoo = rho;
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( i = 0; i < A->rsize; i++ )		// p = r + alpha * p
      A->p[i] = alpha * A->p[i] + A->r[i];

    dtmp = 0;
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j) reduction(+:dtmp)
    #endif
    for ( i = 0; i < A->rsize; i++ )
    {
      A->q[i] = 0.0;
      for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )// q = A * p
         A->q[i] += A->a[j] * A->p[A->ja[j]-1];
      dtmp += A->p[i] * A->q[i];		// dtmp = (p,q)
    }
    alpha = rho / dtmp;				// alpha = rho / (p,q)

    dtmp = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:dtmp)
    #endif
    for ( i = 0; i < A->rsize; i++ )
    {
       A->r[i] -= alpha * A->q[i];		// r = r - alpha * q
       dtmp += A->r[i] * A->r[i];		// compute norm of residual
    }
    dtmp = sqrt(dtmp);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( i = 0; i < A->rsize; i++ )		// x = x + alpha * p
      x[i] += alpha * A->p[i];

    it++;					// increase iteration counter

  } while ( (it < MAX_IT) && (dtmp > STOP_TOL*normb) );


  return( it );	// return number of iterations performed
}

static int cg_pcJacobi( sparsemat *A, double *x, double normb )
// -----------------------------------------------------------------------------
// Routine: cg_pcJacobi - CG solver with Jacobi preconditioner
// Author : J. Bakosi
// -----------------------------------------------------------------------------
// Solves linear system A * x = b with the conjugate gradient method with Jacobi
// preconditioning.
//
// Input parameters:
//   A    - a symmetric positive definite sparse matrix in compressed sparse row
//          (CSR) format,
//  normb - norm of right hand side vector b
//
// Output parameters:
//   x - the solution,
//   returns the number of iterations performed
// -----------------------------------------------------------------------------
{
  int i, j, it;
  double rho, rhoo, alpha, dtmp;

  // perform conjugate gradients iteration with Jacobi preconditioner
  it = 0;
  do
  {
    rho = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:rho)
    #endif
    for ( i = 0; i < A->rsize; i++ )
    {
       A->z[i] = A->r[i] / A->d[i];	// apply Jacobi pc: solve D * z = r
       rho += A->r[i] * A->z[i];	// rho = (r,z)
    }

    if ( it == 0 ) alpha = 0.0; else alpha = rho/rhoo;
    rhoo = rho;
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( i = 0; i < A->rsize; i++ )	// p = z + alpha * p
      A->p[i] = alpha * A->p[i] + A->z[i];

    dtmp = 0;
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j) reduction(+:dtmp)
    #endif
    for ( i = 0; i < A->rsize; i++ )
    {
      A->q[i] = 0.0;
      for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )	// q = A * p
         A->q[i] += A->a[j] * A->p[A->ja[j]-1];
      dtmp += A->p[i] * A->q[i];			// dtmp = (p,q)
    }
    alpha = rho / dtmp;					// alpha = rho / (p,q)

    dtmp = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:dtmp)
    #endif
    for ( i = 0; i < A->rsize; i++ )
    {
       A->r[i] -= alpha * A->q[i];		// r = r - alpha * q
       dtmp += A->r[i] * A->r[i];		// compute norm of residual
    }
    dtmp = sqrt(dtmp);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( i = 0; i < A->rsize; i++ )		// x = x + alpha * p
      x[i] += alpha * A->p[i];

    it++;					// increase iteration counter

  } while ( (it < MAX_IT) && (dtmp > STOP_TOL*normb) );


  return( it );	// return number of iterations performed
}

static int cg_pcSGS( sparsemat *A, double *x, double normb )
// -----------------------------------------------------------------------------
// Routine: cg_pcSGS - CG solver with symmetric Gauss-Seidel preconditioner
// Author : J. Bakosi
// -----------------------------------------------------------------------------
// Solves linear system A * x = b with the conjugate gradient method with
// symmetric Gauss-Seidel preconditioning using level scheduling.
//
// Input parameters:
//   A    - a symmetric positive definite sparse matrix in compressed sparse row
//          (CSR) format,
//  normb - norm of right hand side vector b
//
// Output parameters:
//   x    - the solution,
//   returns the number of iterations performed
// -----------------------------------------------------------------------------
{
  int i, j, k, l, it;
  double rho, rhoo, alpha, dtmp;

  // //////////////////////////////////////////////
  // THIS NEEDS TESTING! In the parallel forward and backwards solve with level
  // scheduling variable 'j' was declared 'shared' and it probably should be
  // private. That could well be the reason why it didn't scale at all. Now it's
  // fixed, but it still needs to be tested whether it gives the right solution
  // and also to see if it scales now.
  // //////////////////////////////////////////////

  // perform conjugate gradients iteration with SGS preconditioner
  it = 0;
  do
  {
    // apply SGS pc: solve M * z = r, where M = (D+L) * D^(-1) * (D+L^T)
    /*for ( i = 0; i < A->rsize; i++ )			// forward solve
    {
      dtmp = 0.0;
      for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )
        if ( A->ja[j] <= i ) dtmp += A->a[j] * A->u[A->ja[j]-1];
      A->u[i] = (A->r[i]-dtmp) / A->d[i];
    }*/
    for ( l = 0; l < A->nlev; l++ )			// forward solve
      #ifdef _OPENMP
      #pragma omp parallel for private(k,i,j,dtmp)
      #endif
      for ( k = A->lev2[l]; k < A->lev2[l+1]; k++ )
      {
	i = A->lev1[k];
        dtmp = 0.0;
        for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )
          if ( A->ja[j] <= i ) dtmp += A->a[j] * A->u[A->ja[j]-1];
        A->u[i] = (A->r[i]-dtmp) / A->d[i];
      }

    /*for ( i = A->rsize-1; i >= 0; i-- )		// backward solve
    {
      dtmp = 0.0;
      for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )
        if ( A->ja[j] > i+1 ) dtmp += A->a[j] * A->z[A->ja[j]-1];
      A->z[i] = A->u[i] - dtmp / A->d[i];
    }*/
    for ( l = A->nlev-1; l >= 0; l-- )			// backward solve
      #ifdef _OPENMP
      #pragma omp parallel for private(k,i,j,dtmp)
      #endif
      for ( k = A->lev2[l]; k < A->lev2[l+1]; k++ )
      {
	i = A->lev1[k];
        dtmp = 0.0;
        for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )
          if ( A->ja[j] > i+1 ) dtmp += A->a[j] * A->z[A->ja[j]-1];
        A->z[i] = A->u[i] - dtmp / A->d[i];
      }

    rho = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:rho)
    #endif
    for ( i = 0; i < A->rsize; i++ )
       rho += A->r[i] * A->z[i];			// rho = (r,z)

    if ( it == 0 ) alpha = 0.0; else alpha = rho/rhoo;
    rhoo = rho;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( i = 0; i < A->rsize; i++ )			// p = z + alpha * p
      A->p[i] = alpha * A->p[i] + A->z[i];
    
    dtmp = 0;
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j) reduction(+:dtmp)
    #endif
    for ( i = 0; i < A->rsize; i++ )
    {
      A->q[i] = 0.0;
      for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )	// q = A * p
         A->q[i] += A->a[j] * A->p[A->ja[j]-1];
      dtmp += A->p[i] * A->q[i];			// dtmp = (p,q)
    }
    alpha = rho / dtmp;					// alpha = rho / (p,q)
    
    dtmp = 0;
    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:dtmp)
    #endif
    for ( i = 0; i < A->rsize; i++ )
    {
       A->r[i] -= alpha * A->q[i];		// r = r - alpha * q
       dtmp += A->r[i] * A->r[i];		// compute norm of residual
    }
    dtmp = sqrt(dtmp);
    
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( i = 0; i < A->rsize; i++ )		// x = x + alpha * p
      x[i] += alpha * A->p[i];
    
    it++;					// increase iteration counter

  } while ( (it < MAX_IT) && (dtmp > STOP_TOL*normb) );

  return( it );	// return number of iterations performed
}

int cg_solve( sparsemat *A, double *b, double *x, int pc )
// -----------------------------------------------------------------------------
// Routine: cg_solve - CG solver
// Author : J. Bakosi
// -----------------------------------------------------------------------------
// Solves linear system A * x = b with the conjugate gradient method.
//
// Input parameters:
//   A  - a symmetric positive definite sparse matrix in compressed sparse row
//        (CSR) format, both triangles stored
//   b  - the right hand side,
//   pc - 0 : don't use preconditioner,
//        1 : use Jacobi preconditioner,
//        2 : use symmetric Gauss-Seidel preconditioner
//
// Output parameters:
//   x  - the solution,
//   returns the number of iterations needed to converge.
// -----------------------------------------------------------------------------
{
  int i, j;
  double normb;

  // compute initial residual and norm of right hand side
  normb = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for private(i,j) reduction(+:normb)
  #endif
  for ( i = 0; i < A->rsize; i++ )			
  {
    normb += b[i] * b[i];		// compute norm of right hand side
    A->r[i] = b[i];			// r = b - A * x
    for ( j = A->ia[i]-1; j < A->ia[i+1]-1; j++ )
       A->r[i] -= A->a[j] * x[A->ja[j]-1];
  }
  normb = sqrt( normb );


  // do cg iterations with different preconditioners
  switch ( pc )
  {
    case PC_NONE :
       pc_none(A);			// don't use preconditioner
       j = cg_pcnone(A,x,normb);	// do cg iteration without pc
       break;
    case PC_JACOBI :
       pc_Jacobi(A);			// extract Jacobi preconditioner
       j = cg_pcJacobi(A,x,normb);	// do cg iteration with Jacobi pc
       break;
    case PC_SGS :
       pc_SGS(A);			// prepare for symmetric Gauss-Seidel pc
       j = cg_pcSGS(A,x,normb);		// do cg iteration with SGS pc
       break;
  }


  // return the number of iterations performed
  return( j );
}
