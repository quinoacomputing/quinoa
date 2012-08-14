//  ------------------------------------------------------------------------------------------------------------
//
//  Copyright 2007 Jozsef Bakosi
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  ------------------------------------------------------------------------------------------------------------
//
//  Functions that implement elliptic relaxation in 2d.
//  For more info see main.cc
//


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mkl.h"
#include "const.h"
#include "macros.h"
#include "sparsemat.h"
#include "elliptic_relaxation.h"
#include "cg.h"



#ifndef WALLFUNCTIONS




static void build_lhs( int nelem, int nbpoin, int nthreads, int *bpg, int *inpoel, int *bpsup1, int *bpsup2,
                       int *bptags, double *Ae, double *dNx, double *dNy, double *u2, double *f,
		       double *wnz, sparsemat *R )
//
// builds lhs for elliptic relaxation step
//
{
  int e, eN, A, B, C, UA, UB, UC, i, j, k, cw;
  double aL2, dtmp;
  double tke[NNODE], l1[NNODE], l2[NNODE], L[NNODE], eps[NNODE];


  // zero matrix R
  dpzero( R->a, R->nnz, nthreads );


  // build lhs of elliptic relaxation step...
  #ifdef _OPENMP
  #pragma omp parallel for private(e,eN,A,B,C,UA,UB,UC,tke,eps,l1,l2,i,L,aL2,dtmp)
  #endif
  for ( e = 0; e < nelem; e++ )		// loop over all elements
  {
    eN = e*NNODE;
    A = inpoel[eN+0];			// get indices of nodes
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;

    // nodal turbulent kinetic energy in element e
    tke[0] = 0.5*(u2[UA+0] + u2[UA+1] + u2[UA+2]);
    tke[1] = 0.5*(u2[UB+0] + u2[UB+1] + u2[UB+2]);
    tke[2] = 0.5*(u2[UC+0] + u2[UC+1] + u2[UC+2]);
    if (tke[0]<EPS) tke[0]=BOUND;
    if (tke[1]<EPS) tke[1]=BOUND;
    if (tke[2]<EPS) tke[2]=BOUND;
    // rate of dissipation of turbulent kinetic energy in element e
    eps[0] = f[A]*(tke[0] + CT2*f[A]/RE);
    eps[1] = f[B]*(tke[1] + CT2*f[B]/RE);
    eps[2] = f[C]*(tke[2] + CT2*f[C]/RE);
    if (eps[0]<EPS) eps[0]=BOUND;
    if (eps[1]<EPS) eps[1]=BOUND;
    if (eps[2]<EPS) eps[2]=BOUND;
    // nodal turbulent-, and Kolmogorov-lengthscales in element e
    l1[0] = sqrt(tke[0]*tke[0]*tke[0])/eps[0];
    l1[1] = sqrt(tke[1]*tke[1]*tke[1])/eps[1];
    l1[2] = sqrt(tke[2]*tke[2]*tke[2])/eps[2];
    l2[0] = 1.0/sqrt(sqrt(eps[0]*RE3));
    l2[1] = 1.0/sqrt(sqrt(eps[1]*RE3));
    l2[2] = 1.0/sqrt(sqrt(eps[2]*RE3));

    // precompute Laplacian with reusing arrays tke and eps
    tke[0] = dNx[eN+0]*dNx[eN+0] + dNy[eN+0]*dNy[eN+0];	// main diagonal
    tke[1] = dNx[eN+1]*dNx[eN+1] + dNy[eN+1]*dNy[eN+1];
    tke[2] = dNx[eN+2]*dNx[eN+2] + dNy[eN+2]*dNy[eN+2];
    eps[0] = dNx[eN+0]*dNx[eN+1] + dNy[eN+0]*dNy[eN+1];	// off-diagonal
    eps[1] = dNx[eN+0]*dNx[eN+2] + dNy[eN+0]*dNy[eN+2];
    eps[2] = dNx[eN+1]*dNx[eN+2] + dNy[eN+1]*dNy[eN+2];

    for ( i = 0; i < RDOF; i++ )
    {
      // nodal characteristic lengthscales in element e
      //
      // THIS IS PROBABLY NOT GENERAL ENOUGH !!!
      //
      if (i==4) dtmp = CXI; else dtmp = 1.0;
      L[0] = CL*fmax(dtmp*l1[0], CETA*l2[0]);
      L[1] = CL*fmax(dtmp*l1[1], CETA*l2[1]);
      L[2] = CL*fmax(dtmp*l1[2], CETA*l2[2]);
      // element-average lengthscale square
      aL2 = (L[0]+L[1]+L[2])/3.0;
      aL2 *= aL2;
    
      // main diagonal
      addmr(R,A,A,i,Ae[e]*(1.0/6.0 + aL2*tke[0]));	// block: (1,1)
      addmr(R,B,B,i,Ae[e]*(1.0/6.0 + aL2*tke[1]));	// block: (2,2)
      addmr(R,C,C,i,Ae[e]*(1.0/6.0 + aL2*tke[2]));	// block: (3,3)
      
      // off-diagonal
      dtmp = Ae[e]*(1.0/12.0 + aL2*eps[0]);
      addmr(R,A,B,i,dtmp);				// block: (1,2)
      addmr(R,B,A,i,dtmp);				// block: (2,1)
      dtmp = Ae[e]*(1.0/12.0 + aL2*eps[1]);
      addmr(R,A,C,i,dtmp);				// block: (1,3)
      addmr(R,C,A,i,dtmp);				// block: (3,1)
      dtmp = Ae[e]*(1.0/12.0 + aL2*eps[2]);
      addmr(R,B,C,i,dtmp);			 	// block: (2,3)
      addmr(R,C,B,i,dtmp);				// block: (3,2)
    }
  }


  // save nonzeros in columns of wall points for later use
  for ( cw=i=0; i < nbpoin; i++ )			// loop over all boundary points
    if (bptags[i*3+1] == 4)	// wall
      for ( j = bpsup2[i]+1; j <= bpsup2[i+1]; j++ )	// loop over nonzeros in block row bpg[i]
        for ( k = 0; k < RDOF; k++ )			// loop over all degrees of freedom
	  wnz[cw++] = getmr(R,bpg[i],bpsup1[j],k);	// save matrix value, increase wall counter
  
  
  #ifdef _OPENMP
  #pragma omp parallel for private(i,k,j)
  #endif
  // modify matrix for Dirichlet conditions for wall
  for ( i = 0; i < nbpoin; i++ )			// loop over all boundary points
    if (bptags[i*3+1] == 4)	// wall
    {
       // zero column
       for ( j = bpsup2[i]+1; j <= bpsup2[i+1]; j++ )	// loop over nonzeros in block row bpg[i]
         for ( k = 0; k < RDOF; k++ )			// loop over all degrees of freedom
           insmr(R,bpsup1[j],bpg[i],k,0.0);		// insert zero
       
       // zero row
       for ( k = 0; k < RDOF; k++ )			// loop over all degrees of freedom
         for ( j = R->ia[bpg[i]*RDOF+k]-1; j < R->ia[bpg[i]*RDOF+k+1]-1; j++ )
           R->a[j] = 0.0;

       // insert 1.0 into main diagonal
       for ( k = 0; k < RDOF; k++ )			// loop over all degrees of freedom
         insmr(R,bpg[i],bpg[i],k,1.0);			// insert zero
    }
  
  
  #ifdef _OPENMP
  #pragma omp parallel for private(i,k,j)
  #endif
  // modify matrix for Dirichlet conditions for top
  // homogeneous Dirichlet for off-diagonal (diagonal remains homogeneous Neumann)
  for ( i = 0; i < nbpoin; i++ )			// loop over all boundary points
    if (bptags[i*3+1] != 4)	// not wall
    {
       // zero column
       for ( j = bpsup2[i]+1; j <= bpsup2[i+1]; j++ )	// loop over nonzeros in block row bpg[i]
         for ( k = 1; k < RDOF-1; k++ ) if (k!=4)	// touch only off-diagonal
           insmr(R,bpsup1[j],bpg[i],k,0.0);		// insert zero
       
       // zero row
       for ( k = 1; k < RDOF-1; k++ ) if (k!=4)		// touch only off-diagonal
         for ( j = R->ia[bpg[i]*RDOF+k]-1; j < R->ia[bpg[i]*RDOF+k+1]-1; j++ )
           R->a[j] = 0.0;
       
       // insert 1.0 into main diagonal
       for ( k = 1; k < RDOF-1; k++ ) if (k!=4)		// touch only off-diagonal
         insmr(R,bpg[i],bpg[i],k,1.0);
    }
}










int rstep( int nelem, int nbpoin, int nthreads,
           int *bpg, int *inpoel, int *wlc, int *bpsup1, int *bpsup2, int *bptags,
	   double *Ae, double *wpnr, double *rrhs, double *rho,
	   double *dNx, double *dNy, double *u2, double *f, double *wnz, double *du,
	   sparsemat *R )
//
// elliptic relaxation step
//
// returns the number of iterations needed to converge for the linear solve
//
{
  int e, eN, A, B, C, UA, UB, UC, RA, RB, RC, i, j, k, cw;
  double ak, ak2, af, det, dtmp, Av, val;
  double tke[NNODE], au2[U2DOF], ba[U2DOF], adu[NDIMN*NDIMN], nr[NDIMN], dc[RDOF];


  // re-build lhs for elliptic relaxation with new L
  build_lhs( nelem, nbpoin, nthreads, bpg, inpoel, bpsup1, bpsup2, bptags,
             Ae, dNx, dNy, u2, f, wnz, R );

  // zero rhs
  dpzero( rrhs, R->rsize, nthreads );

  
  // build rhs for elliptic relaxation step...
  #ifdef _OPENMP
  #pragma omp parallel for private(e,eN,j,A,B,C,UA,UB,UC,RA,RB,RC,tke,ak,ak2,af,au2,ba,det,Av,adu,dtmp)
  #endif
  for ( e = 0; e < nelem; e++ )
  {
    // get node-information for element e...
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;

    // nodal turbulent kinetic energy
    tke[0] = 0.5*(u2[UA+0] + u2[UA+1] + u2[UA+2]);
    tke[1] = 0.5*(u2[UB+0] + u2[UB+1] + u2[UB+2]);
    tke[2] = 0.5*(u2[UC+0] + u2[UC+1] + u2[UC+2]);
    
    // element-average turbulent kinetic energy
    ak = (tke[0] + tke[1] + tke[2])/3.0;
    if (ak<EPS) ak=BOUND;
    ak2 = 2.0*ak;

    // element-average mean frequency
    af = (f[A]+f[B]+f[C])/3.0;

    // element-average Reynolds stress
    #ifdef __INTEL_COMPILER
    #pragma vector always
    #endif
    for ( j = 0; j < 6; j++ )
      au2[j] = (u2[UA+j]+u2[UB+j]+u2[UC+j])/3.0;

    // element-averge Reynolds stress anisotropy
    #ifdef __INTEL_COMPILER
    #pragma vector always
    #endif
    for ( j = 0; j < 6; j++ )
      ba[j] = au2[j]/ak2;
    ba[0] -= 1.0/3.0;
    ba[1] -= 1.0/3.0;
    ba[2] -= 1.0/3.0;
    
    // element-average determinant of Reynolds stress tensor
    det = au2[0]*(au2[1]*au2[2] - au2[5]*au2[5])
         -au2[3]*(au2[3]*au2[2] - au2[4]*au2[5])
         +au2[4]*(au2[3]*au2[5] - au2[4]*au2[1]);
    dtmp = CV*det*3.375/(ak*ak*ak);
    Av = fmin(1.0,dtmp);

    // mean velocity derivatives
    #ifdef __INTEL_COMPILER
    #pragma vector always
    #endif
    for ( j = 0; j < 4; j++ )
      adu[j] = (du[A*4+j] + du[B*4+j] + du[C*4+j])/3.0;
    
    // precompute tensor indices
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    // scatter-add rhs, wiggly-p_11
    dtmp = Ae[e]/3.0*ak*( 0.5*(1.0-C1)*af + C2*Av*adu[0] + GAMMA5*ba[3]*(adu[2]-adu[1]) );
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RA+0] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RB+0] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RC+0] += dtmp;

    // wiggly-p_12
    dtmp = Ae[e]/3.0*ak*( (C2*Av+GAMMA5/3.0)*adu[1] - GAMMA5/3.0*adu[2] + GAMMA5*ba[0]*(adu[1]-adu[2]) );
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RA+1] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RB+1] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RC+1] += dtmp;

    // wiggly-p_21
    dtmp = Ae[e]/3.0*ak*( (C2*Av+GAMMA5/3.0)*adu[2] - GAMMA5/3.0*adu[1] + GAMMA5*ba[1]*(adu[2]-adu[1]) );
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RA+3] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RB+3] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RC+3] += dtmp;

    // wiggly-p_22
    dtmp = Ae[e]/3.0*ak*( 0.5*(1.0-C1)*af + C2*Av*adu[3] + GAMMA5*ba[3]*(adu[1]-adu[2]) );
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RA+4] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RB+4] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RC+4] += dtmp;

    // wiggly-p_31
    dtmp = Ae[e]/3.0*ak*GAMMA5*ba[5]*(adu[2]-adu[1]);
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RA+6] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RB+6] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RC+6] += dtmp;

    // wiggly-p_32
    dtmp = Ae[e]/3.0*ak*GAMMA5*ba[4]*(adu[1]-adu[2]);
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RA+7] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RB+7] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RC+7] += dtmp;
    
    // wiggly-p_33
    dtmp = Ae[e]/3.0*ak*0.5*(1.0-C1)*af;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RA+8] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RB+8] += dtmp;
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    rrhs[RC+8] += dtmp;
  }


  // apply boundary conditions for R and RRHS...
  // applying inhomogeneous Dirichlet bcs:
  //   o find the boundary point index (i),
  //   o substract column (i) of the matrix multiplied by the value of the boundary condition (xi)
  //     from the rhs (ie. for loop over j: rhsj -= aji*xi),
  //   o put 1.0 in the main diagonal,
  //   o zero everything else in the row and column corresponding to the boundary point index (i),
  //   o put the boundary value (xi) into rhs (ie. ri = xi)
  //
  // wall: inhomogeneous Dirichlet for wiggly-p_ij (i,j=1,2), homogeneous for all other components
  for ( cw=i=0; i < nbpoin; i++ )				// loop over boundary points
    if (bptags[i*3+1] == 4) // wall
    {
      nr[0] = wpnr[wlc[bpg[i]]*NDIMN+0];			// get wall normal for wall point i
      nr[1] = wpnr[wlc[bpg[i]]*NDIMN+1];
      val = -4.5*f[bpg[i]]*f[bpg[i]]*CT2/RE;			// boundary value: -4.5*eps
      dc[0] = val*nr[0]*nr[0];					// build Dirichlet condition for tensor
      dc[1] = val*nr[0]*nr[1];
      dc[3] = val*nr[1]*nr[0];
      dc[4] = val*nr[1]*nr[1];
      dc[2]=dc[5]=dc[6]=dc[7]=dc[8]=0.0;

      // substract boundary_value*matrix_element from rhs
      for ( j = bpsup2[i]+1; j <= bpsup2[i+1]; j++ )		// loop over nonzeros in block row bpg[i]
        for ( k = 0; k < RDOF; k++ )				// loop over all degrees of freedom
	  rrhs[bpsup1[j]*RDOF+k] -= dc[k]*wnz[cw++];		// substract from rhs, increase wall counter

      // put boundary value into rhs
      rrhs[RDOF*bpg[i]+0] = val*nr[0]*nr[0];			// (0,0)
      rrhs[RDOF*bpg[i]+1] = val*nr[0]*nr[1];			// (0,1)
      rrhs[RDOF*bpg[i]+3] = val*nr[1]*nr[0];			// (1,0)
      rrhs[RDOF*bpg[i]+4] = val*nr[1]*nr[1];			// (1,1)
      rrhs[RDOF*bpg[i]+2] = rrhs[RDOF*bpg[i]+5] =
      rrhs[RDOF*bpg[i]+6] = rrhs[RDOF*bpg[i]+7] = rrhs[RDOF*bpg[i]+8] = 0.0;
   }
  
  
  // top: homogeneous Dirichlet for off-diagonal (diagonal remains homogeneous Neumann)
  #ifdef _OPENMP
  #pragma omp parallel for private(i,j)
  #endif
  for ( i = 0; i < nbpoin; i++ )	                       	// loop over boundary points
    if (bptags[i*3+1] != 4) // not wall
    {
      // put boundary value into rhs
      j = bpg[i]*RDOF;
      rrhs[j+1] = rrhs[j+2] = rrhs[j+3] =
      rrhs[j+5] = rrhs[j+6] = rrhs[j+7] = 0.0;
    }


  // solve linear system (R * rho = rrhs), use Jacobi preconditioner and return
  return( cg_solve(R,rrhs,rho,PC_JACOBI) );
}




#endif // WALLFUNCTIONS
