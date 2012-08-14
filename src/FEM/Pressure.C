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
//  Functions that implement the solution of the Poisson-equation for the mean pressure
//  for more info see main.cc
//

#include <math.h>
#include <stdio.h>
#include "mkl.h"
#include "const.h"
#include "macros.h"
#include "sparsemat.h"
#include "pressure.h"
#include "particles.h"
#include "cg.h"









void pinit( int nelem, int nbpoin, int nthreads, int *inpoel, int *bpg, int *bpsup1, int *bpsup2,
            double maxx, double *coord,
	    double *Ae, double *dNx, double *dNy, double *ofnz, sparsemat *P )
//
// initializes linear solver for mean pressure solve
// this is only called once during setup
//
{
  int e, eN, i, j, cof, A, B, C;
  double dtmp;


  // zero matrix P
  dpzero( P->a, P->nnz, nthreads );


  // build lhs for Poisson-equation
  #ifdef _OPENMP
  #pragma omp parallel for private(e,eN,A,B,C,dtmp)
  #endif
  for ( e = 0; e < nelem; e++ )		// loop over all elements
  {
    eN = e*NNODE;
    A = inpoel[eN+0];			// get indices of nodes
    B = inpoel[eN+1];
    C = inpoel[eN+2];

    // main diagonal
    addma(P,A,A,Ae[e]*(dNx[eN+0]*dNx[eN+0] + dNy[eN+0]*dNy[eN+0]));	// (1,1)
    addma(P,B,B,Ae[e]*(dNx[eN+1]*dNx[eN+1] + dNy[eN+1]*dNy[eN+1]));	// (2,2)
    addma(P,C,C,Ae[e]*(dNx[eN+2]*dNx[eN+2] + dNy[eN+2]*dNy[eN+2]));	// (3,3)

    // off-diagonal
    dtmp = Ae[e]*(dNx[eN+0]*dNx[eN+1] + dNy[eN+0]*dNy[eN+1]);
    addma(P,A,B,dtmp);							// (1,2)
    addma(P,B,A,dtmp);							// (2,1)
    dtmp = Ae[e]*(dNx[eN+0]*dNx[eN+2] + dNy[eN+0]*dNy[eN+2]);
    addma(P,A,C,dtmp);							// (1,3)
    addma(P,C,A,dtmp);							// (3,1)
    dtmp = Ae[e]*(dNx[eN+1]*dNx[eN+2] + dNy[eN+1]*dNy[eN+2]);
    addma(P,B,C,dtmp);							// (2,3)
    addma(P,C,B,dtmp);							// (3,2)
  }


  // save nonzero values of columns for outflow for later use
  for ( cof=i=0; i < nbpoin; i++ )			// loop over all boundary points
  {
     // extract nonzeros from column corresonding to outflow boundary point i from matrix P
     if ( fabs(coord[bpg[i]*NDIMN+0]-maxx)<EPS )	// outflow
       for ( j = P->ia[bpg[i]]-1; j < P->ia[bpg[i]+1]-1; j++ )// loop over nonzeros in row bpg[i]
 	 ofnz[cof++] = P->a[j];
  }


  // modify matrix for Dirichlet conditions for outflow
  #ifdef _OPENMP
  #pragma omp parallel for private(i,j)
  #endif
  for ( i = 0; i < nbpoin; i++ )			// loop over all boundary points
    if ( fabs(coord[bpg[i]*NDIMN+0]-maxx)<EPS )		// outflow
    {
       // zero column
       for ( j = bpsup2[i]+1; j <= bpsup2[i+1]; j++ )	// loop over nonzeros in row bpg[i]
         insma(P,bpsup1[j],bpg[i],0.0);
       
       // zero row
       for ( j = P->ia[bpg[i]]-1; j < P->ia[bpg[i]+1]-1; j++ )
         P->a[j] = 0.0;
       
       // insert 1.0 into main diagonal
       insma(P,bpg[i],bpg[i],1.0);
    }
}










int pstep( int nwe, int npoin, int nelem, int nbpoin, int nthreads, double dt, double maxx,
           int *inpoel, int *binpoel, int *we, int *weo, int *bpg, int *betags,
	   double *odpn, double *oudt, double *wel, double *wenr, double *prhs, double *du,
	   double *ddu, double *pr, double *dpr, double *u2,
	   double *dNx, double *dNy, double *Ae, double *coord, sparsemat *P, double *u )
//
// pressure-Poisson step
//
// returns the number of iterations needed for the linear solve to converge
//
{
  int e, eN, pN, p4, A, B, C, i, pit, UA, UB, UC;
  double dpn, divu;
  double f[NNODE], w[NNODE], nr[NDIMN];


  // zero rhs
  dpzero( prhs, npoin, nthreads );
  

  // build rhs for Poisson-equation...
  #ifdef _OPENMP
  #pragma omp parallel for private(e,eN,A,B,C,divu,f)
  #endif
  for ( e = 0; e < nelem; e++ )
  {
    // get node-information for element e
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
   
    // element-average divergence of <Ui> / 3.0 / dt
    divu = (du[A*4+0]+du[B*4+0]+du[C*4+0] + du[A*4+3]+du[B*4+3]+du[C*4+3])/9.0/dt;

    // artificial pressure diffusion for stabilization
    f[0] = dNx[eN+0] * (dNx[eN+0]*pr[A] + dNx[eN+1]*pr[B] + dNx[eN+2]*pr[C])
         + dNy[eN+0] * (dNy[eN+0]*pr[A] + dNy[eN+1]*pr[B] + dNy[eN+2]*pr[C]);
    f[1] = dNx[eN+1] * (dNx[eN+0]*pr[A] + dNx[eN+1]*pr[B] + dNx[eN+2]*pr[C])
         + dNy[eN+1] * (dNy[eN+0]*pr[A] + dNy[eN+1]*pr[B] + dNy[eN+2]*pr[C]);
    f[2] = dNx[eN+2] * (dNx[eN+0]*pr[A] + dNx[eN+1]*pr[B] + dNx[eN+2]*pr[C])
         + dNy[eN+2] * (dNy[eN+0]*pr[A] + dNy[eN+1]*pr[B] + dNy[eN+2]*pr[C]);

    // scatter-add element-level right hand side
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    prhs[A] -= Ae[e]*(divu + Cp*f[0]/dt);
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    prhs[B] -= Ae[e]*(divu + Cp*f[1]/dt);
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    prhs[C] -= Ae[e]*(divu + Cp*f[2]/dt);
  }


  // inhomogeneous Neumann condition for wall-elements
  for ( nwe=i=0; i < nbpoin; i++ )	// loop over all boundary elements
    if (betags[i*2+0] == 4) // wall
    {
       pN = NDIMN*binpoel[i*NBNODE+0];	// get the index of one of the nodes of boundary element i
       p4 = 4*binpoel[i*NBNODE+0];
       eN = we[nwe]*NNODE;		// get node information for domain-element corresponding to wall-element nwe
       A = inpoel[eN+0];
       B = inpoel[eN+1];
       C = inpoel[eN+2];
       UA = A*U2DOF;
       UB = B*U2DOF;
       UC = C*U2DOF;
       nr[0] = wenr[nwe*NDIMN+0];	// get normal of wall-element nwe
       nr[1] = wenr[nwe*NDIMN+1];

       // d<p>/dn = (\nu d^2<Ui>/dxjdxj - d<uiuj>/dxj - <Uj> d<Ui>/dxj - d<Ui>/dt) ni
       dpn = (ddu[pN+0]*nr[0]      // nu * (d^2<U>/dx^2 + d^2<U>/dy^2) nx
            + ddu[pN+1]*nr[1])/RE  // nu * (d^2<V>/dx^2 + d^2<V>/dy^2) ny
	    - (dNx[eN+0]*u2[UA+0] + dNx[eN+1]*u2[UB+0] + dNx[eN+2]*u2[UC+0] + // (d<uu>/dx + d<uv>/dy) nx
	       dNy[eN+0]*u2[UA+3] + dNy[eN+1]*u2[UB+3] + dNy[eN+2]*u2[UC+3] ) * nr[0]
	    - (dNx[eN+0]*u2[UA+3] + dNx[eN+1]*u2[UB+3] + dNx[eN+2]*u2[UC+3] + // (d<uv>/dx + d<vv>/dy) ny
	       dNy[eN+0]*u2[UA+1] + dNy[eN+1]*u2[UB+1] + dNy[eN+2]*u2[UC+1] ) * nr[1]
            - (u[pN+0]*du[p4+0] + u[pN+1]*du[p4+1]) * nr[0]  // (<U> d<U>/dx + <V> d<U>/dy) nx
            - (u[pN+0]*du[p4+2] + u[pN+1]*du[p4+3]) * nr[1]  // (<U> d<V>/dx + <V> d<V>/dy) ny
            - oudt[nwe]		// d<Ui>/dt ni
            ;

       // depending on which side of the element is aligned with the wall the integral contributes
       // differently to nodes: node weo[nwe] is not on the wall so it doesn't contribute
       w[0]=w[1]=w[2]=0.5;
       w[weo[nwe]] = 0.0;

       // scatter-add boundary integral contributions
       prhs[A] += wel[nwe]*w[0]*((dpn-odpn[nwe]) + Cp*odpn[nwe]);
       prhs[B] += wel[nwe]*w[1]*((dpn-odpn[nwe]) + Cp*odpn[nwe]);
       prhs[C] += wel[nwe]*w[2]*((dpn-odpn[nwe]) + Cp*odpn[nwe]);

       nwe++;				// increase wall-element counter
    }

  
  // homogeneous Dirichlet at outflow:
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( i = 0; i < nbpoin; i++ )			// loop over all boundary points
    if ( fabs(coord[bpg[i]*NDIMN+0]-maxx)<EPS )		// outflow
      prhs[bpg[i]] = 0.0;

  // initial guess for dpr = 0
  dpzero( dpr, npoin, nthreads );


  // solve linear system: P * dpr = prhs, use Jacobi preconditioner
  pit = cg_solve(P,prhs,dpr,PC_JACOBI);


  // update mean pressure
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( i = 0; i < npoin; i++ ) pr[i] += dpr[i];

  return( pit );
}







void velcorr( int npar, double dt, int *inpoel, int *elp, double *parvel, double *dpr, double *dNx, double *dNy )
//
// corrects velocity of particles after pressure projection
//
{
  int p, p3, eN, A, B, C;

  
  // correct particle velocities
  #ifdef _OPENMP
  #pragma omp parallel for private(p,p3,eN,A,B,C)
  #endif
  for ( p = 0; p < npar; p++ )
  {
    p3 = p*3;
    eN = elp[p]*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];

    // correct particle velocities
    parvel[p3+0] -= dt*(dNx[eN+0]*dpr[A] + dNx[eN+1]*dpr[B] + dNx[eN+2]*dpr[C]);
    parvel[p3+1] -= dt*(dNy[eN+0]*dpr[A] + dNy[eN+1]*dpr[B] + dNy[eN+2]*dpr[C]);
  }
}






void save_old_dpn( int nwe, int nbpoin, int *binpoel, int *betags, double *wenr,
                   double *ddu, double *odpn, double *oudt,
                   int *we, int *inpoel, double *u2, double *dNx, double *dNy,
                   double *u, double *du )
//
// saves the normal derivative of mean pressure (d<p>/dn) in wall-points for next timestep
//
{
  int i, pN, p4, eN, A, B, C, UA, UB, UC;
  double nr[NDIMN];


  for ( nwe=i=0; i < nbpoin; i++ )	// loop over all boundary elements
    if (betags[i*2+0] == 4) // wall
    {
       pN = NDIMN*binpoel[i*NBNODE+0];	// get the index of one of the nodes of boundary element i
       p4 = 4*binpoel[i*NBNODE+0];
       eN = we[nwe]*NNODE;		// get node information for domain-element corresponding to wall-element nwe
       A = inpoel[eN+0];
       B = inpoel[eN+1];
       C = inpoel[eN+2];
       UA = A*U2DOF;
       UB = B*U2DOF;
       UC = C*U2DOF;
       nr[0] = wenr[nwe*NDIMN+0];	// get normal of wall-element nwe
       nr[1] = wenr[nwe*NDIMN+1];

       // d<p>/dn = (\nu d^2<Ui>/dxjdxj - d<uiuj>/dxj - <Uj> d<Ui>/dxj - d<Ui>/dt) ni
       odpn[nwe] = (ddu[pN+0]*nr[0]      // nu * (d^2<U>/dx^2 + d^2<U>/dy^2) nx
                  + ddu[pN+1]*nr[1])/RE  // nu * (d^2<V>/dx^2 + d^2<V>/dy^2) ny
	          - (dNx[eN+0]*u2[UA+0] + dNx[eN+1]*u2[UB+0] + dNx[eN+2]*u2[UC+0] + // (d<uu>/dx + d<uv>/dy) nx
	             dNy[eN+0]*u2[UA+3] + dNy[eN+1]*u2[UB+3] + dNy[eN+2]*u2[UC+3] ) * nr[0]
	          - (dNx[eN+0]*u2[UA+3] + dNx[eN+1]*u2[UB+3] + dNx[eN+2]*u2[UC+3] + // (d<uv>/dx + d<vv>/dy) ny
	             dNy[eN+0]*u2[UA+1] + dNy[eN+1]*u2[UB+1] + dNy[eN+2]*u2[UC+1] ) * nr[1]
                  - (u[pN+0]*du[p4+0] + u[pN+1]*du[p4+1]) * nr[0]  // (<U> d<U>/dx + <V> d<U>/dy) nx
                  - (u[pN+0]*du[p4+2] + u[pN+1]*du[p4+3]) * nr[1]  // (<U> d<V>/dx + <V> d<V>/dy) ny
                  - oudt[nwe]		// d<Ui>/dt ni
                  ;

       nwe++;				// increase wall-element counter
    }
}
