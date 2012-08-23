// -----------------------------------------------------------------------------
// \file    src/SDE/Langevin.C
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Langevin equation
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------
//
//  Functions that implement:
//   * the generalized Langevin model with elliptic relaxation,
//   * the simplified Langevin model with wall-functions,
//   * the IEM/IECM micromixing model.
//
//  Contains the main timestepping loop, the main loop that updates particles
//  and implements particle boundary conditions, in essence, this file contains
//  the core of the whole program.
//
// -----------------------------------------------------------------------------

#include <cstdio>
#include <string.h>
#include <cmath>
#ifdef _OPENMP
#include "omp.h"
#endif
#include "mkl.h"
#include "Const.h"
#include "Macros.h"
#include "SparseMatrix.h"
#include "Matrix3.h"
#include "Random.h"
#include "Langevin.h"
#include "TimeStepping.h"
#include "Particles.h"
#include "Postprocess.h"
#include "Pressure.h"
#include "EllipticRelaxation.h"
#ifndef WALLFUNCTIONS
#include "RandomErrors.h"
#endif

/*
#ifdef MICROMIXING
static double dist_from_source( int p, double *parcoord )
// -----------------------------------------------------------------------------
// Routine: double_dist_from_source - Compute distance of particle p from source
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  register int p2 = p*2;

  return( sqrt((parcoord[p2+0]-XS)*(parcoord[p2+0]-XS) +
               (parcoord[p*2+1]-YS)*(parcoord[p*2+1]-YS)) );
}
#endif
*/

#ifdef WALLFUNCTIONS
static void parwallfunctionbc(int p, double tke, double af, double *au,
                              double *au2, double *dp, double *rBA, double *nro,
                              double *parvel, double *parfreq)
// -----------------------------------------------------------------------------
// Routine: parwallfunctionbc - Enforce wall boundary condition on velocity and
//                              frequency of particle p using wall-functions
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int p3;
  double ui, vi, ur, vr, up, sgn, gammat, us, ue, a, dpw0;
  double nr[NDIMN], /*dpw[NDIMN],*/ auw[NDIMN], au2w[9], m[9], m2[9], t[9];
 
  p3 = p*3;

  // we will operate with inward-pointing normal
  nr[0] = -nro[0];
  nr[1] = -nro[1];

  // compute incident particle velocity in wall-coordinate system
  ui = parvel[p3+0]*rBA[0] + parvel[p3+1]*rBA[1];       // streamwise component
  vi = parvel[p3+0]*nr[0] + parvel[p3+1]*nr[1];         // wall-normal component

  // transform mean pressure gradient into wall-coordinate system
  dpw0 = dp[0]*rBA[0] + dp[1]*rBA[1];
  //dpw[1] = dp[0]*nr[0] + dp[1]*nr[1];

  // transform mean velocity into wall-coordinate system
  auw[0] = au[0]*rBA[0] + au[1]*rBA[1];
  auw[1] = au[0]*nr[0] + au[1]*nr[1];

  // transform Reynolds stress tensor into wall-coordinate system
  m[0] = au2[0];  m[1] = au2[3];  m[2] = au2[4];   // copy out symmetric tensor
  m[3] = au2[3];  m[4] = au2[1];  m[5] = au2[5];
  m[6] = au2[4];  m[7] = au2[5];  m[8] = au2[2];
  t[0] = rBA[0];  t[1] = rBA[1];  t[2] = 0.0;// construct transformation matrix
  t[3] = nr[0];   t[4] = nr[1];   t[5] = 0.0;
  t[6] = 0.0;     t[7] = 0.0;     t[8] = 1.0;
  m3mult(t, m, m2);       // transform Reynolds stress tensor
  m3mult(t, m2, au2w);

  // compute reflected particle velocity in wall-coordinate system
  vr = -vi;     // vr
  up = sqrt(sqrt(CMU))*sqrt(tke);
  if ( au2w[1]*dpw0 > 0.0 ) sgn = 1.0; else sgn = -1.0;
  gammat = fmax(0.0,sgn);
  us = sqrt( sqrt(CMU)*tke + gammat*fabs(YP*dpw0) );
  ue = us/KAPPA * log(E*YP*us*RE);
  a = 2.0*up*up*auw[0]*fabs(auw[0])/(au2w[4]*ue*ue);
  ur = ui + a*vi;       // ur

  // assign reflected particle velocity (inverse-transform ur and vr)
  a = rBA[0]*nr[1] - rBA[1]*nr[0]; // determinant of the transformation matrix
  parvel[p3+0] = (nr[1]*ur - rBA[1]*vr)/a;      // x component
  parvel[p3+1] = (rBA[0]*vr - nr[0]*ur)/a;      // y component

  // compute reflected particle frequency
  parfreq[p] *= exp( BETA*vi/(YP*af) );         // omegar
}
#endif

static void lstep(int stage, int nelem, int npar, int nwe, int nthreads,
                  int *inpoel, int *impcnt, int *elp, int *npeldt, int *wec,
                  int *esuel, int *npel, int *esupel1, int *esupel2, double dt,
                  double minx, double miny, double maxx, double maxy,
                  double *parcoord, double *parvel, double *dNx, double *dNy,
                  double *wrBA, double *wrA, double *wel, double *wenr,
                  double *du, double *Ae, double *pr, double *rg, double *coord,
                  int *wallstruck, int *psel1, int *psel2, double *u,
                  double *parfreq, double *f, double *u2, int *bruteforced
                  #ifndef WALLFUNCTIONS
                  , double *ddu, double *ru, double *avc0, int *wA, double *rho,
                  VSLStreamStatePtr *stream
                  #endif
                  #ifdef MICROMIXING
                  , double *parc, double *partm, double t
                    #ifdef VCIEM
                    , int *cp, double *vcte
                    #else
                    , double *c
                    #endif
                  #endif
                 )
// -----------------------------------------------------------------------------
// Routine: lstep - Advance particles
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int eN, j, jN, p, p2, p3, pNGRPP, A, B, C, NA, NB, NC, ws, bf, myid, UA, UB,
      UC;
  double a1, a2, wdist, l, eps, tke, dtmp, af, Pk, S;
  double dp[NDIMN], dv[3], dW[NGRPP], adu[NDIMN*NDIMN], rBA[NDIMN], nr[NDIMN],
         pc0[NDIMN], au[NDIMN], ge[9], au2[U2DOF];
  #ifndef WALLFUNCTIONS
  int RA, RB, RC;
  double c0, var, fmean;
  double wdist0, cn, redt, sredt;
  double addu[NDIMN];
  #endif

  // advance particles
  #ifdef ME
  dt *= ALPHA[stage];   // select stage if multi-stage timestepping
  #endif

  #ifndef WALLFUNCTIONS
  cn = 0.0;
  redt = 2.0/RE*dt;
  sredt = sqrt(redt);
  #endif
  ws = bf = 0;

  #ifdef _OPENMP
    #ifndef WALLFUNCTIONS
      #pragma omp parallel private(myid,p,p2,p3,pNGRPP,a1,a2,eN,A,B,C,NA,NB,NC,\
                                   UA,UB,UC,RA,RB,RC,au,adu,addu,af,tke,eps,   \
                                   au2,Pk,S,dp,ge,c0,dW,dv,dtmp,j,jN,wdist,    \
                                   wdist0,rBA,nr,l,fmean,var,pc0)              \
                                   reduction(+:cn,ws,bf)
    #else
      #pragma omp parallel private(myid,p,p2,p3,pNGRPP,a1,a2,eN,A,B,C,NA,NB,NC,\
                                   UA,UB,UC,au,adu,af,tke,eps,au2,Pk,S,dp,ge,  \
                                   dW,dv,dtmp,j,jN,wdist,rBA,nr,l,pc0)         \
                                   reduction(+:ws,bf)
    #endif // WALLFUNCTIONS
  #endif // _OPENMP
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<npar; p++)      // loop over all particles
    {
      p2 = p*2;
      p3 = p*3;
      pNGRPP = p*NGRPP;
      eN = elp[p]*NNODE;
      A = inpoel[eN+0];
      B = inpoel[eN+1];
      C = inpoel[eN+2];
      NA = A*NDIMN;
      NB = B*NDIMN;
      NC = C*NDIMN;
      UA = A*U2DOF;
      UB = B*U2DOF;
      UC = C*U2DOF;
      #ifndef WALLFUNCTIONS
      RA = A*RDOF;
      RB = B*RDOF;
      RC = C*RDOF;
      #endif

      // calculate element-average statistics
      au[0] = (u[NA+0]+u[NB+0]+u[NC+0])/3.0;    // <U>
      au[1] = (u[NA+1]+u[NB+1]+u[NC+1])/3.0;    // <V>
      adu[0] = (du[A*4+0]+du[B*4+0]+du[C*4+0])/3.0;     // d<U>/dx
      adu[1] = (du[A*4+1]+du[B*4+1]+du[C*4+1])/3.0;     // d<U>/dy
      adu[2] = (du[A*4+2]+du[B*4+2]+du[C*4+2])/3.0;     // d<V>/dx
      adu[3] = (du[A*4+3]+du[B*4+3]+du[C*4+3])/3.0;     // d<V>/dy
      #ifndef WALLFUNCTIONS
      // d^2<U>/dx^2 + d^2<U>/dy^2
      addu[0] = (ddu[NA+0]+ddu[NB+0]+ddu[NC+0])/3.0;
      // d^2<V>/dx^2 + d^2<V>/dy^2
      addu[1] = (ddu[NA+1]+ddu[NB+1]+ddu[NC+1])/3.0;
      #endif

      dp[0] = dNx[eN+0]*pr[A] + dNx[eN+1]*pr[B] + dNx[eN+2]*pr[C]; // d<p>/dx
      dp[1] = dNy[eN+0]*pr[A] + dNy[eN+1]*pr[B] + dNy[eN+2]*pr[C]; // d<p>/dy

      af = (f[A]+f[B]+f[C])/3.0;                  // mean turbulent frequency
      tke = ( (u2[UA+0]+u2[UB+0]+u2[UC+0])        // turbulent kinetic energy
             +(u2[UA+1]+u2[UB+1]+u2[UC+1])
             +(u2[UA+2]+u2[UB+2]+u2[UC+2]) )/6.0;
      if ( tke < EPS ) tke=BOUND;
      #ifndef WALLFUNCTIONS
      eps = af*(tke + CT2*af/RE);       // dissipation (elliptic relaxation)
      #else
      eps = af*tke;                     // dissipation (wall-functions)
      #endif
      if ( eps < EPS ) eps=BOUND;

      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (j=0; j<6; j++)       // Reynolds stress tensor
        au2[j] = (u2[UA+j]+u2[UB+j]+u2[UC+j])/3.0;
      // production of turbulent kinetic energy
      Pk = -au2[0]*adu[0] - au2[3]*(adu[1]+adu[2]) - au2[1]*adu[3];
      S = COM2 - COM1*Pk/eps;   // source/sink for turbulent frequency
      
      #ifndef WALLFUNCTIONS     // C0 and G using elliptic relaxation
        #ifdef __INTEL_COMPILER
        #pragma vector always
        #endif
        for ( j = 0; j < 9; j++ )
          // tensor rho (wiggly-p_ij)
          ge[j] = (rho[RA+j] + rho[RB+j] + rho[RC+j])/3.0;
        // c0 = -2rho_ij<uiuj>/(3k*eps)
        c0 = -2.0*(ge[0]*au2[0] +
                   ge[4]*au2[1] +
                   ge[8]*au2[2] +
                  (ge[1]+ge[3])*au2[3] +
                  (ge[2]+ge[6])*au2[4] +
                  (ge[5]+ge[7])*au2[5])/(3.0*tke*eps);
        // ge -= 0.5*eps*delta_ij
        ge[0] -= 0.5*eps;
        ge[4] -= 0.5*eps;
        ge[8] -= 0.5*eps;
        #ifdef __INTEL_COMPILER
        #pragma vector always
        #endif
        for (j=0; j<9; j++) ge[j] /= tke;       // ge /= tke
        if ( c0 < EPS ) c0 = BOUND;
        // add up c0 for diagnostics
        cn += c0;
      #else     // C0 and G using wall-functions with SLM
        ge[0] = ge[4] = ge[8] = -(0.5+0.75*C0)*eps/tke;
        ge[1] = ge[2] = ge[3] = ge[5] = ge[6] = ge[7] = 0.0;
      #endif

      // get NGRPP Gaussian random numbers into local array
      memcpy( dW, rg+pNGRPP, NGRPP*sizeof(double) );

      // compute velocity increment
      #ifndef WALLFUNCTIONS     // with elliptic relaxation
      dtmp = sqrt(c0*eps*dt);
      dv[0] = -dp[0]*dt +
              redt*addu[0] +
              sredt*(adu[0]*dW[0] + adu[1]*dW[1])*stage +
              (ge[0]*(parvel[p3+0]-au[0]) +
               ge[1]*(parvel[p3+1]-au[1]) +
               ge[2]*parvel[p3+2])*dt +
              dtmp*dW[2]*stage;
      dv[1] = -dp[1]*dt +
              redt*addu[1] +
              sredt*(adu[2]*dW[0] +
              adu[3]*dW[1])*stage +
              (ge[3]*(parvel[p3+0]-au[0]) +
               ge[4]*(parvel[p3+1]-au[1]) +
               ge[5]*parvel[p3+2])*dt +
              dtmp*dW[3]*stage;
      dv[2] = 0.0 +
              (ge[6]*(parvel[p3+0]-au[0]) +
               ge[7]*(parvel[p3+1]-au[1]) +
               ge[8]*parvel[p3+2])*dt +
              dtmp*dW[4]*stage;
      #else     // with wall-functions
      dtmp = sqrt(C0*eps*dt);
      dv[0] = -dp[0]*dt +
              (ge[0]*(parvel[p3+0]-au[0]) +
               ge[1]*(parvel[p3+1]-au[1]) +
               ge[2]*parvel[p3+2])*dt +
              dtmp*dW[0]*stage;
      dv[1] = -dp[1]*dt +
              (ge[3]*(parvel[p3+0]-au[0]) +
               ge[4]*(parvel[p3+1]-au[1]) +
               ge[5]*parvel[p3+2])*dt +
              dtmp*dW[1]*stage;
      dv[2] = (ge[6]*(parvel[p3+0]-au[0]) +
               ge[7]*(parvel[p3+1]-au[1]) +
               ge[8]*parvel[p3+2])*dt +
              dtmp*dW[2]*stage;
      #endif
      // update particle velocity
      parvel[p3+0] += dv[0];
      parvel[p3+1] += dv[1];
      parvel[p3+2] += dv[2];

      // update particle frequency
      #ifndef WALLFUNCTIONS     // with elliptic relaxation
      parfreq[p] += (-C3*af*(parfreq[p]-af) - S*af*parfreq[p])*dt +
                    af*sqrt(2.0*C3*C4*parfreq[p]*dt)*dW[5]*stage;
      #else     // with wall-functions
      parfreq[p] += (-C3*af*(parfreq[p]-af) - S*af*parfreq[p])*dt +
                    af*sqrt(2.0*C3*C4*parfreq[p]*dt)*dW[3]*stage;
      #endif
      if (parfreq[p]<0.0) parfreq[p]=0.0;

      // save (x,y of) particle position from previous timestep
      pc0[0] = parcoord[p2+0];
      pc0[1] = parcoord[p2+1];

      // update particle position
      #ifndef WALLFUNCTIONS
      parcoord[p2+0] += parvel[p3+0]*dt + sredt*dW[0]*stage;
      parcoord[p2+1] += parvel[p3+1]*dt + sredt*dW[1]*stage;
      #else
      parcoord[p2+0] += parvel[p3+0]*dt;
      parcoord[p2+1] += parvel[p3+1]*dt;
      #endif

      // enforce wall boundary condition on particle
      for (j=0; j<nwe; j++)     // loop over all wall-elements
        if (j == wec[elp[p]]) { // only bother with the closest wall element
        // (this is not so bulletproof, but faster than checking all)
          jN = j*NDIMN;
          // compute distance between new particle position and wall-element j
          wdist = (wrBA[jN+0]*(parcoord[p2+1]-wrA[jN+1]) -
                   wrBA[jN+1]*(parcoord[p2+0]-wrA[jN+0]))/wel[j];

          #ifndef WALLFUNCTIONS
          // compute distance between old particle position and wall-element j
          wdist0 = (wrBA[jN+0]*(pc0[1]-wrA[jN+1]) -
                    wrBA[jN+1]*(pc0[0]-wrA[jN+0]))/wel[j];
          #endif
          
          // test particle wall-trajectories
          if (wdist < 0.0) {    // trajectory 1
            rBA[0] = wrBA[jN+0]/wel[j]; // calculate normalized edge-vector rBA
            rBA[1] = wrBA[jN+1]/wel[j];
            nr[0] = wenr[jN+0];         // get normal of wall-element j
            nr[1] = wenr[jN+1];
            l = 2.0*(rBA[1]*(wrA[jN+0]-parcoord[p2+0]) +
                     rBA[0]*(parcoord[p2+1]-wrA[jN+1])) /
                    (rBA[1]*nr[0] - rBA[0]*nr[1]);
            parcoord[p2+0] += l*nr[0];  // reflect particle
            parcoord[p2+1] += l*nr[1];

            #ifndef WALLFUNCTIONS
            // no slip
            parvel[p3+0] = parvel[p3+1] = parvel[p3+2] = 0.0;
            // sample gamma distribution with mean w1 and C4*w1*w1 variance
            // w1 = f[node_A_of_wall-element], calculated in stat()
            fmean = f[wA[j]];
            var = C4*fmean*fmean;
            CheckVslError( vdRngGamma(GAMMA_METHOD, stream[myid], 1, parfreq+p,
                                      fmean*fmean/var, 0.0, var/fmean) );
            #else
            // specify particle conditions based on wall-functions
            parwallfunctionbc(p, tke, af, au, au2, dp, rBA, nr, parvel,
                              parfreq);
            #endif

            #ifdef MICROMIXING            
            parc[p] = SS;       // heated cylinder
            #endif

            ws++;       // increase number of particles that struck the wall
          }
          #ifndef WALLFUNCTIONS
          else
            if (exp(-wdist0*wdist*RE/dt) > ru[p]) { // trajectory 3
              parvel[p3+0] = parvel[p3+1] = parvel[p3+2] = 0.0; // no slip
              // sample gamma distribution with mean w1 and C4*w1*w1 variance
              // w1 = f[node_A_of_wall-element], calculated in stat()
              fmean = f[wA[j]];
              var = C4*fmean*fmean;
              CheckVslError( vdRngGamma(GAMMA_METHOD, stream[myid], 1,
                                        parfreq+p, fmean*fmean/var, 0.0,
                                        var/fmean) );

              #ifdef MICROMIXING
              parc[p] = SS;     // heated cylinder
              #endif

              ws++;     // increase number of particles that struck the wall
            }
          #endif
        }

      // put back particle if it tries to escape from the domain through illegal
      // means (this is a sanity check (and fix) in case in the above fancy
      // wall-condition something went wrong, or rather the particle nastily
      // sneaked through between two wall elements)
      if ((parcoord[p2+0]*parcoord[p2+0]+parcoord[p2+1]*parcoord[p2+1]) < D*D/4)
      { // put it back where it was in the previous timestep
        parcoord[p2+0] = pc0[0];
        parcoord[p2+1] = pc0[1];
      }

      // enforce top boundary condition on particle (free-slip)
      if (parcoord[p2+1] > maxy) {
         parcoord[p2+1] = 2.0*maxy-parcoord[p2+1]; // put it back
         parvel[p3+1] = -parvel[p3+1];             // with opposite y velocity

         #ifdef MICROMIXING
         parc[p] = 0.0;         // scalar goes out
         #endif
      }
      // enforce bottom boundary condition on particle (free-slip)
      if (parcoord[p2+1] < miny) {
         parcoord[p2+1] = 2.0*miny-parcoord[p2+1]; // put it back
         parvel[p3+1] = -parvel[p3+1];             // with opposite y velocity

         #ifdef MICROMIXING
         parc[p] = 0.0;         // scalar goes out
         #endif
      }
      // enforce inflow/outflow boundary condition on particle
      if (parcoord[p2+0] > maxx) {      // outflow
        parcoord[p2+0] = minx+0.0001;   // put it back at the inflow

        parvel[p3+0] = 1.0;             // prescribed velocity at inflow
        parvel[p3+1] = parvel[p3+2] = 0.0;

        #ifdef MICROMIXING
        parc[p] = 0.0;  // no scalar coming in at inflow
        #endif
      }
      if ( parcoord[p2+0] < minx ) parcoord[p2+0] = maxx-0.0001;

      // find particle
      if (!PARINEL(p,elp[p],a1,a2,coord,inpoel,parcoord,Ae))
        findpar(p, elp[p], myid, nelem, minx, maxx, npeldt, elp, inpoel, esuel,
                esupel1, esupel2, coord, parcoord, Ae, &bf);

      // micromixing
      #ifdef MICROMIXING
      if (t > RETIME) {  // release scalar at RETIME
        // micromixing timescale (tm)
        #ifndef WALLFUNCTIONS
        //partm[p] = fmin(Cs*pow(DSHS/eps,1.0/3.0) + Ct*tke/eps,
        //                fmax(tke/eps,CT/sqrt(RE*eps)));
        partm[p] = fmax( Ct*tke/eps, 0.01 );
        #else
        partm[p] = fmin( Cs*pow(DSHS/eps,1.0/3.0) + Ct*tke/eps, tke/eps );
        #endif

        // scalar increment
        #ifdef VCIEM
        parc[p] += (vcte[elp[p]*CNBI+cp[p]]-parc[p])*dt/partm[p];  // VCIEM
        #else
        parc[p] += ((c[A]+c[B]+c[C])/3.0 - parc[p])*dt/partm[p];   // IEM
        #endif
      }
      #endif
    }
  }
  #ifndef WALLFUNCTIONS
  (*avc0) = cn/npar;    // average c0 over npar particles (for diagnostics)
  #endif
  // number of particles that struck the wall in this timestep (for diagnostics)
  (*wallstruck) = ws;
  // number of particles that had to be brute-forced in this timestep (for
  // diagnostics)
  (*bruteforced) = bf;

  // update psel, npel, elp
  genpsel( nelem, npar, nthreads, npeldt, npel, elp, psel1, psel2 );

  // improve spatial distribution of particles if needed
  improvepar2(nelem, npar, nthreads, impcnt, npeldt, elp, npel, parcoord,
              parvel, psel1, psel2, parfreq
              #ifdef MICROMIXING
              , parc
              #endif
             );

  // regenerate random number tables
  regenrng_tables( nthreads,
                   #ifndef WALLFUNCTIONS
                   ru,
                   #endif
                   rg );
}

static double calcdt(int npoin, double *sqrtAp, double *u, double *f
                     #ifdef MICROMIXING
                     , double *c
                     #endif
                    )
// -----------------------------------------------------------------------------
// Routine: calcdt -  Calculate the size of the next timestep
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int p, pN;
  double vel, dt, mindt;


  // find maximum nodal velocity
  mindt = 1.0;
  #ifdef _OPENMP
  #pragma omp parallel for private(p,pN,vel,dt)
  #endif
  for ( p = 0; p < npoin; p++ ) {
    pN = p*NDIMN;
    vel = sqrt( u[pN+0]*u[pN+0] + u[pN+1]*u[pN+1] + f[p]*f[p]/RE/RE
                #ifdef MICROMIXING
                + c[p]*c[p]
                #endif
              );

    dt = sqrtAp[p]/vel;
    if ( dt < mindt )
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      {
        if ( dt < mindt )
          mindt = dt;
      }
    
  }

  return( CFL*mindt );
}

void langevin(int nelem, int npoin, int nbpoin, int it0, int it, int nwe,
              int npar, int nthreads, int restarted, int maxn, int ec_size,
              double minx, double miny, double maxx, double maxy, double t,
              double *drag, double *lift, double *tdrag, double *tlift,
              int *wec, int *esup1, int *esup2, int *bpg, int *inpoel,
              int *binpoel, int *npel, int *elp, int *npeldt, int *we, int *weo,
              int *esuel, int *betags, int *esupel1, int *esupel2, int *onum,
              int *az1, int *az2, double *r_az, double *dir_az, double *N_az,
              double *odpn, double *oudt, double *wenr, double *sqrtAp,
              double *wrBA, double *wrA, double *wel, double *prhs, double *dpr,
              double *Ae, double *dNx, double *dNy, double *parcoord,
              double *parvel, double *u, double *ue, double *tu, double *tpr,
              double *du, double *ddu, double *coord, double *pr, double *rg,
              sparsemat *P, int *psel1, int *psel2, int ndl, int *es1, int *es2,
              int *ec, double *uprof, double *vprof, double *u2prof,
              double *v2prof, double *u3prof, double *v3prof, double *u4prof,
              double *v4prof, double *uvprof, double *tkeprof, double *epsprof,
              double *dete, double *yprof, double *xprof_c, double *uprof_c,
              double *pprof_c, int angle_size, int *inprof_w, double *angle_w,
              double *parfreq, double *u2, double *u3, double *u4, double *f,
              double *u2e, double *u3e, double *u4e, double *fe, double *tu2,
              double *tu3, double *tu4, double *tf, double *tdu
              #ifndef WALLFUNCTIONS
              , int *bpsup1, int *bpsup2, int *bptags, double *ru,
              VSLStreamStatePtr *stream, int *wA, int *wlc, double *wpnr,
              double *wnz, double *rrhs, sparsemat *R, double *rho
              #endif
              #ifdef MICROMIXING
              , double *parc, double *partm, double *c, double *c2, double *c3,
              double *c4, double *tm, double *uc, double *ce, double *c2e,
              double *c3e, double *c4e, double *tme, double *uce, double *tc,
              double *tc2, double *tc3, double *tc4, double *ttm, double *tuc,
              int npl, int *epdfloc, int *sl, double *tpdf, double *cprof,
              double *c2prof, double *c3prof, double *c4prof, double *tmprof
                #ifdef VCIEM
                 , int *cp, double *vcte
                 #endif
               #endif
             )
// -----------------------------------------------------------------------------
// Routine: langevin - Main timestepping loop using the Langevin equation
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int wallstruck=0, pit=0, minnp=0, sm=10000, ne=0, lne=0, nc=0, impcnt, msc=0,
       nic=0, bf, stage;
  long int hrs2end=0, mins2end=0, secs2end=0, hrs2beg=0, mins2beg=0, secs2beg=0;
  double dt;
  FILE *opos;
  struct timeval start_time;

  #ifndef WALLFUNCTIONS
  int rit=0;
  double avc0=0.0;
  #endif

  // clear contents of temporal results file if not restarted
  if (!restarted) {
    if ( !(opos = fopen(OUT_POSTPROCESS_TIME_FILENAME,"w")) )
      ERR("Cannot open OUT_POSTPROCESS_TIME_FILENAME\n");
    fclose( opos );
  }

  // output particle positions
  //outparpos( npar, t, parcoord );

  gettimeofday( &start_time, (struct timezone*)0 );
  // main timsetepping loop...
  for (; (t<MAXTIME) && (it<MAXTS); t+=dt,it++) {
    // calculate size of timestep
    dt = calcdt( npoin, sqrtAp, u, f
                 #ifdef MICROMIXING
                 , c
                 #endif
               );

    #ifdef ME
    for (stage=0; stage<STAGES; stage++) {  // multi-stage timestepping scheme
    #else
    stage = 1;  // Euler-Maruyama timestepping scheme
    #endif

      #ifndef WALLFUNCTIONS
      // elliptic relaxation
      rit = rstep(nelem, nbpoin, nthreads, bpg, inpoel, wlc, bpsup1, bpsup2,
                  bptags, Ae, wpnr, rrhs, rho, dNx, dNy, u2, f, wnz, du, R);
      #endif

      // advance particles
      lstep(stage, nelem, npar, nwe, nthreads, inpoel, &impcnt, elp, npeldt,
            wec, esuel, npel, esupel1, esupel2, dt, minx, miny, maxx, maxy,
            parcoord, parvel, dNx, dNy, wrBA, wrA, wel, wenr, du, Ae, pr, rg,
            coord, &wallstruck, psel1, psel2, u, parfreq, f, u2, &bf
            #ifndef WALLFUNCTIONS
            , ddu, ru, &avc0, wA, rho, stream
            #endif
            #ifdef MICROMIXING
            , parc, partm, t
              #ifdef VCIEM
              , cp, vcte
              #else
              , c
              #endif
            #endif
           );

      #ifdef ME
      if (stage == STAGES-1) {  // pressure projection only in the last stage
      #endif
        // pressure projection for mean pressure
        pit = pstep(nwe, npoin, nelem, nbpoin, nthreads, dt, maxx, inpoel,
                    binpoel, we, weo, bpg, betags, odpn, oudt, wel, wenr, prhs,
                    du, ddu, pr, dpr, u2, dNx, dNy, Ae, coord, P, u);
        // correct velocities after pressure projection
        velcorr(npar, dt, inpoel, elp, parvel, dpr, dNx, dNy);
      #ifdef ME
      }
      #endif

      // extract statistics from particles
      ne = stat(nelem, npoin, nbpoin, npar, nthreads, dt, npel, esup1, esup2,
                inpoel, elp, nwe, wenr, ue, u, du, ddu, parvel, dNx, dNy, &nc,
                betags, binpoel, oudt, u2e, u3e, u4e, fe, u2, u3, u4, f,
                parfreq
                #ifndef WALLFUNCTIONS
                , bptags, bpg, we, weo
                #endif
                #ifdef MICROMIXING
                , ce, c2e, c3e, c4e, tme, uce, c, c2, c3, c4, tm, uc, parc,
                partm
                  #ifdef VCIEM
                  , cp, vcte, psel1, psel2
                  #endif
                #endif
               );

      #ifdef ME
      if (stage == STAGES-1) {  // pressure projection only in the last stage
      #endif
        // save normal derivative of mean pressure at wall-elements for next
        // timestep
        save_old_dpn(nwe, nbpoin, binpoel, betags, wenr, ddu, odpn, oudt, we,
                     inpoel, u2, dNx, dNy, u, du);
      #ifdef ME
      }
      #endif

    #ifdef ME
    } // end of for-loop for multi-stage timestepping scheme
    #endif

    if ((minnp=minnpel(nelem,npel)) < sm) sm = minnp;   // find smallest minnpel
    if (impcnt > msc) msc = impcnt;               // find largest improve count
    if (impcnt != 0) nic++;     // count number of iterations needed correction
    if (ne > lne) lne = ne;     // find largest number of empty element
    
    // dump results and restartfile at every ST2DUMPth iteration
    if (!(it%ST2DUMP)) {
       // output instantaneous 3d fields
       out_3dstep( nelem, npoin, inpoel, coord, u, pr, u2, u3, u4, f, du,
                   /*#ifndef WALLFUNCTIONS
                    rho,
                   #endif*/
                   #ifdef MICROMIXING
                   c, c2, c3, c4, tm,
                   #endif
                   OUT_POSTPROCESS_3D_FILENAME, INST );

       // output time-averaged fields
       if (t > AVTIME) {
         out_3dstep( nelem, npoin, inpoel, coord, tu, tpr, tu2, tu3, tu4, tf,
                     tdu,
                     /*#ifndef WALLFUNCTIONS
                     rho,
                     #endif*/
                     #ifdef MICROMIXING
                     tc, tc2, tc3, tc4, ttm,
                     #endif
                     OUT_POSTPROCESS_3DTAV_FILENAME, TAV );
       }

       // output particle positions
       //outparpos(npar, t, parcoord);

       // output restartfile
       #ifdef SAVERESTART
       saverestartpar(nthreads, npoin, npar, *onum, it, t, parcoord, parvel, tu,
                      tdu, parfreq, tu2, tu3, tu4, tf, tpr
                      #ifndef WALLFUNCTIONS
                      , stream
                      #endif
                      #ifdef MICROMIXING
                      , parc, tc, tc2, tc3, tc4, tuc
                      #endif
                     );
       #endif
    }

    //outparpos(npar, t, parcoord);

    // output 'instantaneous' fields, spatially-averaged fields and
    // compute time-averaged fields (if time-averaging started)
    out_inst_dist(npoin, nbpoin, minnp, pit, maxn, t, dt, drag, lift, tdrag,
                  tlift, betags, binpoel, onum, u, wenr, wel, du, ndl, maxx,
                  maxy, es1, es2, ec_size, ec, az1, az2, xprof_c, uprof_c,
                  pprof_c, inpoel, coord, uprof, vprof, u2prof, v2prof, u3prof,
                  v3prof, u4prof, v4prof, uvprof, tkeprof, epsprof, dete, yprof,
                  angle_size, inprof_w, angle_w, r_az, dir_az, N_az,
                  #ifndef WALLFUNCTIONS
                  rit,
                  #endif
                  #ifdef MICROMIXING
                  c, c2, c3, c4, tm, uc, tc, tc2, tc3, tc4, ttm, tuc, npl,
                  epdfloc, esupel1, esupel2, c2e, psel1, psel2, parc, ce, sl,
                  tpdf, cprof, c2prof, c3prof, c4prof, tmprof,
                  #endif
                  pr, tu, u2, u3, u4, f, tu2, tu3, tu4, tf, tpr, tdu );

    // calculate estimated time left
    calctimes(&start_time, it, it0, t, dt,
              &hrs2beg, &mins2beg, &secs2beg, &hrs2end, &mins2end, &secs2end );

    // output timestep info
    printf("it = %d, t = %.5g\tdt = %g"
           #ifdef VCIEM
           "  c:%.3g"
           #endif
           " ["
           #ifndef WALLFUNCTIONS
           "R:%d "
           #endif
           "P:%d] %d,%d,%d"
           " b:%d w:%d"
           #ifndef WALLFUNCTIONS
           " c0 = %g"
           #endif
           "  %ld:%ld:%ld\t%ld:%ld:%ld\n",
           it, t, dt,
           #ifdef VCIEM
           (double)nc/nelem,
           #endif
           #ifndef WALLFUNCTIONS
           rit,
           #endif
           pit, minnp, impcnt, ne, bf, wallstruck,
           #ifndef WALLFUNCTIONS
           avc0,
           #endif
           hrs2beg, mins2beg, secs2beg, hrs2end, mins2end, secs2end);
    fflush(stdout);
  }
  printf("-----------------------------------------------------------------\n");
  printf("Run stats:\n");
  printf(" * smallest minnpel: %d\n",sm);
  printf(" * most severe correction: %d\n",msc);
  printf(" * number of iterations corrected: %d (%.3g%%)\n",nic,100.0*nic/it);
  printf(" * largest number of empty elements: %d (%.3g%%)\n",lne,
                                                              100.0*lne/nelem);
  fflush(stdout);

  // output time-averaged statistics on the cylinder surface
  out_surface_stats(tdrag, tlift);

  // output instantaneous 3d fields
  out_3dstep(nelem, npoin, inpoel, coord, u, pr, u2, u3, u4, f, du,
             /*#ifndef WALLFUNCTIONS
             rho,
             #endif*/
             #ifdef MICROMIXING
             c, c2, c3, c4, tm,
             #endif
             OUT_POSTPROCESS_3D_FILENAME, INST);

  // output time-averaged 3d fields
  out_3dstep(nelem, npoin, inpoel, coord, tu, tpr, tu2, tu3, tu4, tf, tdu,
             /*#ifndef WALLFUNCTIONS
             rho,
             #endif*/
             #ifdef MICROMIXING
             tc, tc2, tc3, tc4, ttm,
             #endif
             OUT_POSTPROCESS_3DTAV_FILENAME, TAV );

  // output restartfile
  #ifdef SAVERESTART
  saverestartpar(nthreads, npoin, npar, *onum, it, t, parcoord, parvel, tu, tdu,
                 parfreq, tu2, tu3, tu4, tf, tpr
                 #ifndef WALLFUNCTIONS
                 , stream
                 #endif
                 #ifdef MICROMIXING
                 , parc, tc, tc2, tc3, tc4, tuc
                 #endif
                );
  #endif
}
