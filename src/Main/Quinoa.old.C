// -----------------------------------------------------------------------------
// \file    src/Main/Quinoa.C
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Quinoa main
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <cstdio>
#include <cmath>
#ifdef _OPENMP
#include "omp.h"
#endif
#include "mkl.h"
#include "Macros.h"
#include "Const.h"
#include "Mesh.h"
#include "SparseMatrix.h"
#include "Random.h"
#ifndef WALLFUNCTIONS
#include "EllipticRelaxation.h"
#endif
#include "Pressure.h"
#include "Langevin.h"
#include "Particles.h"
#include "Postprocess.h"

// global, static variables and pointers to dynamically allocated arrays
static int nelem, npoin, nbpoin, npar, nifp, nofp, nwnzr, nwnzc, nwe, it, it0,
           nthreads, restarted, onum, ndl, maxn, ec_size, angle_size;
static int *inpoel, *binpoel, *bpg, *esup1, *esup2, *psup1, *psel1, *psel2,
           *es1, *es2, *ec, *psup2, *bpsup1, *bpsup2, *elp, *npel, *esupel1,
           *esupel2, *esuel, *ifl, *ofl, *we, *weo, *wA, *wec, *npeldt, *bptags,
           *betags, *wlc, *wnzr, *inprof_w, *az1, *az2;
static double minsqrtAp, t, t0, minx, miny, maxx, maxy;
static double drag[2], lift[2], tdrag[2], tlift[2], dir_az[MAXAZ*NDIMN];
static double *Ae, *sqrtAp, *dNx, *dNy, *dete, *coord, *parcoord, *parvel,
              *ifpp, *ofpp, *ifnz, *ofnz, *wpnr, *wenr, *wel, *wrBA, *wrA,
              *odpn, *tu, *tdu, *u, *du, *ddu, *ue, *pr, *prhs, *dpr, *tpr, *rg,
              *oudt, *uprof, *vprof, *u2prof, *v2prof, *u3prof, *v3prof,
              *u4prof, *v4prof, *uvprof, *tkeprof, *epsprof, *yprof, *xprof_c,
              *uprof_c, *pprof_c, *angle_w, *r_az, *N_az, *parfreq, *u2, *u3,
              *u4, *f, *fe, *u2e, *u3e, *u4e, *tu2, *tu3, *tu4, *tf;
static sparsemat P;
static VSLStreamStatePtr *stream;

#ifdef MICROMIXING
static int npl;
static int *sl, *epdfloc;
static double *parc, *partm, *tc, *tc2, *tc3, *tc4, *tuc, *c, *c2, *c3, *c4,
              *ce, *c2e, *c3e, *c4e, *uc, *uce, *tm, *tme, *ttm, *cprof,
              *c2prof, *c3prof, *c4prof, *tmprof, *tpdf;
#ifdef VCIEM
static int *cp;
static double *vcte;
#endif
#endif

#ifndef WALLFUNCTIONS
static double *ru, *rrhs, *wnz, *rho;
static sparsemat R;
#endif

static void createvecmat()
// -----------------------------------------------------------------------------
// Routine: creatyevecemat - Create vectors and matrices that store physical
//                           variables and coefficient matrices, all these are
//                           destroyed in destroyvecmat()
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  ////////////////////////////
  // POINT-BASED QUANTITIES
  #ifndef WALLFUNCTIONS
  // elliptic relaxation tensor and and rhs of elliptic relaxation solve
  if ( !(rho = (double*)calloc(npoin*RDOF,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(rrhs = (double*)calloc(npoin*RDOF,sizeof(double))) )
    ERR("can't allocate memory!");
  #endif
  // mean velocity and its second derivatives
  if ( !(u = (double*)calloc(npoin*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(ddu = (double*)calloc(npoin*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  // mean turbulent frequency, mean pressure, rhs for mean pressure,
  // difference of mean pressure, unconditioned mean, variance, skewness,
  // flatness of scalar, micromixing timescale
  if ( !(f = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(pr = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(prhs = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(dpr = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  #ifdef MICROMIXING
  if ( !(c = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(c2 = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(c3 = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(c4 = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(tm = (double*)calloc(npoin,sizeof(double))) )
    ERR("can't allocate memory!");
  #endif
  // first derivatives of mean velocity
  if ( !(du = (double*)calloc(npoin*NDIMN*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  // Reynolds stress tensor
  if ( !(u2 = (double*)calloc(npoin*U2DOF,sizeof(double))) )
    ERR("can't allocate memory!");
  // skewness and flatness of streamwise and cross-stream velocity
  if ( !(u3 = (double*)calloc(npoin*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(u4 = (double*)calloc(npoin*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  // velocity-scalar correlation vector
  #ifdef MICROMIXING
  if ( !(uc = (double*)calloc(npoin*3,sizeof(double))) )
    ERR("can't allocate memory!");
  #endif

  //////////////////////////////
  // ELEMENT-BASED QUANTITIES
  // mean velocity
  if ( !(ue = (double*)calloc(nelem*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  // mean turbulent frequency, unconditioned mean, variance, skewness,
  // flatness of scalar, micromixing timescale
  if ( !(fe = (double*)calloc(nelem,sizeof(double))) )
    ERR("can't allocate memory!");
  #ifdef MICROMIXING
  if ( !(ce = (double*)calloc(nelem,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(c2e = (double*)calloc(nelem,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(c3e = (double*)calloc(nelem,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(c4e = (double*)calloc(nelem,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(tme = (double*)calloc(nelem,sizeof(double))) )
    ERR("can't allocate memory!");
  #endif
  // Reynolds stress tensor
  if ( !(u2e = (double*)calloc(nelem*U2DOF,sizeof(double))) )
    ERR("can't allocate memory!");
  // skewness and flatness of streamwise and cross-stream velocity
  if ( !(u3e = (double*)calloc(nelem*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(u4e = (double*)calloc(nelem*NDIMN,sizeof(double))) )
    ERR("can't allocate memory!");
  #ifdef MICROMIXING
  #ifdef VCIEM
  // velocity-conditioned mean of scalar
  if ( !(vcte = (double*)calloc(nelem*CNBI,sizeof(double))) )
    ERR("can't allocate memory!");
  #endif
  // velocity-scalar correlation vector
  if ( !(uce = (double*)calloc(nelem*3,sizeof(double))) )
    ERR("can't allocate memory!");
  #endif

  ///////////////////
  // sparse matrices
  // for elliptic relaxation
  #ifndef WALLFUNCTIONS
  create_sparsemat( &R, RDOF, npoin, psup1, psup2 );
  #endif
  // for mean pressure solve
  create_sparsemat( &P, 1, npoin, psup1, psup2 );
}

static void Solve( void )
// -----------------------------------------------------------------------------
// Routine: Solve - Initialize linear solver for mean pressure, print out some
//                  info, call main timestepping loop
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  // initialize linear solver for mean pressure solve
  pinit(nelem, nbpoin, nthreads, inpoel, bpg, bpsup1, bpsup2,
        maxx, coord, Ae, dNx, dNy, ofnz, &P);
  
  printf("------------------------------------------------------------------\n");
  printf(" * times: t0 = %.3g   "
         #ifdef MICROMIXING
         "RETIME = %.3g   "
         #endif
         "AVTIME = %.3g   MAXTIME = %.3g\n"
         " * timesteps: it0 = %d   MAXTS = %d\n"
         " * iterative solver: STOP_TOL = %.3g   MAX_IT = %d\n"
         " * CFL = %.3g   Cp = %.3g   NBI = %d   ST2DUMP = %d\n",
         t0,
         #ifdef MICROMIXING
         RETIME,
         #endif
         AVTIME, MAXTIME,
         it0, MAXTS, STOP_TOL, MAX_IT, CFL, Cp, NBI, ST2DUMP );

  printf(" * solver: turbulent (");
  #ifdef WALLFUNCTIONS
  printf("wall-functions)\n");
  #else
  printf("elliptic relaxation)\n");
  #endif // WALLFUNCTIONS

  #ifdef ME
  printf(" * timestepping scheme: 2-stage modified Euler\n");
  #else
  printf(" * timestepping scheme: Euler-Maruyama\n");
  #endif

  #ifdef MICROMIXING
    #ifdef VCIEM
      #ifdef PROJECTION
      printf(" * micromixing: VCIEM, CNBI = %d, projection\n", CNBI);
      #else
      printf(" * micromixing: VCIEM, CNBI = %d, no projection\n", CNBI);
      #endif
    #else
    printf(" * micromixing: IEM\n" );
    #endif
  #else
  printf(" * no micromixing\n" );
  #endif

  printf(" * results: DL = %d, AZ = %d", ndl, MAXAZ);
  #ifdef MICROMIXING
  printf(", PL = %d\n", npl);
  #else
  printf("\n");
  #endif
  printf("----------------------------------------------------------------\n" );
  fflush(stdout);

  // main timestepping loop...
  langevin(nelem, npoin, nbpoin, it0, it, nwe, npar, nthreads, restarted, maxn,
           ec_size, minx, miny, maxx, maxy, t, drag, lift, tdrag, tlift, wec,
           esup1, esup2, bpg, inpoel, binpoel, npel, elp, npeldt, we, weo,
           esuel, betags, esupel1, esupel2, &onum, az1, az2, r_az, dir_az, N_az,
           odpn, oudt, wenr, sqrtAp, wrBA, wrA, wel, prhs, dpr, Ae, dNx, dNy,
           parcoord, parvel, u, ue, tu, tpr, du, ddu, coord, pr, rg, &P, psel1,
           psel2, ndl, es1, es2, ec, uprof, vprof, u2prof, v2prof, u3prof,
           v3prof, u4prof, v4prof, uvprof, tkeprof, epsprof, dete, yprof,
           xprof_c, uprof_c, pprof_c, angle_size, inprof_w, angle_w, parfreq,
           u2, u3, u4, f, u2e, u3e, u4e, fe, tu2, tu3, tu4, tf, tdu
           #ifndef WALLFUNCTIONS
           , bpsup1, bpsup2, bptags, ru, stream, wA, wlc, wpnr, wnz, rrhs, &R,
           rho
           #endif
           #ifdef MICROMIXING
           , parc, partm, c, c2, c3, c4, tm, uc, ce, c2e, c3e, c4e, tme, uce,
           tc, tc2, tc3, tc4, ttm, tuc, npl, epdfloc, sl, tpdf, cprof, c2prof,
           c3prof, c4prof, tmprof
             #ifdef VCIEM
             , cp, vcte
             #endif
           #endif
          );
}

static void calcbounds( void )
// -----------------------------------------------------------------------------
// Routine: caclbounds - Calculate bounds of domain in x and y directions
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i, iN;

  minx = miny = +1.0e+10;
  maxx = maxy = -1.0e+10;
  for ( i = 0; i < npoin; i++ ) {
    iN = i*NDIMN;
    if (coord[iN+0] < minx) minx = coord[iN+0];
    if (coord[iN+0] > maxx) maxx = coord[iN+0];
    if (coord[iN+1] < miny) miny = coord[iN+1];
    if (coord[iN+1] > maxy) maxy = coord[iN+1];
  }
}

static void create_bc_structs( void )
// -----------------------------------------------------------------------------
// Routine: create_bc_structs - Create derived data structures for setting
//                              boundary conditions
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i, j, jN, A, B, cif, cof, k, haveit;
  double dtmp;
  double nr[NNODE];
  int cw;

  ////////////////
  // auxiliary arrays for setting dynamically changing inhomogeneous Dirichlet
  // conditions at inflow/outflow for the pressure Poisson equation
  //
  // arrays that store the global-to-local mapping of boundary node indices
  // along the inflow and outflow
  if ( !(ifl = (int*)calloc(npoin,sizeof(int))) )
    ERR("Can't allocate memory!");
  if ( !(ofl = (int*)calloc(npoin,sizeof(int))) )
    ERR("Can't allocate memory!");
  // count up number of boundary points at inflow and outflow and generate
  // global-to-local mapping arrays
  for (nifp=nofp=i=0; i < nbpoin; i++) {
     if (fabs(coord[bpg[i]*NDIMN+0]-minx) < EPS)
       ifl[bpg[i]] = nifp++;    // inflow
     if (fabs(coord[bpg[i]*NDIMN+0]-maxx) < EPS)
       ofl[bpg[i]] = nofp++;    // outflow
  }
  // allocate arrays that store values of Dirichlet conditions for
  // inflow/outflow
  if ( !(ifpp = (double*)calloc(nifp,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(ofpp = (double*)calloc(nofp,sizeof(double))) )
    ERR("Can't allocate memory!");
  // count up total number of nonzero matrix elements in boundary columns of
  // matrix P along inflow/outflow
  for (cif=cof=i=0; i<nbpoin; i++) {    // loop over all boundary points
     // add number of points surrounding inflow boundary point i (+ main diag)
     if (fabs(coord[bpg[i]*NDIMN+0]-minx) < EPS)
       cif += bpsup2[i+1]-bpsup2[i]+1;
     // add number of points surrounding outflow boundary point i (+ main diag)
     if (fabs(coord[bpg[i]*NDIMN+0]-maxx) < EPS)
       cof += bpsup2[i+1]-bpsup2[i]+1;
  }
  // allocate arrays to store all nonzero values of columns for inflow/outflow
  // points of matrix P
  if ( !(ifnz = (double*)calloc(cif,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(ofnz = (double*)calloc(cof,sizeof(double))) )
    ERR("Can't allocate memory!");

  ////////////////
  // auxiliary arrays for setting inhomogeneous Dirichlet conditions
  // at the wall for the elliptic relaxation equation
  // count up total number of nonzero matrix elements in boundary columns of
  // matrix R along the wall
  for ( cw=i=0; i<nbpoin; i++) {        // loop over all boundary points
     // add number of points surrounding wall boundary point i (+ main diag)
     if (bptags[i*3+1] == 4)    // wall
       cw += bpsup2[i+1]-bpsup2[i]+1;
  }
  #ifndef WALLFUNCTIONS
  // allocate array to store all nonzero values and row indices of columns for
  // wall points of matrix R
  if ( !(wnz = (double*)calloc(cw*RDOF,sizeof(double))) )
    ERR("Can't allocate memory!");
  #endif

  ////////////////
  // auxiliary arrays for setting inhomogeneous Dirichlet conditions at the wall
  // for elliptic relaxation and for setting inhomogeneous Neumann conditions at
  // the wall for the mean-pressure Poisson equation and for setting
  // wall-boundary conditions for particles
  //
  // arrays that store the global-to-local mapping of column indices of
  // block-matrix that stores nonzero matrix elements of wall-boundary points of
  // matrix R
  if ( !(wlc = (int*)calloc(npoin,sizeof(int))) )
    ERR("Can't allocate memory!");
  // determine dimensions of block-matrix:
  // column-dimension (nwnzc) = total number of wall points
  // row-dimension (nwnzr) = total number of nonrecurring indices of points
  // surrounding wall points
  for (nwnzr=nwnzc=i=0; i<nbpoin; i++)  // loop over all boundary points
     if (bptags[i*3+1] == 4) {          // wall
       nwnzc++;                         // found a wall point
       // add up number of points surrounding wall point i
       nwnzr += bpsup2[i+1]-bpsup2[i];
     }
  // allocate arrays to store row/column indices of nonzero matrix elements of
  // wall-points of matrix R
  if ( !(wnzr = (int*)calloc(nwnzr,sizeof(int))) )
    ERR("Can't allocate memory!");
  // populate row-, and column-index arrays of nonzeros and generate
  // global-to-local arrays (after the following loop nwnzr will be smaller,
  // since it will only account for the nonrecurring indices)
  for (nwnzr=nwnzc=i=0; i<nbpoin; i++)  // loop over all boundary points
     if (bptags[i*3+1] == 4) {          // wall
       wlc[bpg[i]] = nwnzc++;           // store wall-local column-index
       // loop over points surrounding wall point i
       for (j=bpsup2[i]+1; j<=bpsup2[i+1]; j++) {
          // check if it has been encountered before
          for (haveit=k=0; k<nwnzr; k++)
            if (wnzr[k] == RDOF*bpsup1[j])
              haveit = 1;
          if (!haveit) {
            // store row index if it has not been encountered
            wnzr[nwnzr] = RDOF*bpsup1[j];
            nwnzr++;// incrase number of surrounding points around wall point i
          }
       }
     }

  ////////////////
  // count up number of wall-elements
  for (nwe=i=0; i<nbpoin; i++)  // loop over all boundary elements
    if (betags[i*2+0] == 4)     // wall
       nwe++;

  // allocate arrays to store wall-normals in wall-points and in wall-elements
  // and to store d<p>/dn in wall-elements in the previous timestep and to
  // store (d<Ui>/dt ni) in wall-elements and to store the indices of domain
  // elements at the wall and to store the node-index of the non-wall (outlier)
  // node of wall domain-elements and to store length of wall-elements and to
  // store edge-vector rBA of wall-elements and to store vector rA of
  // wall-elements and to store global index of node A of wall-elements and to
  // store closest wall-element index of domain-elements
  if ( !(wpnr = (double*)calloc(nwnzc*NDIMN,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(wenr = (double*)calloc(nwe*NDIMN,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(odpn = (double*)calloc(nwe,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(oudt = (double*)calloc(nwe,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(we = (int*)calloc(nwe,sizeof(int))) )
    ERR("Can't allocate memory!");
  if ( !(weo = (int*)calloc(nwe,sizeof(int))) )
    ERR("Can't allocate memory!");
  if ( !(wel = (double*)calloc(nwe,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(wrBA = (double*)calloc(nwe*NDIMN,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(wrA = (double*)calloc(nwe*NDIMN,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(wA = (int*)calloc(nwe,sizeof(int))) )
    ERR("Can't allocate memory!");
  if ( !(wec = (int*)calloc(nelem,sizeof(int))) )
    ERR("Can't allocate memory!");

  // initialize array wec with -1, meaning the element is not close to any wall
  for (i=0; i<nelem; i++) wec[i] = -1;

  // calculate wall-normals in wall-points and in wall-elements,
  // find domain-element indices corresponding to wall-elements,
  // find the node-index of the non-wall node of wall domain-elements,
  // calculate the length of wall-elements,
  // store edge-vector rBA of wall-elements,
  // store vector rA of wall-elements,
  // store index of node A of wall-elements,
  // fill up array that stores the closest wall-element of domain elements
  // (sharing two points with wall)
  for (nwe=i=0; i<nbpoin; i++)  // loop over all boundary elements
    if (betags[i*2+0] == 4) {   // wall
       // get node info and store index of node A
       A = wA[nwe] = binpoel[i*NBNODE+0];
       B = binpoel[i*NBNODE+1];

       odpn[nwe] = 0.0; // initialize old d<p>/dn as zero

       // store vector rA of wall-element
       wrA[nwe*NDIMN+0] = coord[A*NDIMN+0];
       wrA[nwe*NDIMN+1] = coord[A*NDIMN+1];
       // store edge-vector rBA of wall-element
       wrBA[nwe*NDIMN+0] = coord[B*NDIMN+0]-coord[A*NDIMN+0];
       wrBA[nwe*NDIMN+1] = coord[B*NDIMN+1]-coord[A*NDIMN+1];
       // calculate length of wall-element
       wel[nwe] = sqrt(wrBA[nwe*NDIMN+0]*wrBA[nwe*NDIMN+0] +
                       wrBA[nwe*NDIMN+1]*wrBA[nwe*NDIMN+1]);

       // calculate normal of boundary element
       nr[0] = coord[B*NDIMN+1] - coord[A*NDIMN+1];
       nr[1] = coord[A*NDIMN+0] - coord[B*NDIMN+0];
       dtmp = sqrt(nr[0]*nr[0] + nr[1]*nr[1]);
       nr[0] /= dtmp;
       nr[1] /= dtmp;

       wenr[nwe*NDIMN+0] = nr[0];       // store wall normal for wall-element
       wenr[nwe*NDIMN+1] = nr[1];
       // store wall normal for both points
       wpnr[wlc[A]*NDIMN+0] = wpnr[wlc[B]*NDIMN+0] = nr[0];
       wpnr[wlc[A]*NDIMN+1] = wpnr[wlc[B]*NDIMN+1] = nr[1];

       // search for domain-elements sharing two points with boundary element i
       for (j=0; j<nelem; j++) {
          jN = j*NNODE;

          // test if both A and B are points of element j 
          if ( ((A==inpoel[jN+0]) || (A==inpoel[jN+1]) || (A==inpoel[jN+2])) &&
               ((B==inpoel[jN+0]) || (B==inpoel[jN+1]) || (B==inpoel[jN+2])) ) {
            // connect domain-element to wall-element
            we[nwe] = j;
            // store wall-local index of domain-element (domain-element j is the
            // closest one to wall-element nwe)
            wec[j] = nwe;
            //printf("%d, %d\n",j,wec[j]);

            // find which node is not at the wall
            if (((A==inpoel[jN+0]) && (B==inpoel[jN+1])) ||
                 ((B==inpoel[jN+0]) && (A==inpoel[jN+1])))
              weo[nwe] = 2;
            if (((A==inpoel[jN+0]) && (B==inpoel[jN+2])) ||
                 ((B==inpoel[jN+0]) && (A==inpoel[jN+2])))
              weo[nwe] = 1;
            if (((A==inpoel[jN+1]) && (B==inpoel[jN+2])) ||
                 ((B==inpoel[jN+1]) && (A==inpoel[jN+2])))
              weo[nwe] = 0;

            j = nelem;  // get out
          }
       }

       nwe++;   // increase number of wall elements
     }
}

#ifdef MICROMIXING
static int findel( double xs, double ys )
// -----------------------------------------------------------------------------
// Routine: findel - returns index of 2d element of location (xs,ys)
//                   returns -1 if location not found in any element
//                   (this is only used for finding source location, etc. (for
//                   stuff that is run only once),
//                   a faster and highly optimized search is used for particle
//                   tracking)
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int e, A, B, C;
  double NA, NB, NC;
  double x[NNODE], y[NNODE];

  for (e=0; e<nelem; e++) {     // loop over all elements
     A = inpoel[e*NNODE+0];     // get node information of element e
     B = inpoel[e*NNODE+1];
     C = inpoel[e*NNODE+2];
     x[0] = coord[A*NDIMN+0];  y[0] = coord[A*NDIMN+1];
     x[1] = coord[B*NDIMN+0];  y[1] = coord[B*NDIMN+1];
     x[2] = coord[C*NDIMN+0];  y[2] = coord[C*NDIMN+1];

     // evaluate shapefunctions at particle location (Cramer's rule)
     NA = (xs*(y[1]-y[2]) + x[1]*(y[2]-ys) + x[2]*(ys-y[1]))/dete[e];
     NB = (x[0]*(ys-y[2]) + xs*(y[2]-y[0]) + x[2]*(y[0]-ys))/dete[e];
     NC = (x[0]*(y[1]-ys) + x[1]*(ys-y[0]) + xs*(y[0]-y[1]))/dete[e];

     if ((fmin(NA,1-NA)>0) && (fmin(NB,1-NB)>0) && (fmin(NC,1-NC)>0))
       return(e);
  }

  return(-1);
}
#endif

static void create_time_averaged_quantities( void )
// -----------------------------------------------------------------------------
// Routine: create_time_averaged_quantities - Allocate time-averaged quantities
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  // arrays for storing time-averaged quantities over the whole domain
  if ( !(tu = (double*)calloc(npoin*NDIMN,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tu2 = (double*)calloc(npoin*U2DOF,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tu3 = (double*)calloc(npoin*NDIMN,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tu4 = (double*)calloc(npoin*NDIMN,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tf = (double*)calloc(npoin,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tpr = (double*)calloc(npoin,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tdu = (double*)calloc(npoin*4,sizeof(double))) )
    ERR("Can't allocate memory!");
  #ifdef MICROMIXING
  if ( !(tc = (double*)calloc(npoin,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tc2 = (double*)calloc(npoin,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tc3 = (double*)calloc(npoin,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tc4 = (double*)calloc(npoin,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(ttm = (double*)calloc(npoin,sizeof(double))) )
    ERR("Can't allocate memory!");
  if ( !(tuc = (double*)calloc(npoin*3,sizeof(double))) )
    ERR("Can't allocate memory!");
  #endif
}

#ifdef MICROMIXING
static void allocate_space_for_pdfs( void )
// -----------------------------------------------------------------------------
// Routine: allocate_space_for_pdfs - Allocate space for pdfs at given locations
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int p, j;

  // count up number of locations to save scalar pdf and conditional statistics
  // at npl <= MAXPL, since the domain size could be smaller then the values in
  // PL, ie. not all locations specified in PL might be found on the grid
  for (npl=p=0; p<MAXPL; p++)
    if (findel(PL[p*2+0],PL[p*2+1]) != -1)
      npl++;

  // array to store the index of locations found
  if ( !(sl = (int*)calloc(npl,sizeof(int))) )
    ERR("Can't allocate memory!");
  // array to store element indices of pdf locations
  if ( !(epdfloc = (int*)calloc(npl,sizeof(int))) )
    ERR("Can't allocate memory!");
  // array to store time-averaged pdfs of temperature at different locations
  if ( !(tpdf = (double*)calloc(npl*NBI,sizeof(double))) )
    ERR("Can't allocate memory!");

  // find and store 2d elements of locations where pdfs and conditional
  // statistics of scalar will be saved
  for (npl=p=0; p<MAXPL; p++)
     if ((j=findel(PL[p*2+0],PL[p*2+1])) != -1) {
        epdfloc[npl] = j;       // store 2d element index
        sl[npl] = p;            // store location index (in array PL)
        npl++;                  // increase number of locations to save pdf at
     }
     else
       printf("PDF location not found in any element: %g,%g\n",PL[p*2+0],
                                                               PL[p*2+1]);
}
#endif

static void allocate_space_for_profiles( void )
// -----------------------------------------------------------------------------
// Routine: allocate_space_for_profiles - Allocate space for points to store
//                                        downstream profiles
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i, e, eN, p, A, B, C, n;

  for (maxn=ndl=i=0; i<MAXDL; i++) {    // loop over all downstream locations
                      // where cross-stream profiles are supposed to be saved
    for (n=e=0; e<nelem; e++) {     // loop over all elements
       eN = e*NNODE;                // get node-information
       A = inpoel[eN+0];
       B = inpoel[eN+1];
       C = inpoel[eN+2];

       // order nodes according to non-decreasing x coordinate in index-order
       // A, C, B
       if ( coord[B*NDIMN+0] < coord[A*NDIMN+0] ) { p=A; A=B; B=p; }
       if ( coord[C*NDIMN+0] < coord[A*NDIMN+0] ) { p=A; A=C; C=p; }
       if ( coord[B*NDIMN+0] < coord[C*NDIMN+0] ) { p=B; B=C; C=p; }

       if ((DL[i] > coord[A*NDIMN+0]) && (DL[i] < coord[B*NDIMN+0])) {
         n++;   // if in bounding box, it will have at least one point

         // also if an edge is aligned with the upper wall, it will have an
         // additional point
         if (fabs(coord[A*NDIMN+1]-maxy)<EPS) {
           if ((fabs(coord[B*NDIMN+1]-maxy)<EPS) ||
               (fabs(coord[C*NDIMN+1]-maxy)<EPS))
             n++;
         }
         else
           if ((fabs(coord[B*NDIMN+1]-maxy)<EPS) &&
               (fabs(coord[C*NDIMN+1]-maxy)<EPS))
             n++;
       }
    }

    if (n) {    // if domain contains the desired downstream location
      ndl++;    // increase number of downstream locations to be saved
                // find the largest number of points to be saved in one
                // downstream location
      if (n > maxn)
        maxn = n;
    }
  }

  // allocate arrays to store the points of cross-stream profiles
  if ( !(uprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!");  // <U>
  if ( !(vprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!");  // <V>
  if ( !(u2prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // <uu>
  if ( !(v2prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // <vv>
  if ( !(u3prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // skewness u
  if ( !(v3prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // skewness v
  if ( !(u4prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // flatness u
  if ( !(v4prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // flatness v
  if ( !(uvprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // <uv>
  if ( !(tkeprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!");// tke
  if ( !(epsprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!");// eps
  if ( !(yprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!");  // y coordinate
  #ifdef MICROMIXING
  if ( !(cprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!");  // scalar mean
  if ( !(c2prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // scalar variance
  if ( !(c3prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // scalar skewness
  if ( !(c4prof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // scalar kurtosis
  if ( !(tmprof = (double*)calloc(ndl*maxn,sizeof(double))) )
    ERR("can't allocate memory!"); // micromixing timescale
  #endif
}

static void allocate_space_for_centerline_profiles( void )
// -----------------------------------------------------------------------------
// Routine: allocate_space_for_centerline_profiles - Allocates space for points
//          to store centerline streamwise profiles downstream profiles
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int e, eN, p, A, B, C, n;

  for (n=e=0; e<nelem; e++) {    // loop over all elements
     eN = e*NNODE;               // get node-information
     A = inpoel[eN+0];
     B = inpoel[eN+1];
     C = inpoel[eN+2];

     // order nodes according to non-decreasing y coordinate in index-order
     // A, C, B
     if ( coord[B*NDIMN+1] < coord[A*NDIMN+1] ) { p=A; A=B; B=p; }
     if ( coord[C*NDIMN+1] < coord[A*NDIMN+1] ) { p=A; A=C; C=p; }
     if ( coord[B*NDIMN+1] < coord[C*NDIMN+1] ) { p=B; B=C; C=p; }

     if ((0.0 > coord[A*NDIMN+1]) && (0.0 < coord[B*NDIMN+1])) {
       n++;     // if in bounding box, it will have at least one point

     // also if an edge is aligned with the outflow, it will have an additional
     // point
     if (fabs(coord[A*NDIMN+0]-maxx)<EPS) {
       if ((fabs(coord[B*NDIMN+0]-maxx)<EPS) ||
           (fabs(coord[C*NDIMN+0]-maxx)<EPS))
         n++;
     } else
       if ((fabs(coord[B*NDIMN+0]-maxx)<EPS) &&
           (fabs(coord[C*NDIMN+0]-maxx)<EPS))
         n++;
     }
  }

  // allocate arrays to store the points of centerline streamwise profiles
  if ( !(uprof_c = (double*)calloc(n,sizeof(double))) )
    ERR("can't allocate memory!");  // <U>
  if ( !(pprof_c = (double*)calloc(n,sizeof(double))) )
    ERR("can't allocate memory!");  // 2<P>
  if ( !(xprof_c = (double*)calloc(n,sizeof(double))) )
    ERR("can't allocate memory!");  // x coordinate
}

static void gen_wall( void )
// -----------------------------------------------------------------------------
// Routine: gen_wall - Allocates space for point indices and angles along the
//                     cylinder wall
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int p, pN, e, itmp;
  double r, dtmp;

  // determine size of arrays
  for (angle_size=p=0; p<npoin; p++ ) { // loop over all points
    pN = p*NDIMN;

    if (coord[pN+1] > -EPS) {           // consider upper half only
       // compute square of radius
       r = sqrt(coord[pN+0]*coord[pN+0] + coord[pN+1]*coord[pN+1]);

       if (r-D/2 < EPS)        // if point is alongthe cylinder wall
         angle_size++;         // increase number of wall points
    }
  }

  // allocate arrays to store wall point indices, their angles and profiles
  if ( !(inprof_w = (int*)calloc(angle_size,sizeof(int))) )
    ERR("can't allocate memory!");      // point index
  if ( !(angle_w = (double*)calloc(angle_size,sizeof(double))) )
    ERR("can't allocate memory!");      // angle

  // fill up arrays
  for (angle_size=p=0; p<npoin; p++) {  // loop over all points
    pN = p*NDIMN;

    if (coord[pN+1] > -EPS) {           // consider upper half only
      // compute square of radius
      r = sqrt(coord[pN+0]*coord[pN+0] + coord[pN+1]*coord[pN+1]);

      if (r-D/2 < EPS) {                // if point is along the cylinder wall
        inprof_w[angle_size] = p;       // store point index
        // store angle
        angle_w[angle_size] = atan(coord[pN+1]/coord[pN+0])*180.0/PI;
        if (coord[pN+0] > 0) angle_w[angle_size] = 180-angle_w[angle_size];
        if (coord[pN+0] < 0) angle_w[angle_size] = -angle_w[angle_size];

        angle_size++;                   // increase number of wall points
      }
    }
  }

  // (bubble-)sort point indices of increasing angle
  for ( p=1; p<angle_size; p++)
    for (e=0; e<(angle_size-p); e++) {
      if (angle_w[e] > angle_w[e+1]) {
        SWAP(angle_w[e], angle_w[e+1], dtmp);
        SWAP(inprof_w[e], inprof_w[e+1], itmp);
      }
    }
}

static void gen_az( void )
// -----------------------------------------------------------------------------
// Routine: gen_az - Allocate linked lists that store the element indices
//                   intersecting azimuthal lines
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int p, e, eN, n, nN, m, A, B, C, NA, NB, NC, z, zN, itmp;
  double dtmp;
  double g[NDIMN], x[NNODE], y[NNODE];

  // determine size of linked list
  for (n=z=0; z<MAXAZ; z++) {   // loop over all azimuthal lines
    zN = z*NDIMN;

    // compute and store direction vector based on azimuth angle
    dir_az[zN+0] = -sin(PI/2 - AZ[z]*PI/180.0);
    dir_az[zN+1] = cos(PI/2 - AZ[z]*PI/180.0);

    for (e=0; e<nelem; e++) {   // loop over all elements
      eN = e*NNODE;             // get element information
      A = inpoel[eN+0];
      B = inpoel[eN+1];
      C = inpoel[eN+2];
      NA = A*NDIMN;
      NB = B*NDIMN;
      NC = C*NDIMN;

      // only consider upper half
      if ((coord[NA+1] > 0) && (coord[NB+1] > 0) && (coord[NC+1] > 0))
        // if any of the element edges intersects the azimuthal line
        if ((intersect(coord+NA, coord+NB, dir_az+zN, g)) ||
            (intersect(coord+NB, coord+NC, dir_az+zN, g)) ||
            (intersect(coord+NC, coord+NA, dir_az+zN, g)))
          n++;  // increase number of elements intersected
    }
  }

  // allocate linked list to store element indices intersected by azimuthal
  // lines, radial coordinates, evaluated shapefunctions at intersections and
  // profiles
  if ( !(az1 = (int*)malloc(n*sizeof(int))) )
    ERR("can't allocate memory!");
  if ( !(az2 = (int*)malloc((MAXAZ+1)*sizeof(int))) )
    ERR("can't allocate memory!");
  if ( !(r_az = (double*)malloc(n*sizeof(double))) )
    ERR("can't allocate memory!");
  if ( !(N_az = (double*)malloc(n*NNODE*sizeof(double))) )
    ERR("can't allocate memory!");

  // fill up linked list with element indices that intersect the azimuthal lines
  for (az2[0]=-1,n=z=0; z<MAXAZ; z++) { // loop over all azimuthal lines
    zN = z*NDIMN;

    for (m=e=0; e<nelem; e++) {         // loop over all elements
      eN = e*NNODE;                     // get element information
      A = inpoel[eN+0];
      B = inpoel[eN+1];
      C = inpoel[eN+2];
      NA = A*NDIMN;
      NB = B*NDIMN;
      NC = C*NDIMN;

      // only consider upper half
      if ((coord[NA+1] > 0) && (coord[NB+1] > 0) && (coord[NC+1] > 0))
        // if any of the element edges intersects the azimuthal line
        if ((intersect(coord+NA, coord+NB, dir_az+zN, g)) ||
            (intersect(coord+NB, coord+NC, dir_az+zN, g)) ||
            (intersect(coord+NC, coord+NA, dir_az+zN, g))) {
          // store radial coordinate of intersection point
          r_az[n] = sqrt(g[0]*g[0] + g[1]*g[1]);
          // store element index
          az1[n] = e;

          // extract node coordinates
          x[0] = coord[NA+0];  y[0] = coord[NA+1];
          x[1] = coord[NB+0];  y[1] = coord[NB+1];
          x[2] = coord[NC+0];  y[2] = coord[NC+1];

          // evaluate and store shapefunctions at intersection (Cramer's rule)
          nN = n*NNODE;
          N_az[nN+0] = (g[0]*(y[1]-y[2]) +
                        x[1]*(y[2]-g[1]) +
                        x[2]*(g[1]-y[1]))/dete[e];
          N_az[nN+1] = (x[0]*(g[1]-y[2]) +
                        g[0]*(y[2]-y[0]) +
                        x[2]*(y[0]-g[1]))/dete[e];
          N_az[nN+2] = (x[0]*(y[1]-g[1]) +
                        x[1]*(g[1]-y[0]) +
                        g[0]*(y[0]-y[1]))/dete[e];

          n++;// increase total number of elements intersected (for all angles)
          m++;// increase number of elements intersected (for angle z)
        }
    }
    az2[z+1] = n-1;     // store next-azimuth index

    // (bubble-)sort profile-points to increasing order in radial coordinate
    // (for angle z)
    for (p=1; p<m; p++)
      for (e=0; e<(m-p); e++)
        if (r_az[n-m+e] > r_az[n-m+e+1]) {
          SWAP(r_az[n-m+e], r_az[n-m+e+1], dtmp); // swap radial coordinate
          SWAP(az1[n-m+e], az1[n-m+e+1], itmp);   // update element index
          // swap shapefunctions
          SWAP(N_az[(n-m+e)*NNODE+0], N_az[(n-m+e+1)*NNODE+0], dtmp);
          SWAP(N_az[(n-m+e)*NNODE+1], N_az[(n-m+e+1)*NNODE+1], dtmp);
          SWAP(N_az[(n-m+e)*NNODE+2], N_az[(n-m+e+1)*NNODE+2], dtmp);
        }
  }
}

static void out_wall_meshinfo( void )
// -----------------------------------------------------------------------------
// Routine: out_wall_meshinfo - Output mesh info about the wall resolution
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i, nwe;
  double minsize_az, maxsize_az, minsize_ra, maxsize_ra, size_ra;
  double wrCA[NDIMN];

  minsize_az = minsize_ra = 1.0e+5;
  maxsize_az = maxsize_ra = 0.0;
  for (nwe=i=0; i<nbpoin; i++)  // loop over all boundary elements
    if (betags[i*2+0] == 4) {   // wall
      // compute edge vector CA
      wrCA[0] = coord[inpoel[we[nwe]*NNODE+2]*NDIMN+0] -
                coord[inpoel[we[nwe]*NNODE+0]*NDIMN+0];
      wrCA[1] = coord[inpoel[we[nwe]*NNODE+2]*NDIMN+1] -
                coord[inpoel[we[nwe]*NNODE+0]*NDIMN+1];

      // compute radial size of wall element
      size_ra = fabs(wrBA[nwe*NDIMN+0]*wrCA[1] -
                     wrBA[nwe*NDIMN+1]*wrCA[0]) / wel[nwe];

      if ( wel[nwe] < minsize_az )      // find smallest azimuthal length
        minsize_az = wel[nwe];          // store smallest azimuthal length

      if ( size_ra < minsize_ra )       // find smallest radial length
        minsize_ra = size_ra;           // store smallest radial length

      if ( wel[nwe] > maxsize_az )      // find largeest azimuthal length
        maxsize_az = wel[nwe];          // store largest azimuthal length

      if ( size_ra > maxsize_ra )       // find largest radial length
        maxsize_ra = size_ra;           // store largest radial length

      nwe++;    // increase number of wall-elements considered
    }

  // output: number of wall elements,
  //         smallest length in azimuthal and radial direction
  //         largest length in azimuthal and radial direction
  printf( " * wall elements: nwe = %d\n", nwe );
  printf( "                  minsize_az = %g\tminsize_ra = %g\n",
          minsize_az, minsize_ra );
  printf( "                  maxsize_az = %g\tmaxsize_ra = %g\n",
          maxsize_az, maxsize_ra );
}

static void Initialize()
// -----------------------------------------------------------------------------
// Routine: Initialize - Initialize
// Author : J. Bakosi
// -----------------------------------------------------------------------------
// - finds out the number of threads we have,
// - prepares unstructured mesh,
// - finds the element index of the source location,
// - finds out the element index of each boundary-measurement locations where
//   concentrations are sampled,
// - calculates bounds of domain,
// - creates vectors and sparse matrices,
// - creates arrays to store time-averaged quantities,
// - creates auxiliary structures for efficient access of boundary data,
// - check whether proper restartfiles exist for a restart
// - prepares random number generators and tables,
// - initially generates particles and allocates necessary data structures,
// - computes initial Eulerian statistics
// -----------------------------------------------------------------------------
{
  int itmp, k, samenthreads;
  char filename[STRLEN];
  FILE *restartfile;

  // query number of threads available
  #ifdef _OPENMP
  nthreads = omp_get_max_threads();
  #else
  nthreads = 1;
  #endif

  // prepare unstructured grid
  prepmsh(&npoin, &nbpoin, &nelem, &coord, &bpg, &binpoel, &inpoel, &esup1,
          &esup2, &psup1, &psup2, &bpsup1, &bpsup2, &esupel1, &esupel2, &esuel,
          &bptags, &betags, &Ae, &dNx, &dNy, &dete, &sqrtAp, &minsqrtAp );

  printf(" * number of threads: %d\n",nthreads);
  fflush(stdout);

  // calculate bounds of domain
  calcbounds();

  // create vectors and matrices for storing physical variables
  createvecmat();

  // allocate and initialize space for time-averaged quantities
  create_time_averaged_quantities();

  // create auxiliary global dynamic arrays
  create_bc_structs();

  // output some more mesh info
  out_wall_meshinfo();

  // check if restarted
  if ((restartfile = fopen(RESTART_FILENAME,"r"))) {
    restarted = 1;
    fclose(restartfile);

    // check for the rest of restartfiles
    for (itmp=k=0; k<MAXNTHREADS; k++) {
      #ifndef WALLFUNCTIONS
      // construct filename for uniform stream used in tables
      sprintf( filename, "%s.u.%d", RESTART_FILENAME, k );
      // check if file exists
      if ( (restartfile = fopen(filename,"r")) ) {
        itmp++;
        fclose(restartfile);
      }
      #endif

      // construct filename for Gaussian stream used in tables
      sprintf(filename, "%s.g.%d", RESTART_FILENAME, k);
      // check if file exists
      if ((restartfile = fopen(filename,"r"))) {
        itmp++;
        fclose(restartfile);
      }

      #ifndef WALLFUNCTIONS
      // construct filename for stream used for a few numbers at a time
      sprintf(filename, "%s.f.%d", RESTART_FILENAME, k);
      // check if file exists
      if ((restartfile = fopen(filename,"r"))) {
        itmp++;
        fclose(restartfile);
      }
      #endif
    }

    #ifndef WALLFUNCTIONS
    if (itmp != 3*nthreads) {
    #else
    if (itmp != nthreads) {
    #endif
      samenthreads = 0;
      printf(" o warning: nthreads differs from number of restartfiles\n");
      fflush(stdout);
    } else samenthreads = 1;
  }
  else restarted = 0;

  // initialize random number generator streams for sampling a few at a time
  preprng_streams(nthreads,
                  &stream      // <- allocated
                  #ifndef WALLFUNCTIONS
                  , restarted, samenthreads
                  #endif
                 );

  #ifdef MICROMIXING
  // allocate space for scalar pdfs at selected locations
  allocate_space_for_pdfs();
  #endif

  // allocate memory for profiles at downstream locations
  allocate_space_for_profiles();

  // allocate memory for profiles at streamwise centerline
  allocate_space_for_centerline_profiles();

  // generate linked lists that store the element indices of each stripe around
  // downstream locations
  genes( ndl, nelem, miny, maxy, inpoel, coord,
         &es1, &es2 );  // <- allocated

  // generate array that stores the element indices of the stripe around the
  // streamwise centerline
  genec( nelem, maxx, inpoel, coord,
         &ec,           // <- allocated
         &ec_size );    // <- modified

  // allocate array to store and extract point indices along the cylinder wall
  gen_wall();

  // allocate linked lists that store the element indices intersecting azimuthal
  // lines
  gen_az();

  // generate particles
  initgenpar( nelem, nthreads, inpoel, esuel, esupel1, esupel2, minx, maxx, Ae,
              coord, stream, &npar, &onum, &elp, &npel, &npeldt, &parcoord,
              &parvel, npoin, tu, tdu, tpr, &it, &it0, &t, &t0, &psel1, &psel2,
              &parfreq, tu2, tu3, tu4, tf
              #ifdef MICROMIXING
              , &parc, &partm,
              tc, tc2, tc3, tc4, tuc
              #ifdef VCIEM
              , &cp
              #endif
              #endif
            );

  // initialize random number generator streams for sampling whole tables at a
  // time
  preprng_tables(NGRPP*npar, nthreads, restarted, samenthreads,
                 #ifndef WALLFUNCTIONS
                 NURPP*npar, &ru,
                 #endif
                 &rg );

  // extract (initial) statistics from particles (throw away diagnostics)
  stat(nelem, npoin, nbpoin, npar, nthreads, 1.0, npel, esup1, esup2, inpoel,
       elp, nwe, wenr, ue, u, du, ddu, parvel, dNx, dNy, &itmp, betags, binpoel,
       oudt, u2e, u3e, u4e, fe, u2, u3, u4, f, parfreq
       #ifndef WALLFUNCTIONS
       , bptags, bpg, we, weo
       #endif
       #ifdef MICROMIXING
       , ce, c2e, c3e, c4e, tme, uce, c, c2, c3, c4, tm, uc, parc, partm
         #ifdef VCIEM
         , cp, vcte, psel1, psel2
         #endif
       #endif
      );
}

static void Finalize( void )
// -----------------------------------------------------------------------------
// Routine: Finalize - Finalize
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  free( coord );
  free( binpoel );
  free( inpoel );
  free( bptags );
  free( betags );
  free( esup1 );
  free( esup2 );
  free( psup1 );
  free( psup2 );
  free( bpsup1 );
  free( bpsup2 );
  free( bpg );
  free( esupel1 );
  free( esupel2 );
  free( esuel );
  free( Ae );
  free( dNx );
  free( dNy );
  free( dete );
  free( sqrtAp );

  #ifndef WALLFUNCTIONS
  free( rho );
  free( rrhs );
  #endif
  free( ddu );
  free( u );
  free( f );
  free( pr );
  free( prhs );
  free( dpr );

  #ifdef MICROMIXING
  free( c );
  free( c2 );
  free( c3 );
  free( c4 );
  free( tm );
  free( uc );
  free( ce );
  free( c2e );
  free( c3e );
  free( c4e );
  free( tme );
  free( uce );
  #ifdef VCIEM
  free( vcte );
  #endif
  #endif

  free( du );
  free( ue );
  free( u2 );
  free( u3 );
  free( u4 );
  free( fe );
  free( u2e );
  free( u3e );
  free( u4e );
  #ifndef WALLFUNCTIONS
  destroy_sparsemat( &R );
  #endif
  destroy_sparsemat( &P );

  free( tu );
  free( tu2 );
  free( tu3 );
  free( tu4 );
  free( tf );
  free( tpr );
  free( tdu );
  #ifdef MICROMIXING
  free( tc );
  free( tc2 );
  free( tc3 );
  free( tc4 );
  free( ttm );
  free( tuc );
  #endif

  destroyrng_tables( nthreads,
                     #ifndef WALLFUNCTIONS
                     &ru,
                     #endif
                     &rg );
  destroyrng_streams( nthreads, &stream );

  free( ifl );
  free( ofl );
  free( ifpp );
  free( ofpp );
  free( ifnz );
  free( ofnz );
  free( wlc );
  #ifndef WALLFUNCTIONS
  free( wnz );
  #endif
  free( wnzr );
  free( wpnr );
  free( wenr );
  free( odpn );
  free( oudt );
  free( we );
  free( weo );
  free( wel );
  free( wrBA );
  free( wrA );
  free( wA );
  free( wec );

  free( parcoord );
  free( parvel );
  free( parfreq );
  #ifdef MICROMIXING
  free( parc );
  free( partm );
  #ifdef VCIEM
  free( cp );
  #endif
  #endif
  free( elp );
  free( npel );
  free( npeldt );
  free( psel1 );
  free( psel2 );
  destroy_tmp_par();

  #ifdef MICROMIXING
  free( sl );
  free( epdfloc );
  free( tpdf );
  free( cprof );
  free( c2prof );
  free( c3prof );
  free( c4prof );
  free( tmprof );
  #endif
  free( es1 );
  free( es2 );
  free( ec );
  free( uprof );
  free( vprof );
  free( u2prof );
  free( v2prof );
  free( u3prof );
  free( v3prof );
  free( u4prof );
  free( v4prof );
  free( uvprof );
  free( tkeprof );
  free( epsprof );
  free( yprof );

  free( xprof_c );
  free( uprof_c );
  free( pprof_c );

  free( inprof_w );
  free( angle_w );

  free( az1 );
  free( az2 );
  free( r_az );
  free( N_az );
}

int main( void )
// -----------------------------------------------------------------------------
// Routine: Main - Main
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  Initialize();
  Solve();
  Finalize();

  return(0);
}
