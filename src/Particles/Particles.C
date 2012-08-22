// -----------------------------------------------------------------------------
// \file    src/Particles/Particles.h
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Functions related to particles in 1D and 2D
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------
//
//  Functions related to Lagrangian particles in 1d and 2d.
//  
//  Together with its header-file this file contains:
//   *  assignment of initial particle properties,
//   *  initial generation of particles into 1d and 2d elements,
//   *  several particle search algorithms for 1d and 2d grids,
//   *  diagnostic functions,
//   *  routines for the (re-)generation of dynamically changing linked lists,
//   *  a function to generate an additional derived data structure (a linked
//      list) that store 2d elements surrounding nodes of a 2d element (this is
//      used in the efficient access of memory in the particle search algorithm),
//   *  functions to intially generate particles,
//   *  functions that implement the different restart capabilities
//      (saving/loading of restartfiles),
//   *  functions that generate additional data structures for the easy access
//      of stripe-related information (exploiting downstream homogeneity),
//   *  functions that compute statistics for velocity and scalar (including
//      velocity-conditioned mean for VCIEM using the projection method of Fox)
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <cmath>
#include "mkl.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Const.h"
#include "Macros.h"
#include "Particles.h"
#include "SparseMatrix.h"
#include "Matrix3.h"
#include "Random.h"
#include "RandomErrors.h"
#include "Overlap.h"
#include "Sort.h"

static int /*_parid, *_elp;*/ *_npels, *_npeli;
static double /**_parcoord, *_parvel, *_parfreq, *_parc, *_dist,*/ *_v1e;
#ifdef MICROMIXING
static double *_v2e;
#endif
static FILE *restartfile;

static void assign_initial_stats_from_scratch(int i, int myid, double *parvel,
                                              double *parfreq,
                                              VSLStreamStatePtr *stream
	                                      #ifdef MICROMIXING
					      , double *parc
					      #endif
					     )
// -----------------------------------------------------------------------------
// Routine: assign_initial_stats_from_scratch - Assign initial statistics from
//                                              scratch to particle i
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int i3 = i*3;
  double mean, var;
  double cov[9], x[3];

  // assign particle velocity with given mean and covariance matrix
  mean = 1.0;
  cov[0] = 0.01/3.0;
  cov[4] = 0.01/3.0;
  cov[8] = 0.01/3.0;
  cov[1] = cov[3] = 0.0;
  cov[2] = cov[5] = cov[6] = cov[7] = 0.0;
  m3cholesky( cov );
  // obtaing three standard Gaussian random numbers from table
  CheckVslError( vdRngGaussian(GAUSSIAN_METHOD, stream[myid], 3, x, 0.0, 1.0) );
  parvel[i3+0] = mean + cov[0]*x[0] + cov[1]*x[1] + cov[2]*x[2];
  parvel[i3+1] = cov[3]*x[0] + cov[4]*x[1] + cov[5]*x[2];
  parvel[i3+2] = cov[6]*x[0] + cov[7]*x[1] + cov[8]*x[2];

  // assign turbulent frequency: gamma with prescribed mean (1.0) and
  // variance (0.25)
  mean = 1.0;  var = 0.25;
  CheckVslError( vdRngGamma(GAMMA_METHOD, stream[myid], 1, parfreq+i,
                 mean*mean/var, 0.0, var/mean) );

  #ifdef MICROMIXING
  // assign particle scalar
  parc[i] = 0.0;
  #endif
}

static void genpar(int e, int i, int myid, int *inpoel, double *parcoord,
                   double *parvel, double *parfreq,
                   double *coord, double *Ae, VSLStreamStatePtr *stream
	           #ifdef MICROMIXING
		   , double *parc
                   #endif
		  )
// -----------------------------------------------------------------------------
// Routine: genpar - Generate a particle (whose index will be i) into element e
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int success=0, A, B, C;
  double pos[NDIMN], NA;

  // get node information of element e
  A = inpoel[e*NNODE+0];
  B = inpoel[e*NNODE+1];
  C = inpoel[e*NNODE+2];

  do {
    // get two uniformly distributed random numers between [0...1)
    CheckVslError(vdRngUniform(UNIFORM_METHOD, stream[myid], 2, pos, 0.0, 1.0));

     NA = 1.0-pos[0]-pos[1];	// evaluate local shapefunction of node A

     if ( fmin(NA,1-NA) > 0 ) {	// if (xp,yp) falls into element, keep it...
        success = 1;		// get out of while cycle

        // store global coordinates
	parcoord[i*2+0] = (1.0-pos[0]-pos[1])*coord[A*NDIMN+0]
			  + pos[0]*coord[B*NDIMN+0]
			  + pos[1]*coord[C*NDIMN+0];
	parcoord[i*2+1] = (1.0-pos[0]-pos[1])*coord[A*NDIMN+1]
			  + pos[0]*coord[B*NDIMN+1]
			  + pos[1]*coord[C*NDIMN+1];

        // assign initial statistics to particle
        assign_initial_stats_from_scratch(i, myid, parvel, parfreq, stream
	                                  #ifdef MICROMIXING
					  , parc
					  #endif
					 );
     }
  } while (!success);

  // test if particle i is indeed in element e
  // (this should not be necessary, however, for some weird reason, I
  // encountered situations when the generated particle positions would pass the
  // above loop (coming out with success=1) and the following test still would
  // fail, this is rare, eg. a single particle among 7M particles and only in
  // parallel and only with the VSL_BRNG_MCG59 generator, using VSL_BRNG_MCG31
  // has been okay so far (which is the default), but let's just keep this test
  // here, it's inexpensive and only done once during initialization anyway)
  double a1, a2;
  if (!PARINEL(i,e,a1,a2,coord,inpoel,parcoord,Ae)) printf("N");
}

static int parinelx(int p, int *e, int nelem, int *nexte, double minx,
                    double maxx, int *inpoel, int *esuel, int *esupel1,
                    int *esupel2, double *coord, double *parcoord, double *Ae,
		    int cry, int *bf )
// -----------------------------------------------------------------------------
// Routine: parinelx - particle search based on cross-products
// Author : J. Bakosi
// -----------------------------------------------------------------------------
// returns 1 if particle p is in element e, 0 if not
// increases bf if particle had to be brute-forced (only for diagnostics)
//
// if particle not found in element e, in nexte it returns the neighboring
// element of e in the direction of the particle (ie. the nexte element is
// the hint to check as next)
// -----------------------------------------------------------------------------
{
  int A, B, C, AN, BN, CN, eN, i;
  double NA, NC, xp, yp, a1, a2;
  double x[NNODE], y[NNODE];

  // get node-, and particle position-information
  eN = (*e)*NNODE;
  A = inpoel[eN+0];
  B = inpoel[eN+1];
  C = inpoel[eN+2];
  AN = A*NDIMN;
  BN = B*NDIMN;
  CN = C*NDIMN;
  x[0] = coord[AN+0];  y[0] = coord[AN+1];
  x[1] = coord[BN+0];  y[1] = coord[BN+1];
  x[2] = coord[CN+0];  y[2] = coord[CN+1];
  xp = parcoord[p*2+0];  yp = parcoord[p*2+1];
  
  // evaluate two cross products
  // that's enough to decide whether the particle is in the element or not
  // this is the same as evaluating the shapefunctions at the particle location
  NA = (x[2]-x[1])*(yp-y[1]) - (xp-x[1])*(y[2]-y[1]);
  NC = (xp-x[1])*(y[0]-y[1]) - (x[0]-x[1])*(yp-y[1]);

  if ( (NA>0) && (NC>0) && (0.5*(NA+NC)<Ae[*e]) )
    return(1);				// particle found
  else					// particle not found
  {
    double NB = Ae[*e]-0.5*(NA+NC);	// compute third shapefunction
    
    if ( NA < NB )			// find out which direction to go next
    {
      if ( NA < NC )
        *nexte = esuel[eN+0];
      else
        *nexte = esuel[eN+2];
    }
    else
    {
      if ( NB < NC )
        *nexte = esuel[eN+1];
      else
        *nexte = esuel[eN+2];
    }

    // if we are at the boundary, things are not so cool there,
    // so we make sure the particle is found, even with brute-force
    // if nothing else helps
    if ( *nexte == -1 ) {
      // look for particle in the elements surrounding the nodes of element e
      for ( i = esupel2[*e]+1; i <= esupel2[(*e)+1]; i++ )
        if (PARINEL(p,esupel1[i],a1,a2,coord,inpoel,parcoord,Ae)) {
          // now that we will get out as the particle found,
          // make sure e will contain the element of the particle
	  *e = esupel1[i];
	  return(1);
	}

      // if we got here, the particle still hasn't been found,
      // so fall back to brute-force
      for ( i = 0; i < nelem; i++ )
        if (PARINEL(p,i,a1,a2,coord,inpoel,parcoord,Ae)) {
          if (cry && (fabs(xp-minx-0.0001)>EPS) && (fabs(xp-maxx+0.0001)>EPS))
            (*bf)++;	// increase number of brute-forced particles
	  *e = i;	// now that we will get out as the particle found,
	  return(1);	// make sure e will contain the element of the particle
	}

      printf("%d,%d: %g,%g\n",*e,*nexte,xp,yp);
      ERR("Particle not found!");
    }

    return(0);
  }
}

int findpar(int p, int e, int myid, int nelem, double minx, double maxx,
            int *npeldt, int *elp, int *inpoel, int *esuel, int *esupel1,
            int *esupel2, double *coord, double *parcoord, double *Ae, int *bf )
// -----------------------------------------------------------------------------
// Routine: findpar - Find particle in 2d grid and update elp, npeldt
// Author : J. Bakosi
// -----------------------------------------------------------------------------
// this should only be called if particle p hasn't been found in its element
// p - particle number
// e - element where the particles has last been seen
// if particle found, updates elp, npeldt and returns 1, otherwise returns 0
// -----------------------------------------------------------------------------
{
  int nexte;

  npeldt[myid*nelem+e]--;	// register leaving particle from element e
  
  // look for particle until it's found
  while (!(parinelx(p, &e, nelem, &nexte, minx, maxx, inpoel, esuel, esupel1,
                    esupel2, coord, parcoord, Ae, 1, bf))) {
    e = nexte;                  // try next best element
  }

  elp[p] = e;                   // store new element index
  npeldt[myid*nelem+e]++;       // register arriving particle in element e

  return(1);
}

int minnpel( int nelem, int *npel )
// -----------------------------------------------------------------------------
// Routine: minnpel - Return the minimum number of particles / element
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int e, min;

  // find minimum number of particles/element
  min = npel[0];
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( e = 1; e < nelem; e++ )
    if ( npel[e] < min )
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      {
        if ( npel[e] < min )
          min = npel[e];
      }

  return( min );
}

/* not needed right now
static int maxnpel( int nelem, int *npel )
// -----------------------------------------------------------------------------
// Routine: minnpel - Return the maxmum number of particles / element
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int e, max;

  // find maximum number of particles/element
  max = npel[0];
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( e = 1; e < nelem; e++ )
    if ( npel[e] > max )
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      {
        if ( npel[e] > max )
          max = npel[e];
      }

  return( max );
}*/

void genpsel(int nelem, int npar, int nthreads, int *npeldt, int *npel,
             int *elp, int *psel1, int *psel2)
// -----------------------------------------------------------------------------
// Routine: genpsel - Generate psel1, psel2, static npel based on dynamic
//                    npeldt and elp,
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  int e, p, np;

  // generate psel2
  psel2[0] = -1;
  for (e=0; e<nelem; e++) {             // loop over all elements
     for (np=p=0; p < nthreads; p++)    // collect changes from all threads
       np += npeldt[p*nelem+e];

     psel2[e+1] = psel2[e] + np;	// store element index
  }

  // zero static npel
  ipzero(npel, nelem, nthreads);

  // generate psel1 and static npel
  for (p=0; p<npar; p++) {      // loop over all particles
    e = elp[p];                 // get current element number of particle p and
    npel[e]++;                  // increase number of particles/element e
    psel1[psel2[e]+npel[e]] = p;// store particle index in new location of psel1
  }
}

static void moveapar(int i, int j, int nelem, int npar, int nthreads,
                     int *npeldt, int *elp, int *npel, double *parcoord,
                     double *parvel, int *psel1, int *psel2, double *parfreq
                     #ifdef MICROMIXING
	             , double *parc
	             #endif
	            )
// -----------------------------------------------------------------------------
// Routine: moveapar - Moves a particle from element i to element j
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{   
  int pi, pj;

  // get the index of the last particle from element i
  pi = psel1[psel2[i+1]];
  // get a particle from element j (can be any, get the first one)
  pj = psel1[psel2[j]+1];

  // copy particle pj to pi
  parcoord[pi*2+0] = parcoord[pj*2+0];
  parcoord[pi*2+1] = parcoord[pj*2+1];
  parvel[pi*3+0] = parvel[pj*3+0];
  parvel[pi*3+1] = parvel[pj*3+1];
  parvel[pi*3+2] = parvel[pj*3+2];
  parfreq[pi] = parfreq[pj];
  #ifdef MICROMIXING
  parc[pi] = parc[pj];
  #endif

  npeldt[i]--;          // take out particle from element i
  npeldt[j]++;          // put particle into element j
  elp[pi] = j;          // store particle's new element number

  // update psel, npel
  genpsel( nelem, npar, nthreads, npeldt, npel, elp, psel1, psel2 );
}

void improvepar( int nelem, int npar, int nthreads, int *cnt, int *npeldt, int *elp, int *npel,
                 double *parcoord, double *parvel, int *psel1, int *psel2, double *parfreq
                 #ifdef MICROMIXING
	         , double *parc
	         #endif
	       )
//  
// improves spatial distribution of particles by taking out from the densest area
// and putting into the sparsest 
//
// returns (in cnt) the number of improvement cycles called
//
{
  int e, mine, minnp, maxe, maxnp;

  
  (*cnt) = 0;
  do
  {
    // find element containing the smallest and biggest number of particles
    minnp = maxnp = npel[0];
    mine = maxe = 0;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( e = 1; e < nelem; e++ )
    {
      if ( npel[e] < minnp )
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
          if ( npel[e] < minnp )
	  {
	    minnp = npel[e];		// store smallest number of particles/element and
	    mine = e;			// store its element index
	  }
	}

      if ( npel[e] > maxnp )
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
          if ( npel[e] > maxnp )
	  {
	    maxnp = npel[e];		// store largest number of particles/ element and
	    maxe = e;			// store its element index
	  }
	}
    }

    if ((minnp < MINNPEL) && (minnp!=maxnp))
    {
       // increase number of improvement cycles (only for diagnostics)
       (*cnt)++;
  
       // move a particle from element maxe to element mine
       moveapar( maxe, mine, nelem, npar, nthreads, npeldt, elp, npel,
                 parcoord, parvel, psel1, psel2, parfreq
                 #ifdef MICROMIXING
	         , parc
	         #endif
	       );
    }

  } while ((minnp < MINNPEL) && (minnp!=maxnp));
}








void improvepar2( int nelem, int npar, int nthreads, int *cnt, int *npeldt, int *elp, int *npel,
                  double *parcoord, double *parvel, int *psel1, int *psel2, double *parfreq
                  #ifdef MICROMIXING
	          , double *parc
	          #endif
	        )
//  
// improves spatial distribution of particles by taking out from the densest area
// and putting them into the sparsest 
//
// returns (in cnt) the number of particles moved
//
{
  int e, p, i, j, pi, pj, nce, c, myid;


  // make a copy of array npel and its indices
  // also add up number of critical elements (the ones that have less particles than MINNPEL)
  nce = 0;
  #ifdef _OPENMP
  #pragma omp parallel for private(e) reduction(+:nce)
  #endif
  for ( e = 0; e < nelem; e++ )
  {
   if ( npel[e] < MINNPEL ) nce++;	// add up number of critical elements
    _npels[e] = npel[e];		// copy array element
    _npeli[e] = e;			// store index
  }


  // sort array npels and drag along the indices
  quickSort_i1( _npels, _npeli, nelem );

  // now we have the elements with the least number of particles at the bottom of array _npels
  // and the elements with the most number of particles at the top


  // redistribute particles from the top to the bottom
  c = 0;
  #ifdef _OPENMP
  #pragma omp parallel private(myid,e,p,i,j,pi,pj) reduction(+:c)
  #endif
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( e = 0; e < nce; e++ )	// loop over critical elements from the bottom up until we reach MINNPEL
    {
      for ( p = 0; p < MINNPEL-_npels[e]; p++ )	// loop over the number of missing particles in each critical element
      {
        // get element index from the top where particle will be moved from
        i = _npeli[nelem-e-1];
        // get a particle index from the top (this one will be moved to the bottom)
        pi = psel1[psel2[i]+1+p];		// starting from the first one, get next

        // get element index in the bottom where particle will be moved to
        j = _npeli[e];
        // get a particle index in the bottom (whose properties will be used to initialize the newly arriwing particle)
        if (_npels[e] > 0)			// if there is at least one particle
          pj = psel1[psel2[j]+1+(p%_npels[e])];	// starting from the first get the next one, restart if exhausted
        else					// if there are no particles at all
          pj = pi;				// keep the properties of the newly arriving particle
      
        parcoord[pi*2+0] = parcoord[pj*2+0];	// copy particle properties
        parcoord[pi*2+1] = parcoord[pj*2+1];
        parvel[pi*3+0] = parvel[pj*3+0];
        parvel[pi*3+1] = parvel[pj*3+1];
        parvel[pi*3+2] = parvel[pj*3+2];
        parfreq[pi] = parfreq[pj];
        #ifdef MICROMIXING
        parc[pi] = parc[pj];
        #endif

        elp[pi] = j;				// store particle's new element number
      }

      npeldt[myid*nelem+i] -= p;		// take out p particles from top element
      npeldt[myid*nelem+j] += p;		// put p particles into bottom element
      c += p;					// increase total number of particles moved by p
    }
  }

  (*cnt) = c;					// return total number of particles moved

  // update psel, npel
  genpsel( nelem, npar, nthreads, npeldt, npel, elp, psel1, psel2 );
}









/*
static void reorderpar( int npar, int *elp, double *parcoord, double *parvel, double *parfreq,
                        double *parc )
//
// reorder particles to close spatial proximity
//
{
  int p, p2, p3;
  double dx, dy;
 

  printf( "reorderpar()...\n" );
  fflush(stdout);
  
  // make copies of particle properties
  #ifdef _OPENMP
  #pragma omp parallel for private(p,p2,p3)
  #endif
  for ( p = 1; p < npar; p++ )		// loop over all particles (the first one will stay)
  {
    p2 = p*2;
    p3 = p*3;
    _parid[p-1] = p;			// store particle index
    _elp[p] = elp[p];			// store particle element number
    _parcoord[p2+0] = parcoord[p2+0];	// store particle properties
    _parcoord[p2+1] = parcoord[p2+1];
    _parvel[p3+0] = parvel[p3+0];
    _parvel[p3+1] = parvel[p3+1];
    _parvel[p3+2] = parvel[p3+2];
    _parfreq[p] = parfreq[p];
    _parc[p] = parc[p];
  }

  
  // compute particle distances from first particle
  #ifdef _OPENMP
  #pragma omp parallel for private(p,p2,dx,dy)
  #endif
  for ( p = 1; p < npar; p++ )		// loop over all particles in the first bin
  {
    p2 = p*2;
    dx = parcoord[p2+0]-parcoord[0];
    dy = parcoord[p2+1]-parcoord[1];
    _dist[p-1] = dx*dx + dy*dy;
  }

  
  // sort distances of particles (and their indices)
  quickSort( _dist, _parid, npar-1 );

  
  // reorder particles
  #ifdef _OPENMP
  #pragma omp parallel for private(p,p2,p3)
  #endif
  for ( p = 1; p < npar; p++ )
  {
    p2 = p*2;
    p3 = p*3;
    parcoord[p2+0] = _parcoord[_parid[p-1]*2+0];
    parcoord[p2+1] = _parcoord[_parid[p-1]*2+1];
    parvel[p3+0] = _parvel[_parid[p-1]*3+0];
    parvel[p3+1] = _parvel[_parid[p-1]*3+1];
    parvel[p3+2] = _parvel[_parid[p-1]*3+2];
    parfreq[p] = _parfreq[_parid[p-1]];
    parc[p] = _parc[_parid[p-1]];
    elp[p] = _elp[_parid[p-1]];
  }
}
*/









static void initpar( int nelem, int *inpoel, double *Ae, VSLStreamStatePtr *stream,
                     double *coord, int *npar, int *onum, double **parcoord, double **parvel, double **parfreq
		     #ifdef MICROMIXING
		     , double **parc
		     #endif
		   )
//
// generates particles without restartfile
// allocates necessary data structures
//
{
  int e, p, myid;
 
  
  // not restarted from this simulation, we have no samples in time-averaging yet
  *onum = 0;

  // calculate total number of particles
  *npar = nelem*NPEL;

  printf(" * not using restartfile (NPEL = %d, npar = %d)\n", NPEL, *npar);

  // array for storing particle coordinates, velocities, turbulent frequencies and scalar...
  if ( !(*parcoord = (double*)malloc((*npar)*NDIMN*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(*parvel = (double*)malloc((*npar)*3*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(*parfreq = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  #ifdef MICROMIXING
  if ( !(*parc = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  #endif

  printf(" * generating particles...\n");
  fflush(stdout);


  // generate particles into each element
  #ifdef _OPENMP
  #pragma omp parallel private(myid,e,p)
  #endif
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for schedule(guided,100)
    #endif
    for ( e = 0; e < nelem; e++ )	// loop over all elements
      for ( p = 0; p < NPEL; p++ )	// generate NPEL particles into element e
        genpar( e, e*NPEL+p, myid, inpoel, *parcoord, *parvel, *parfreq, coord, Ae, stream 
	        #ifdef MICROMIXING
		, *parc
		#endif
	      );
  }
}







static void restartpar( int npoin, double *tu, double *tdu, double *tpr, double *tu2, double *tu3, double *tu4,
                        double *tf, double **parfreq,
                        #ifdef MICROMIXING
			double *tc, double *tc2, double *tc3, double *tc4, double *tuc,
			double **parc,
			#endif
			int *npar, int *onum, int *it, int *it0,
			double *t, double *t0, double **parcoord, double **parvel
		      )
//
// assign particle properties from restartfile
//
{
  char s[STRLEN];
  FILE *ftav;


  // read in npar, it0, t0
  fscanf( restartfile, "%d\t%d\t%lf\n", npar, it, t );

  // array for storing particle coordinates, velocities, turbulent frequencies and scalar...
  if ( !(*parcoord = (double*)malloc((*npar)*2*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(*parvel = (double*)malloc((*npar)*3*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(*parfreq = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  #ifdef MICROMIXING
  if ( !(*parc = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  #endif

  printf(" * full restart: %s (NPEL = n/a, npar = %d)\n",RESTART_FILENAME,*npar);
     
  *it0 = *it;				// save it0
  *t0 = *t;				// save t0
     
  if (*t0 > MAXTIME) printf(" * (t0=%g) > (MAXTIME=%g), WILL NOT DO ANY TIMESTEP!\n",*t,MAXTIME);
  if (*it0 > MAXTS) printf(" * (it0=%d) > (MAXTS=%d), WILL NOT DO ANY TIMESTEP!\n",*it0,MAXTS);
  fflush(stdout);
     
  fread( *parcoord, sizeof(double), (*npar)*2, restartfile );
  fread( *parvel, sizeof(double), (*npar)*3, restartfile );
  fread( *parfreq, sizeof(double), *npar, restartfile );
  #ifdef MICROMIXING
  fread( *parc, sizeof(double), *npar, restartfile );
  #endif

  // if they exist, read in time-averaged quantities as well
  // (if they don't, they'll start from zero)
  sprintf( s, "%s.tav", RESTART_FILENAME );
  if ( (ftav = fopen(s,"r")) )
  {
    fread( onum, sizeof(int), 1, ftav );
    printf(" * also restarting time-averaging after %d samples\n",*onum);
    fflush(stdout);
    fread( tu, sizeof(double), npoin*NDIMN, ftav );
    fread( tu2, sizeof(double), npoin*U2DOF, ftav );
    fread( tu3, sizeof(double), npoin*NDIMN, ftav );
    fread( tu4, sizeof(double), npoin*NDIMN, ftav );
    fread( tf, sizeof(double), npoin, ftav );
    fread( tpr, sizeof(double), npoin, ftav );
    fread( tdu, sizeof(double), npoin*4, ftav );
    #ifdef MICROMIXING
    fread( tc, sizeof(double), npoin, ftav );
    fread( tc2, sizeof(double), npoin, ftav );
    fread( tc3, sizeof(double), npoin, ftav );
    fread( tc4, sizeof(double), npoin, ftav );
    fread( tuc, sizeof(double), npoin*3, ftav );
    #endif
    fclose(ftav);
  }
  else	// not restarted from this simulation, we have no samples in time-averaging yet
    *onum = 0;

  fclose( restartfile );
}








void saverestartpar( int nthreads, int npoin, int npar, int onum, int it, double t, double *parcoord, double *parvel,
                     double *tu, double *tdu, double *parfreq, double *tu2, double *tu3, double *tu4, double *tf,
                     double *tpr
                     #ifndef WALLFUNCTIONS
                     , VSLStreamStatePtr *stream
                     #endif
		     #ifdef MICROMIXING
		     , double *parc, double *tc, double *tc2, double *tc3, double *tc4, double *tuc
		     #endif
		   )
//
// saves particle properties to restartfile
//
{
  char s[STRLEN];
 

  printf("saverestartpar()...\n");
  fflush(stdout);

  if (!(restartfile = fopen(RESTART_FILENAME,"w"))) ERR("cannot open RESTART_FILENAME");
  
  fprintf( restartfile, "%d\t%d\t%lf\n", npar, it, t );

  // output particle properties
  fwrite( parcoord, sizeof(double), npar*2, restartfile );
  fwrite( parvel, sizeof(double), npar*3, restartfile );
  fwrite( parfreq, sizeof(double), npar, restartfile );
  #ifdef MICROMIXING
  fwrite( parc, sizeof(double), npar, restartfile );
  #endif

  fclose( restartfile );


  // save time-averaged quantities as well (if time-averaging has been started)
  if (t > AVTIME)
  {
    sprintf( s, "%s.tav", RESTART_FILENAME );
    if (!(restartfile = fopen(s,"w"))) ERR("cannot open RESTART_FILENAME");
  
    fwrite( &onum, sizeof(int), 1, restartfile );
    fwrite( tu, sizeof(double), npoin*NDIMN, restartfile );
    fwrite( tu2, sizeof(double), npoin*U2DOF, restartfile );
    fwrite( tu3, sizeof(double), npoin*NDIMN, restartfile );
    fwrite( tu4, sizeof(double), npoin*NDIMN, restartfile );
    fwrite( tf, sizeof(double), npoin, restartfile );
    fwrite( tpr, sizeof(double), npoin, restartfile );
    fwrite( tdu, sizeof(double), npoin*4, restartfile );
    #ifdef MICROMIXING
    fwrite( tc, sizeof(double), npoin, restartfile );
    fwrite( tc2, sizeof(double), npoin, restartfile );
    fwrite( tc3, sizeof(double), npoin, restartfile );
    fwrite( tc4, sizeof(double), npoin, restartfile );
    fwrite( tuc, sizeof(double), npoin*3, restartfile );
    #endif

    fclose(restartfile);
  }


  // save random number streams
  saverng_streams( nthreads
                   #ifndef WALLFUNCTIONS
                   , stream
                   #endif
                 );
}






void genes( int ndl, int nelem, double miny, double maxy, int *inpoel, double *coord, int **es1, int **es2 )
//
// generates linked lists that store element indices of each stripe that surrounds a
// downstream location where profiles are saved
// allocates necessary data structures, these are freed in Finalize()
//
{
  int i, j, e, eN, A, B, C, NA, NB, NC;
  double p1[2], q1[2], r1[2], p2[2], q2[2], r2[2], pe[2], qe[2], re[2];


  // determine the size of es1
  for ( i=j=0; j < ndl; j++ )	// loop over all downstream locations
  {
     // calculate extent of stripes around a downstream location
     p1[0] = DL[j]-0.1;      p1[1] = miny;		// lower triangle vertices
     q1[0] = DL[j]+0.1;      q1[1] = miny;
     r1[0] = p1[0];          r1[1] = maxy;
     p2[0] = q1[0];          p2[1] = q1[1];		// upper triangle vertices
     q2[0] = q1[0];          q2[1] = maxy;
     r2[0] = r1[0];          r2[1] = r1[1];

     for ( e = 0; e < nelem; e++ )	// loop over all 2d elements
     {
	eN = e*NNODE;
        A = inpoel[eN+0];		// get node-information of element e
        B = inpoel[eN+1];
        C = inpoel[eN+2];
	NA = NDIMN*A;
	NB = NDIMN*B;
	NC = NDIMN*C;
	pe[0] = coord[NA+0];  pe[1] = coord[NA+1];
        qe[0] = coord[NB+0];  qe[1] = coord[NB+1];
        re[0] = coord[NC+0];  re[1] = coord[NC+1];
	
        if ( (tri_tri_overlap_test_2d(pe,qe,re,p1,q1,r1)) || (tri_tri_overlap_test_2d(pe,qe,re,p2,q2,r2)))
	  i++;	// increase element counter
     }
  }
  
  // linked list to store the element indices in each stripe
  if ( !(*es1 = (int*)malloc(i*sizeof(int))) ) ERR("Can't allocate memory!");
  if ( !(*es2 = (int*)malloc((ndl+1)*sizeof(int))) ) ERR("Can't allocate memory!");

  // fill linked lists es1, es2
  for ( (*es2)[0]=-1, i=j=0; j < ndl; j++ )	// loop over all stripes
  {
     // calculate extent of stripes around a downstream location
     p1[0] = DL[j]-0.1;      p1[1] = miny;	// lower triangle vertices
     q1[0] = DL[j]+0.1;      q1[1] = miny;
     r1[0] = p1[0];          r1[1] = maxy;
     p2[0] = q1[0];          p2[1] = q1[1];	// upper triangle vertices
     q2[0] = q1[0];          q2[1] = maxy;
     r2[0] = r1[0];          r2[1] = r1[1];

     for ( e = 0; e < nelem; e++ )	// loop over all 2d elements
     {
	eN = e*NNODE;
        A = inpoel[eN+0];		// get node-information of element e
        B = inpoel[eN+1];
        C = inpoel[eN+2];
	NA = NDIMN*A;
	NB = NDIMN*B;
	NC = NDIMN*C;
        pe[0] = coord[NA+0];  pe[1] = coord[NA+1];
        qe[0] = coord[NB+0];  qe[1] = coord[NB+1];
        re[0] = coord[NC+0];  re[1] = coord[NC+1];

	// if element e overlaps either lower or upper triangle of stripe j
        if ( (tri_tri_overlap_test_2d(pe,qe,re,p1,q1,r1)) || (tri_tri_overlap_test_2d(pe,qe,re,p2,q2,r2)))
	  (*es1)[i++] = e;		// store element
     }

     (*es2)[j+1] = i-1;			// store next-stripe index
  }
}









void genec( int nelem, double maxx, int *inpoel, double *coord, int **ec, int *size )
//
// generates array that stores element indices of the stripe that surrounds the streamwise centerline
// where profiles are saved
// allocates array, which is freed in Finalize()
//
{
  int e, eN, A, B, C, NA, NB, NC;
  double p1[2], q1[2], r1[2], p2[2], q2[2], r2[2], pe[2], qe[2], re[2];


  // determine the size of ec
  *size = 0;
  
  // calculate extent of stripes around the centerline
  p1[0] = maxx;      p1[1] = -0.1;	// lower triangle vertices
  q1[0] = maxx;      q1[1] = +0.1;
  r1[0] = 0.0;       r1[1] = -0.1;
  p2[0] = q1[0];     p2[1] = q1[1];	// upper triangle vertices
  q2[0] = r1[0];     q2[1] = q1[1];
  r2[0] = r1[0];     r2[1] = r1[1];

  for ( e = 0; e < nelem; e++ )		// loop over all 2d elements
  {
     eN = e*NNODE;
     A = inpoel[eN+0];			// get node-information of element e
     B = inpoel[eN+1];
     C = inpoel[eN+2];
     NA = NDIMN*A;
     NB = NDIMN*B;
     NC = NDIMN*C;
     pe[0] = coord[NA+0];  pe[1] = coord[NA+1];
     qe[0] = coord[NB+0];  qe[1] = coord[NB+1];
     re[0] = coord[NC+0];  re[1] = coord[NC+1];
	
     if ( (tri_tri_overlap_test_2d(pe,qe,re,p1,q1,r1)) || (tri_tri_overlap_test_2d(pe,qe,re,p2,q2,r2)))
       (*size)++;	// increase element counter
  }
  
  // array to store the element indices in the stripe
  if ( !(*ec = (int*)malloc((*size)*sizeof(int))) ) ERR("Can't allocate memory!");

  // calculate extent of stripes around the centerline
  p1[0] = maxx;      p1[1] = -0.1;	// lower triangle vertices
  q1[0] = maxx;      q1[1] = +0.1;
  r1[0] = 0.0;       r1[1] = -0.1;
  p2[0] = q1[0];     p2[1] = q1[1];	// upper triangle vertices
  q2[0] = r1[0];     q2[1] = q1[1];
  r2[0] = r1[0];     r2[1] = r1[1];

  *size = 0;

  for ( e = 0; e < nelem; e++ )		// loop over all 2d elements
  {
    eN = e*NNODE;
    A = inpoel[eN+0];			// get node-information of element e
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = NDIMN*A;
    NB = NDIMN*B;
    NC = NDIMN*C;
    pe[0] = coord[NA+0];  pe[1] = coord[NA+1];
    qe[0] = coord[NB+0];  qe[1] = coord[NB+1];
    re[0] = coord[NC+0];  re[1] = coord[NC+1];

    // if element e overlaps either lower or upper triangle of the stripe
    if ( (tri_tri_overlap_test_2d(pe,qe,re,p1,q1,r1)) || (tri_tri_overlap_test_2d(pe,qe,re,p2,q2,r2)))
      (*ec)[(*size)++] = e;		// store element index
  }
}










void initgenpar( int nelem, int nthreads, int *inpoel, int *esuel, int *esupel1, int *esupel2,
                 double minx, double maxx, double *Ae, double *coord, VSLStreamStatePtr *stream,
		 int *npar, int *onum, int **elp, int **npel, int **npeldt, double **parcoord,
		 double **parvel, int npoin, double *tu, double *tdu, double *tpr, int *it,
                 int *it0, double *t, double *t0, int **psel1, int **psel2, double **parfreq,
                 double *tu2, double *tu3, double *tu4, double *tf
		 #ifdef MICROMIXING
		 , double **parc, double **partm,
		 double *tc, double *tc2, double *tc3, double *tc4, double *tuc
		 #ifdef VCIEM
		 , int **cp
		 #endif
		 #endif
               )
//
// initially generates particles,
// assigns initial conditions,
// allocates particle-data structures
//
{
  int e, p, myid, nexte, bf;
 
  if ( (restartfile = fopen(RESTART_FILENAME,"r")) )	// if restartfile exists
    // assign particle properties from restartfile (full restart)
    restartpar( npoin, tu, tdu, tpr, tu2, tu3, tu4, tf,
                parfreq,				// <- allocated
                #ifdef MICROMIXING
		tc, tc2, tc3, tc4, tuc,
		parc,					// <- allocated
		#endif
                npar, onum, it, it0, t, t0,		// <- modified
		parcoord, parvel			// <- allocated
	      );
  else							// if restartfile doesn't exist
    // assign given initial properties to particles (dry start)
    initpar( nelem, inpoel, Ae, stream, coord,
             npar, onum,				// <- modified
	     parcoord, parvel, parfreq			// <- allocated
             #ifdef MICROMIXING
	     , parc					// <- allocated
	     #endif
	   );
  
  #ifdef MICROMIXING
    #ifdef VCIEM
    // array to store particle conditioning pointer to bins of velocity space
    if ( !(*cp = (int*)malloc((*npar)*sizeof(int))) ) ERR("Can't allocate memory!");
    #endif
  
  // array to store micromixing timescale for each particle
  if ( !(*partm = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  #endif
  
  // array for storing element-number of particles
  if ( !(*elp = (int*)malloc((*npar)*sizeof(int))) ) ERR("Can't allocate memory!");

  // arrays for storing number of particles in elements
  // npel: static array, that does not change during one timestep
  // npeldt: dynamic array, that reflects the changes of particle movements instantaneously
  if ( !(*npel = (int*)malloc(nelem*sizeof(int))) ) ERR("Can't allocate memory!");
  if ( !(*npeldt = (int*)malloc(nelem*nthreads*sizeof(int))) ) ERR("Can't allocate memory!");

  
  printf(" * initial search of particles...\n");
  fflush(stdout);
  
  // initial search of particles
  #ifdef _OPENMP
  #pragma omp parallel private(myid,e,nexte,p)
  #endif
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif

    // zero element counters
    memset( *npeldt + myid*nelem, 0, nelem*sizeof(int) );

    e = 0;					// start search from element 0
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( p = 0; p < *npar; p++ )		// loop over all particles
      if ( parinelx(p,&e,nelem,&nexte,minx,maxx,inpoel,esuel,esupel1,esupel2,coord,*parcoord,Ae,0,&bf) )
      {						// if particle p is in element e
         (*elp)[p] = e;				// store element index
         (*npeldt)[myid*nelem+e]++;		// increase number of particles/element e
						// (static npel will be generated from this in genpsel)
      }
      else { e = nexte; p--; }			// if particle not found, retry in next best element
  }

  // linked lists for storing the index of particles in each element
  if ( !(*psel1 = (int*)malloc((*npar)*sizeof(int))) ) ERR("Can't allocate memory!");
  if ( !(*psel2 = (int*)malloc((nelem+1)*sizeof(int))) ) ERR("Can't allocate memory!");
  
  // element-based temporary arrays for computation of statistics
  if ( !(_v1e = (double*)malloc(nthreads*10*nelem*sizeof(double))) ) ERR("can't allocate memory!");
  #ifdef MICROMIXING
  if ( !(_v2e = (double*)malloc(nthreads*6*nelem*sizeof(double))) ) ERR("can't allocate memory!");
  #endif

  // temporary arrays for particle redistribution
  if ( !(_npels = (int*)malloc(nelem*sizeof(int))) ) ERR("can't allocate memory!");
  if ( !(_npeli = (int*)malloc(nelem*sizeof(int))) ) ERR("can't allocate memory!");
  
  /*// temporary arrays for reordering particles
  if ( !(_parcoord = (double*)malloc((*npar)*NDIMN*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(_parvel = (double*)malloc((*npar)*3*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(_parfreq = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(_parc = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(_dist = (double*)malloc((*npar)*sizeof(double))) ) ERR("Can't allocate memory!");
  if ( !(_elp = (int*)malloc((*npar)*sizeof(int))) ) ERR("Can't allocate memory!");
  if ( !(_parid = (int*)malloc((*npar)*sizeof(int))) ) ERR("Can't allocate memory!");
  
  // reorder particles for better locality
  reorderpar( *npar, *elp, *parcoord, *parvel, *parfreq, *parc );*/

  // generate (initial) psel and npel
  genpsel( nelem, *npar, nthreads, *npeldt, *npel, *elp, *psel1, *psel2 );
}








#ifdef MICROMIXING

#ifdef VCIEM

#ifdef PROJECTION

static void compute_conditioned_mean_proj( int e, int *npel, int *cp, double *u2e, double *uce,
                                           double *parvel, double *parc, double *vcte,
                                           int *psel1, int *psel2 )
//
// computes velocity conditioned mean of scalar in element e
// using Fox's projection method
//
{
  int i, j, n, eU, e3;
  int np[CNBI], paridp[npel[e]];
  double det;
  double r[U2DOF], proj[3], parvelp[npel[e]];


  eU = e*U2DOF;
  e3 = e*3;

  // compute velocity-space pojection vector for 2d element e
  // rho_ij is the normalized Reynolds stress tensor
  det = (1.0 - u2e[eU+5]*u2e[eU+5])				// det(rho_ij)
        + u2e[eU+3]*(u2e[eU+4]*u2e[eU+5] - u2e[eU+3])
	+ u2e[eU+4]*(u2e[eU+3]*u2e[eU+5] - u2e[eU+4]);

  r[0] = (1.0 - u2e[eU+5]*u2e[eU+5])/det;			// inv(rho_ij)
  r[1] = (1.0 - u2e[eU+4]*u2e[eU+4])/det;
  r[2] = (1.0 - u2e[eU+3]*u2e[eU+3])/det;
  r[3] = (u2e[eU+4]*u2e[eU+5] - u2e[eU+3])/det;
  r[4] = (u2e[eU+3]*u2e[eU+5] - u2e[eU+4])/det;
  r[5] = (u2e[eU+4]*u2e[eU+3] - u2e[eU+5])/det;

  proj[0] = r[0]*uce[e3+0] + r[3]*uce[e3+1] + r[4]*uce[e3+2];	// projection vector
  proj[1] = r[3]*uce[e3+0] + r[1]*uce[e3+1] + r[5]*uce[e3+2];
  proj[2] = r[4]*uce[e3+0] + r[5]*uce[e3+1] + r[2]*uce[e3+2];

  // project velocity-space
  for ( n=0, j=psel2[e]+1; j <= psel2[e+1]; j++, n++ )
  {
    i = psel1[j]*3;
    paridp[n] = psel1[j];					// store index of particle
    parvelp[n] = proj[0]*parvel[i+0] + proj[1]*parvel[i+1] + proj[2]*parvel[i+2];
  }

  quickSort_d1( parvelp, paridp, npel[e] );		// sort projected velocities (and their indices)

  // compute (projected-)velocity-conditioned mean of scalar
  for ( i = 0; i < CNBI; i++ ) { np[i]=0; vcte[e*CNBI+i]=0.0; }
  for ( n=0, j=psel2[e]+1; j <= psel2[e+1]; j++, n++ )	// loop over particles in element e
  {
     i = CNBI*n/npel[e];				// compute bin index
     np[i]++;						// incrase number of particles in bin i
     vcte[e*CNBI+i] += parc[paridp[n]];			// add scalar to bin i
     cp[paridp[n]] = i;					// store conditioning pointer for particle
  }
  
  // finish calculating conditional mean in all CNBI bins, cry out if an empty bin is found (if this happens,
  // the number of particles in element e is too small compared to the number of bins we are trying to use,
  // solution: decrease CNBI and/or increase NPEL)
  for ( j=i=0; i < CNBI; i++ )
    if (np[i]!=0) vcte[e*CNBI+i] /= np[i]; else printf("*");
}

#else // PROJECTION

static int compute_conditioned_mean( int e, int *npel, int *cp, double *parvel, double *parc, double *vcte,
                                     int *psel1, int *psel2 )
//
// computes velocity conditioned mean of scalar in element e
// WITHOUT using Fox's projection method
//
// returns the number of conditioning bins actually used (cnbi)
//
{
  int i, j, n, cnbi;
  int paridc[npel[e]];
  double parvelc[3][npel[e]], np[CNBI];


  // extract particle velocities and indices into separate arrays (to make it easier to work with them)
  for ( n=0, j=psel2[e]+1; j <= psel2[e+1]; j++, n++ )
  {
    i = psel1[j]*3;
    parvelc[0][n] = parvel[i+0];			// copy out particle velocity
    parvelc[1][n] = parvel[i+1];
    parvelc[2][n] = parvel[i+2];
    paridc[n] = psel1[j];				// store index of particle for three dimensions
  }

  if ( npel[e] > BINd*MINNPBI )// see if we have enough particles to do the division in the 1st dimension
  {
    cnbi = BINd;
    // sort all particle indices in element e according to their x velocity
    quickSort_d3( // sort according to this array
                  parvelc[0],
		  // drag these arrays along when sorting (index, y, z coordinates)
		  paridc, parvelc[1], parvelc[2],
		  // the number of array elements to sort (size of array)
		  npel[e] );

    if ( npel[e] > BINd*BINd*MINNPBI )// see if we have enough particles to do the division in the 2nd dimension
    {
      cnbi = BINd*BINd;
      // sort particle indices in all BINd x-parts according to their y velocity
      for ( i = 0; i < BINd; i++ )	// loop over all x-parts
        quickSort_d3( // sort according to this sub-array
                      parvelc[1]+i*npel[e]/BINd,
		      // drag these sub-arrays along when sorting (index, x, z coordinates)
                      paridc+i*npel[e]/BINd, parvelc[0]+i*npel[e]/BINd, parvelc[2]+i*npel[e]/BINd,
		      // the number of array-elements to sort (size of sub-array)
		      npel[e]/BINd );

      if ( npel[e] > BINd*BINd*BINd*MINNPBI )// see if we have enough particles to do the division in the 3rd dimension
      {
        cnbi = BINd*BINd*BINd;
        // sort particle indices in all BINd x BINd y-parts according to their z velocity
        for ( i = 0; i < BINd*BINd; i++ )	// loop over all y-parts in each x-part
          quickSort_d3( // sort according to this sub-array
                        parvelc[2]+i*npel[e]/(BINd*BINd),
	                // drag these sub-arrays along when sorting (index, x, y coordinates)
                        paridc+i*npel[e]/(BINd*BINd), parvelc[0]+i*npel[e]/(BINd*BINd), parvelc[1]+i*npel[e]/(BINd*BINd),
	                // the number of array-elements to sort (size of sub-array)
		        npel[e]/(BINd*BINd) );
      }
    }
  } else cnbi = 1;	// couldn't even do the division in the 1st dimension -> IEM

  // at this point in cnbi we have the number of conditioning bins that is safe to use
  // if we want at least MINNPBI particles per bin (cnbi <= CNBI always)

  // compute velocity-conditioned mean of scalar
  for ( i = 0; i < CNBI; i++ ) { np[i] = vcte[e*CNBI+i] = 0.0; }
  for ( n=0, j=psel2[e]+1; j <= psel2[e+1]; j++, n++ )	// loop over particles in element e
  {
     i = cnbi*n/npel[e];				// compute bin index
     np[i] += 1.0;					// incrase number of particles in bin i
     vcte[e*CNBI+i] += parc[paridc[n]];			// add scalar to bin i
     cp[paridc[n]] = i;					// store conditioning pointer for particle
  }
  
  // finish calculating conditional mean in all cnbi bins, cry out if an empty bin is found (if this happens,
  // the number of particles in element e is too small compared to the number of bins we are trying to use,
  // solution: decrease BINd (or MINNPBI) and/or increase NPEL)
  for ( i = 0; i < cnbi; i++ )
    if (np[i] > EPS) vcte[e*CNBI+i] /= np[i]; else printf("*");

  // return the number of conditioning bins actually used
  return(cnbi);
}

#endif // PROJECTION

#endif // VCIEM

#endif // MICROMIXING














void destroy_tmp_par( void )
//
// destroys temporary arrays for computation of statistics
{
  // element-based temporary arrays
  free( _v1e );
  #ifdef MICROMIXING
  free( _v2e );
  #endif

  // temporary array for particle redistribution
  free( _npels );
  free( _npeli );

  /*// temporary arrays for reordering particles
  free( _parcoord );
  free( _parvel );
  free( _parfreq );
  free( _parc );
  free( _dist );
  free( _elp );
  free( _parid );*/
}







int stat( int nelem, int npoin, int nbpoin, int npar, int nthreads, double dt,
          int *npel, int *esup1, int *esup2, int *inpoel, int *elp,
	  int nwe, double *wenr, double *ue, double *u,
	  double *du, double *ddu, double *parvel,
	  double *dNx, double *dNy, int *nc_, int *betags, int *binpoel, double *oudt,
          double *u2e, double *u3e, double *u4e, double *fe, double *u2, double *u3,
          double *u4, double *f, double *parfreq
	  #ifndef WALLFUNCTIONS
          , int *bptags, int *bpg, int *we, int *weo
	  #endif
          #ifdef MICROMIXING
	  , double *ce, double *c2e, double *c3e, double *c4e, double *tme, double *uce,
	  double *c, double *c2, double *c3, double *c4, double *tm, double *uc, double *parc,
	  double *partm
	    #ifdef VCIEM
            , int *cp, double *vcte, int *psel1, int *psel2
	    #endif
          #endif
        )
//
// extract statistics from particles
//
// returns in nc_ the number of velocity-conditioned 2d elements
//
// returns the number of empty elements (containing no particles)
//
{
  int myid, e, en, eN, eU, e1N, p, pN, pU, p3, p4, i, n, ne, nc, A, B, C, NA, NB, NC;
  double dtmp;
  double nr[NDIMN];
  #ifdef MICROMIXING
  int e3, e2N;
  #endif
  #ifndef WALLFUNCTIONS
  int UA, UB, UC;
  double w1;
  double tke[NNODE];
  #endif


  ne = nc = 0;
  #ifdef _OPENMP
    #ifdef MICROMIXING
      #pragma omp parallel private(myid,p,p3,e,en,eN,e1N,e2N,eU,e3,dtmp,i) reduction(+:ne,nc)
    #else
      #pragma omp parallel private(myid,p,p3,e,en,eN,e1N,eU,dtmp) reduction(+:ne,nc)
    #endif // MICROMIXING
  #endif // _OPENMP
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif


    //////////////////////////// FIRST ORDER STATISTICS ////////////////////////////
    // each thread zeros its own portion of temporary array
    memset( _v1e + myid*nelem*5, 0, nelem*5*sizeof(double) );

    // compute first order statistics into temporary array
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( p = 0; p < npar; p++ )		// loop over all particles
    {
      p3 = p*3;
      eN = (myid*nelem + elp[p])*5;
      _v1e[eN+0] += parvel[p3+0];		// mean x velocity
      _v1e[eN+1] += parvel[p3+1];		// mean y velocity
      _v1e[eN+2] += parfreq[p];			// mean turbulent frequency
      #ifdef MICROMIXING
      _v1e[eN+3] += parc[p];			// unconditioned mean of scalar
      _v1e[eN+4] += partm[p];			// mean micromixing timescale
      #endif
    }

    // collect first order statistics from temporary array from all threads
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( e = 0; e < nelem; e++ )		// loop over all elements
      if (npel[e]==0) ne++;			// increase number of empty elements
      else
      {
        en = e*NDIMN;
        ue[en+0]=ue[en+1]=fe[e]=0.0;
        #ifdef MICROMIXING
	ce[e]=tme[e]=0.0;
	#endif
        for ( p = 0; p < nthreads; p++ )	// collect data from all threads
        {
          eN = (p*nelem + e)*5;
          ue[en+0] += _v1e[eN+0];
          ue[en+1] += _v1e[eN+1];
          fe[e] += _v1e[eN+2];
          #ifdef MICROMIXING
          ce[e] += _v1e[eN+3];
          tme[e] += _v1e[eN+4];
	  #endif
        }
        ue[en+0] /= npel[e];
        ue[en+1] /= npel[e];
        fe[e] /= npel[e];
        #ifdef MICROMIXING
        ce[e] /= npel[e];
        tme[e] /= npel[e];
	#endif
      }


    //////////////////////////// HIGHER ORDER STATISTICS ////////////////////////////
    // each thread zeros its own portion of temporary arrays
    memset( _v1e + myid*nelem*10, 0, nelem*10*sizeof(double) );
    #ifdef MICROMIXING
    memset( _v2e + myid*nelem*6, 0, nelem*6*sizeof(double) );
    #endif

    // compute higher order statistics into temporary arrays
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( p = 0; p < npar; p++ )					// loop over all particles
    {
      p3 = p*3;
      e = elp[p];
      en = e*NDIMN;

      e1N = (myid*nelem + e)*10;
      _v1e[e1N+0] += (parvel[p3+0]-ue[en+0])*(parvel[p3+0]-ue[en+0]);	// <uu> Reynolds stress
      _v1e[e1N+1] += (parvel[p3+1]-ue[en+1])*(parvel[p3+1]-ue[en+1]);	// <vv>
      _v1e[e1N+2] += parvel[p3+2]*parvel[p3+2];				// <ww>
      _v1e[e1N+3] += (parvel[p3+0]-ue[en+0])*(parvel[p3+1]-ue[en+1]);	// <uv>
      _v1e[e1N+4] += (parvel[p3+0]-ue[en+0])*parvel[p3+2];		// <uw>
      _v1e[e1N+5] += (parvel[p3+1]-ue[en+1])*parvel[p3+2];		// <vw>
      dtmp = (parvel[p3+0]-ue[en+0]);
      _v1e[e1N+6] += dtmp*dtmp*dtmp;					// skewness u
      _v1e[e1N+7] += dtmp*dtmp*dtmp*dtmp;				// flatness u
      dtmp = (parvel[p3+1]-ue[en+1]);
      _v1e[e1N+8] += dtmp*dtmp*dtmp;					// skewness v
      _v1e[e1N+9] += dtmp*dtmp*dtmp*dtmp;				// flatness v

      #ifdef MICROMIXING
      e2N = (myid*nelem + e)*6;
      dtmp = parc[p]-ce[e];						// scalar fluctuation
      _v2e[e2N+0] += dtmp*dtmp;						// scalar variance
      _v2e[e2N+1] += dtmp*dtmp*dtmp;					// scalar skewness
      _v2e[e2N+2] += dtmp*dtmp*dtmp*dtmp;				// scalar kurtosis
      _v2e[e2N+3] += dtmp*(parvel[p3+0]-ue[en+0]);			// velocity-scalar correlations
      _v2e[e2N+4] += dtmp*(parvel[p3+1]-ue[en+1]);
      _v2e[e2N+5] += dtmp*parvel[p3+2];
      #endif
    }

    // collect higher order statistics from temporary arrays from all threads
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( e = 0; e < nelem; e++ )		// loop over all elements
      if (npel[e])
      {
        eU = e*U2DOF;
        eN = e*NDIMN;
        u2e[eU+0]=u2e[eU+1]=u2e[eU+2]=u2e[eU+3]=u2e[eU+4]=u2e[eU+5]=
        u3e[eN+0]=u3e[eN+1]=u4e[eN+0]=u4e[eN+1]=0.0;
        #ifdef MICROMIXING
        e3 = e*3;
        c2e[e]=c3e[e]=c4e[e]=uce[e3+0]=uce[e3+1]=uce[e3+2]=0.0;
	#endif
        for ( p = 0; p < nthreads; p++ )	// collect data from all threads
        {
          e1N = (p*nelem + e)*10;
          u2e[eU+0] += _v1e[e1N+0];
          u2e[eU+1] += _v1e[e1N+1];
          u2e[eU+2] += _v1e[e1N+2];
          u2e[eU+3] += _v1e[e1N+3];
          u2e[eU+4] += _v1e[e1N+4];
          u2e[eU+5] += _v1e[e1N+5];
          u3e[eN+0] += _v1e[e1N+6];
          u4e[eN+0] += _v1e[e1N+7];
          u3e[eN+1] += _v1e[e1N+8];
          u4e[eN+1] += _v1e[e1N+9];
          #ifdef MICROMIXING
          e2N = (p*nelem + e)*6;
          c2e[e] += _v2e[e2N+0];
          c3e[e] += _v2e[e2N+1];
          c4e[e] += _v2e[e2N+2];
          uce[e3+0] += _v2e[e2N+3];
          uce[e3+1] += _v2e[e2N+4];
          uce[e3+2] += _v2e[e2N+5];
	  #endif
        }
	u2e[eU+0] /= npel[e];
	u2e[eU+1] /= npel[e];
	u2e[eU+2] /= npel[e];
	u2e[eU+3] /= npel[e];
	u2e[eU+4] /= npel[e];
	u2e[eU+5] /= npel[e];
	u3e[eN+0] /= npel[e];
	u3e[eN+1] /= npel[e];
	u4e[eN+0] /= npel[e];
	u4e[eN+1] /= npel[e];
        // finish calculation of skewness and flatness
        if (u2e[eU+0] > EPS)
        {
          u3e[eN+0] /= sqrt(u2e[eU+0]*u2e[eU+0]*u2e[eU+0]);	// skewness u
          u3e[eN+1] /= sqrt(u2e[eU+1]*u2e[eU+1]*u2e[eU+1]);	// skewness v
          u4e[eN+0] /= u2e[eU+0]*u2e[eU+0];			// flatness u
          u4e[eN+1] /= u2e[eU+1]*u2e[eU+1];			// flatness u
        }

	#ifdef MICROMIXING
        c2e[e] /= npel[e];
        c3e[e] /= npel[e];
        c4e[e] /= npel[e];
        uce[e3+0] /= npel[e];
        uce[e3+1] /= npel[e];
        uce[e3+2] /= npel[e];
	// compute vector of correlation coefficients
        if (c2e[e] > 1.0e-4)
        {
          if (u2e[eU+0] > EPS) uce[e3+0] /= sqrt(u2e[eU+0]*c2e[e]);
          if (u2e[eU+1] > EPS) uce[e3+1] /= sqrt(u2e[eU+1]*c2e[e]);
          if (u2e[eU+2] > EPS) uce[e3+2] /= sqrt(u2e[eU+2]*c2e[e]);
        }
        // finish calculation of skewness and kurtosis
        if (c2e[e] > EPS)
        {
          c3e[e] /= sqrt(c2e[e]*c2e[e]*c2e[e]);		// skewness
          c4e[e] /= c2e[e]*c2e[e];			// kurtosis
        }
	#endif

      }


    //////////////////////////// VELOCITY CONDITIONED MEAN ////////////////////////////
    #ifdef MICROMIXING
    #ifdef VCIEM
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for ( e = 0; e < nelem; e++ )		// loop over all elements
    {
      #ifdef PROJECTION	// use Fox's projection method to compute the velocity-conditioned mean
      e3 = e*3;
      // compute conditioned mean of scalar if not degenerate projection
      if ((fabs(uce[e3+0])>EPS) || (fabs(uce[e3+1])>EPS) || (fabs(uce[e3+2])>EPS))
      {
         nc++;		// in nc we add up the number of conditioned elements
         compute_conditioned_mean_proj( e, npel, cp, u2e, uce, parvel, parc, vcte, psel1, psel2 );
      }
      else
      // degenerate projection: velocity and scalar uncorrelated (locally homogeneous
      // scalar with no scalar gradient or locally laminar flow), revert to unconditioned
      // mean using the first bin (ie. IEM)
      {
         vcte[e*CNBI+0] = ce[e];			// put unconditioned mean into first bin
         for ( i = psel2[e]+1; i <= psel2[e+1]; i++ )	// loop over particles in element e
           cp[psel1[i]] = 0;				// store 0 as conditioning pointer
      }
      #else	// don't use Fox's projection method
      		// in nc we compute the average number of conditioning bins used
      nc += compute_conditioned_mean( e, npel, cp, parvel, parc, vcte, psel1, psel2 );
      #endif	// PROJECTION
    }
    #endif	// VCIEM
    #endif	// MICROMIXING
  }
  (*nc_) = nc;



  // save the previous timestep's (<Ui> ni) on the wall-boundary (needed for computing its time-derivative later)
  for ( nwe=i=0; i < nbpoin; i++ )	// loop over all boundary elements
    if (betags[i*2+0] == 4) // wall
    {
       pN = NDIMN*binpoel[i*NBNODE+0];	// get the index of one of the nodes of boundary element i
       nr[0] = wenr[nwe*NDIMN+0];	// get normal of wall-element nwe
       nr[1] = wenr[nwe*NDIMN+1];

       oudt[nwe] = u[pN+0]*nr[0] + u[pN+1]*nr[1];

       nwe++;				// increase wall-element counter
    }


  // transfer statistics from elements to points
  #ifdef _OPENMP
    #ifdef MICROMIXING
      #pragma omp parallel for private(p,pN,pU,p3,n,i,e,eN,eU,e3)
    #else
      #pragma omp parallel for private(p,pN,pU,p3,n,i,e,eN,eU)
    #endif // MICROMIXING
  #endif // _OPENMP
  for ( p = 0; p < npoin; p++ )
  {
     pN = p*NDIMN;
     pU = p*U2DOF;
     p3 = p*3;

     u[pN+0]=u[pN+1]=u2[pU+0]=u2[pU+1]=u2[pU+2]=u2[pU+3]=u2[pU+4]=u2[pU+5]=
     f[p]=u3[pN+0]=u3[pN+1]=u4[pN+0]=u4[pN+1]=0.0;
     #ifdef MICROMIXING
     c[p]=c2[p]=c3[p]=c4[p]=uc[p3+0]=uc[p3+1]=uc[p3+2]=tm[p]=0.0;
     #endif
     for ( n=0, i = esup2[p]+1; i <= esup2[p+1]; i++, n++ )
     {
        e = esup1[i];
        eN = e*NDIMN;
        eU = e*U2DOF;
	#ifdef MICROMIXING
	e3 = e*3;
	#endif
  
        u[pN+0] += ue[eN+0];
	u[pN+1] += ue[eN+1];
        u2[pU+0] += u2e[eU+0];
        u2[pU+1] += u2e[eU+1];
        u2[pU+2] += u2e[eU+2];
        u2[pU+3] += u2e[eU+3];
        u2[pU+4] += u2e[eU+4];
        u2[pU+5] += u2e[eU+5];
        u3[pN+0] += u3e[eN+0];
        u3[pN+1] += u3e[eN+1];
        u4[pN+0] += u4e[eN+0];
        u4[pN+1] += u4e[eN+1];
        f[p] += fe[e];
        #ifdef MICROMIXING
        c[p] += ce[e];
        c2[p] += c2e[e];
        c3[p] += c3e[e];
        c4[p] += c4e[e];
        uc[p3+0] += uce[e3+0];
        uc[p3+1] += uce[e3+1];
        uc[p3+2] += uce[e3+2];
        tm[p] += tme[e];
	#endif
     }
     u[pN+0] /= n;
     u[pN+1] /= n;
     u2[pU+0] /= n;
     u2[pU+1] /= n;
     u2[pU+2] /= n;
     u2[pU+3] /= n;
     u2[pU+4] /= n;
     u2[pU+5] /= n;
     u3[pN+0] /= n;
     u3[pN+1] /= n;
     u4[pN+0] /= n;
     u4[pN+1] /= n;
     f[p] /= n;
     #ifdef MICROMIXING
     c[p] /= n;
     c2[p] /= n;
     c3[p] /= n;
     c4[p] /= n;
     uc[p3+0] /= n;
     uc[p3+1] /= n;
     uc[p3+2] /= n;
     tm[p] /= n;
     #endif
  }


  // enforce boundary conditions on statistics in points at wall
  #ifndef WALLFUNCTIONS
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( i = 0; i < nbpoin; i++ )			// loop over all boundary points
     if (bptags[i*3+1] == 4)				// wall
       u[bpg[i]*NDIMN+0] = u[bpg[i]*NDIMN+1] =				// <Ui> = 0
       u2[bpg[i]*U2DOF+0] = u2[bpg[i]*U2DOF+1] = u2[bpg[i]*U2DOF+2] = 	// <uiuj> = 0
       u2[bpg[i]*U2DOF+3] = u2[bpg[i]*U2DOF+4] = u2[bpg[i]*U2DOF+5] = 0.0;
  #endif


  // set boundary condition on mean turbulent frequency at wall
  #ifndef WALLFUNCTIONS
  for ( e = 0; e < nwe; e++ )           // loop over all wall elements
  {
     eN = we[e]*NNODE;                   // get element node information for domain-element we[e]
     A = inpoel[eN+0];
     B = inpoel[eN+1];
     C = inpoel[eN+2];
     UA = A*U2DOF;
     UB = B*U2DOF;
     UC = C*U2DOF;
     // compute sqrt(2k) in wall element
     tke[0] = sqrt(u2[UA+0] + u2[UA+1] + u2[UA+2]);
     tke[1] = sqrt(u2[UB+0] + u2[UB+1] + u2[UB+2]);
     tke[2] = sqrt(u2[UC+0] + u2[UC+1] + u2[UC+2]);
     // w1 = d[sqrt(2k)]/dn / CT at the wall
     w1 = -((dNx[eN+0]*tke[0]+dNx[eN+1]*tke[1]+dNx[eN+2]*tke[2])*wenr[e*NDIMN+0] +
            (dNy[eN+0]*tke[0]+dNy[eN+1]*tke[1]+dNy[eN+2]*tke[2])*wenr[e*NDIMN+1])/CT;
     if ( w1 < EPS ) w1=BOUND;
     switch (weo[e])
     {
        case 0  : f[B] = f[C] = w1;  break;
        case 1  : f[A] = f[C] = w1;  break;
        case 2  : f[A] = f[B] = w1;  break;
     }
  }
  #endif


  // compute normal component of mean-velocity time derivative: (d<Ui>/dt ni) on the wall-boundary
  // (needed for wall-boundary condition on the mean pressure, in unsteady case)
  for ( nwe=i=0; i < nbpoin; i++ )	// loop over all boundary elements
    if (betags[i*2+0] == 4) // wall
    {
       pN = NDIMN*binpoel[i*NBNODE+0];	// get the index of one of the nodes of boundary element i
       nr[0] = wenr[nwe*NDIMN+0];	// get normal of wall-element nwe
       nr[1] = wenr[nwe*NDIMN+1];

       oudt[nwe] = (u[pN+0]*nr[0] + u[pN+1]*nr[1] - oudt[nwe]) / dt;

       nwe++;				// increase wall-element counter
    }


  // spatial derivatives
  // 1st pass: first derivative, element based
  #ifdef _OPENMP
  #pragma omp parallel for private(e,eN,A,B,C,NA,NB,NC)
  #endif
  for ( e = 0; e < nelem; e++ )
  {
     eN = e*NNODE;
     A = inpoel[eN+0];   // get node-information for element e
     B = inpoel[eN+1];
     C = inpoel[eN+2];
     NA = A*NDIMN;
     NB = B*NDIMN;
     NC = C*NDIMN;
     
     // d<U>/dx
     _v1e[e*4+0] = dNx[eN+0]*u[NA+0] + dNx[eN+1]*u[NB+0] + dNx[eN+2]*u[NC+0];
     // d<U>/dy
     _v1e[e*4+1] = dNy[eN+0]*u[NA+0] + dNy[eN+1]*u[NB+0] + dNy[eN+2]*u[NC+0];
     // d<V>/dx
     _v1e[e*4+2] = dNx[eN+0]*u[NA+1] + dNx[eN+1]*u[NB+1] + dNx[eN+2]*u[NC+1];
     // d<V>/dy
     _v1e[e*4+3] = dNy[eN+0]*u[NA+1] + dNy[eN+1]*u[NB+1] + dNy[eN+2]*u[NC+1];
  }


  // 2nd pass: transfer first derivative from elements to points
  #ifdef _OPENMP
  #pragma omp parallel for private(p,p4,n,i,e)
  #endif
  for ( p = 0; p < npoin; p++ )
  {
     p4 = p*4;
     du[p4+0]=du[p4+1]=du[p4+2]=du[p4+3]=0.0;
     for ( n=0, i = esup2[p]+1; i <= esup2[p+1]; i++, n++ )
     {
        e = esup1[i];
        du[p4+0] += _v1e[e*4+0];
        du[p4+1] += _v1e[e*4+1];
        du[p4+2] += _v1e[e*4+2];
        du[p4+3] += _v1e[e*4+3];
     }
     du[p4+0] /= n;
     du[p4+1] /= n;
     du[p4+2] /= n;
     du[p4+3] /= n;
  }


  // 3rd pass: Laplacian of mean velocity, element based
  #ifdef _OPENMP
  #pragma omp parallel for private(e,eN,A,B,C)
  #endif
  for ( e = 0; e < nelem; e++ )
  {
     eN = e*NNODE;
     A = inpoel[eN+0];   // get node-information for element e
     B = inpoel[eN+1];
     C = inpoel[eN+2];
     // d^2<U>/dx^2 + d^2<U>/dy^2
     _v1e[e*4+0] = dNx[eN+0]*du[A*4+0] + dNx[eN+1]*du[B*4+0] + dNx[eN+2]*du[C*4+0]
                  +dNy[eN+0]*du[A*4+1] + dNy[eN+1]*du[B*4+1] + dNy[eN+2]*du[C*4+1];
     // d^2<V>/dx^2 + d^2<V>/dy^2
     _v1e[e*4+1] = dNx[eN+0]*du[A*4+2] + dNx[eN+1]*du[B*4+2] + dNx[eN+2]*du[C*4+2]
                  +dNy[eN+0]*du[A*4+3] + dNy[eN+1]*du[B*4+3] + dNy[eN+2]*du[C*4+3];
  }

  // 4th pass: transfer Laplacian from elements to points
  #ifdef _OPENMP
  #pragma omp parallel for private(p,pN,n,i,eN)
  #endif
  for ( p = 0; p < npoin; p++ )
  {
     pN = p*NDIMN;
     ddu[pN+0]=ddu[pN+1]=0.0;
     for ( n=0, i = esup2[p]+1; i <= esup2[p+1]; i++, n++ )
     {
        eN = esup1[i]*4;
        ddu[pN+0] += _v1e[eN+0];
        ddu[pN+1] += _v1e[eN+1];
     }
     ddu[pN+0] /= n;
     ddu[pN+1] /= n;
  }

  return( ne );	// return the number of empty elements encountered
}
