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
//  externally callable functions defined in particles.cc
//

int minnpel( int nelem, int *npel );

int findpar( int p, int e, int myid, int nelem, double minx, double maxx, int *npeldt, int *elp,
             int *inpoel, int *esuel, int *esupel1, int *esupel2, double *coord, double *parcoord,
             double *Ae, int *bf );

void improvepar( int nelem, int npar, int nthreads, int *cnt, int *npeldt, int *elp, int *npel,
                 double *parcoord, double *parvel, int *psel1, int *psel2, double *parfreq
                 #ifdef MICROMIXING
	         , double *parc
	         #endif
	       );
void improvepar2( int nelem, int npar, int nthreads, int *cnt, int *npeldt, int *elp, int *npel,
                  double *parcoord, double *parvel, int *psel1, int *psel2, double *parfreq
                  #ifdef MICROMIXING
	          , double *parc
	          #endif
	        );

void genpsel( int nelem, int npar, int nthreads, int *npeldt, int *npel, int *elp, int *psel1, int *psel2 );

void saverestartpar( int nthreads, int npoin, int npar, int onum, int it, double t, double *parcoord, double *parvel,
                     double *tu, double *tdu, double *parfreq, double *tu2, double *tu3, double *tu4, double *tf,
                     double *tpr
	             #ifndef WALLFUNCTIONS
                     , VSLStreamStatePtr *stream
		     #endif
		     #ifdef MICROMIXING
		     , double *parc, double *tc, double *tc2, double *tc3, double *tc4, double *tuc
		     #endif
		   );

void initgenpar( int nelem, int nthreads, int *inpoel, int *esuel, int *esupel1, int *esupel2,
                 double minx, double maxx, double *Ae, double *coord, VSLStreamStatePtr *stream, 
		 int *npar, int *onum, int **elp, int **npel, int **npeldt,
                 double **parcoord, double **parvel,
		 int npoin, double *tu, double *tdu, double *tpr, int *it, int *it0, double *t, double *t0,
                 int **psel1, int **psel2, double **parfreq, double *tu2, double *tu3, double *tu4, double *tf
		 #ifdef MICROMIXING
		 , double **parc, double **partm,
		 double *tc, double *tc2, double *tc3, double *tc4, double *tuc
		 #ifdef VCIEM
		 , int **cp
		 #endif
		 #endif
               );

void genes( int ndl, int nelem, double miny, double maxy, int *inpoel, double *coord, int **es1, int **es2 );
void genec( int nelem, double maxx, int *inpoel, double *coord, int **ec, int *size );

void destroy_tmp_par( void );

int stat( int nelem, int npoin, int nbpoin, int npar, int nthreads, double dt,
          int *npel, int *esup1, int *esup2, int *inpoel, int *elp,
          int nwe, double *wenr, double *ue, double *u, double *du, double *ddu, double *parvel,
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
        );





#define PARINEL(p,e,a1,a2,coord,inpoel,parcoord,Ae)					\
 (((a1=(coord[inpoel[e*NNODE+2]*NDIMN+0]-coord[inpoel[e*NNODE+1]*NDIMN+0])*		\
       (parcoord[p*2+1]-coord[inpoel[e*NNODE+1]*NDIMN+1])-				\
       (parcoord[p*2+0]-coord[inpoel[e*NNODE+1]*NDIMN+0])*				\
       (coord[inpoel[e*NNODE+2]*NDIMN+1]-coord[inpoel[e*NNODE+1]*NDIMN+1]))>0) &&	\
  ((a2=(parcoord[p*2+0]-coord[inpoel[e*NNODE+1]*NDIMN+0])*				\
       (coord[inpoel[e*NNODE+0]*NDIMN+1]-coord[inpoel[e*NNODE+1]*NDIMN+1])-		\
       (coord[inpoel[e*NNODE+0]*NDIMN+0]-coord[inpoel[e*NNODE+1]*NDIMN+0])*		\
       (parcoord[p*2+1]-coord[inpoel[e*NNODE+1]*NDIMN+1]))>0) && (0.5*(a1+a2)<Ae[e]))
/*
 * The above PARINEL macro as a function:
 * particle search based on vector products
 * cheapest I found, doesn't give a directional hint if particle not found
 *
int parinel( int p, int e )
//
// returns 1 if particle p is in 2d element e, 0 if not
//
{
  int A, B, C;
  double NA, NB, xp, yp;
  double x[NNODE], y[NNODE];


  // get node-, and particle position-information
  A = inpoel[e*NNODE+0];
  B = inpoel[e*NNODE+1];
  C = inpoel[e*NNODE+2];
  x[0] = coord[A*NDIMN+0];  y[0] = coord[A*NDIMN+1];
  x[1] = coord[B*NDIMN+0];  y[1] = coord[B*NDIMN+1];
  x[2] = coord[C*NDIMN+0];  y[2] = coord[C*NDIMN+1];
  xp = parcoord[p*2+0];  yp = parcoord[p*2+1];
  
  // evaluate shapefunctions at particle location (Cramer's rule)
  NA = (x[2]-x[1])*(yp-y[1]) - (xp-x[1])*(y[2]-y[1]);
  NB = (xp-x[1])*(y[0]-y[1]) - (x[0]-x[1])*(yp-y[1]);

  if ( (NA>0) && (NB>0) && (0.5*(NA+NB)<Ae[e]) )
    return(1);
  else
    return(0);
}
*/




