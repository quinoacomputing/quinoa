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
//  Postprocessing functions from/to gmsh format
//  for more info see main.cc
//


#include <stdio.h>
#include <math.h>
#include "macros.h"
#include "const.h"
#include "matrix3.h"
#include "postprocess.h"








void outparpos( int npar, double t, double *parcoord )
//
// output particle positions as a new view into OUT_POSTPROCESS_PARTICLES_FILENAME
// using t as the viewname
//
{
  int p;
  FILE *opos;
  

  printf("outparpos: %.3g\n", t);

  // erase OUT_POSTPROCESS_PARTICLES_FILENAME if this is the first timestep...
  if ( t < EPS )
  {
    opos = fopen( OUT_POSTPROCESS_PARTICLES_FILENAME, "w" );
    fclose( opos );
  }

  if ( !(opos = fopen(OUT_POSTPROCESS_PARTICLES_FILENAME,"a")) )
    ERR("Cannot open OUT_POSTPROCESS_PARTICLES_FILENAME");
  
  // write out particle positions...
  fprintf( opos, "\nView \"%.3g\" {", t );
  for ( p = 0; p < npar; p++ )
    fprintf( opos, "SP(%g,%g,%g) {%d};\n", parcoord[p*2+0], parcoord[p*2+1], 0.0, 1 );
  fprintf( opos, "};\n\n" );

  fclose( opos );
}






void out_3dstep( int nelem, int npoin, int *inpoel, double *coord, double *u,
                 double *pr, double *u2, double *u3, double *u4, double *f, double *du,
		 /*#ifndef WALLFUNCTIONS
		   double *rho,
		 #endif*/
		 #ifdef MICROMIXING
		 double *c, double *c2, double *c3, double *c4, double *tm,
		 #endif
		 const char *filename, int av )
//
// outputs results
//
{
  int e, eN, p, A, B, C, NA, NB, NC, UA, UB, UC;
  FILE *opos;
  double tke[NNODE], eps[NNODE];
  /*#ifndef WALLFUNCTIONS
  int RA, RB, RC;
  #endif*/

  
  if ( av )
    printf("out_3dtavstep()...\n");
  else
    printf("out_3dstep()...\n");
  
  if (!(opos = fopen(filename,"w"))) ERR("Cannot open postprocess file");

  // write out x mean velocity, <U>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_MEANVELOCITY_U );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u[NA+0],
	   coord[NB+0], coord[NB+1], u[NB+0],
	   coord[NC+0], coord[NC+1], u[NC+0],
	   u[NA+0], u[NB+0], u[NC+0] );
  }
  fprintf( opos, "};\n\n" );


  // write out y mean velocity, <V>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_MEANVELOCITY_V );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u[NA+1],
	   coord[NB+0], coord[NB+1], u[NB+1],
	   coord[NC+0], coord[NC+1], u[NC+1],
	   u[NA+1], u[NB+1], u[NC+1] );
  }
  fprintf( opos, "};\n\n" );


  // write out Reynolds stress component <uu>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_REYNOLDS11 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u2[UA+0],
	   coord[NB+0], coord[NB+1], u2[UB+0],
	   coord[NC+0], coord[NC+1], u2[UC+0],
	   u2[UA+0], u2[UB+0], u2[UC+0] );
  }
  fprintf( opos, "};\n\n" );


  // write out Reynolds stress component <vv>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_REYNOLDS22 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u2[UA+1],
	   coord[NB+0], coord[NB+1], u2[UB+1],
	   coord[NC+0], coord[NC+1], u2[UC+1],
	   u2[UA+1], u2[UB+1], u2[UC+1] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out Reynolds stress component <ww>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_REYNOLDS33 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;

    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u2[UA+2],
	   coord[NB+0], coord[NB+1], u2[UB+2],
	   coord[NC+0], coord[NC+1], u2[UC+2],
	   u2[UA+2], u2[UB+2], u2[UC+2] );
  }
  fprintf( opos, "};\n\n" );


  // write out Reynolds stress component <uv>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_REYNOLDS12 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u2[UA+3],
	   coord[NB+0], coord[NB+1], u2[UB+3],
	   coord[NC+0], coord[NC+1], u2[UC+3],
	   u2[UA+3], u2[UB+3], u2[UC+3] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out Reynolds stress component <uw>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_REYNOLDS13 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u2[UA+4],
	   coord[NB+0], coord[NB+1], u2[UB+4],
	   coord[NC+0], coord[NC+1], u2[UC+4],
	   u2[UA+4], u2[UB+4], u2[UC+4] );
  }
  fprintf( opos, "};\n\n" );


  // write out Reynolds stress component <vw>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_REYNOLDS23 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u2[UA+5],
	   coord[NB+0], coord[NB+1], u2[UB+5],
	   coord[NC+0], coord[NC+1], u2[UC+5],
	   u2[UA+5], u2[UB+5], u2[UC+5] );
  }
  fprintf( opos, "};\n\n" );

 
  // write out skewness of streamwise velocity...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_SKEWNESS11 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;

    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u3[NA+0],
	   coord[NB+0], coord[NB+1], u3[NB+0],
	   coord[NC+0], coord[NC+1], u3[NC+0],
	   u3[NA+0], u3[NB+0], u3[NC+0] );
  }
  fprintf( opos, "};\n\n" );


  // write out skewness of cross-stream velocity...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_SKEWNESS22 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;

    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u3[NA+1],
	   coord[NB+0], coord[NB+1], u3[NB+1],
	   coord[NC+0], coord[NC+1], u3[NC+1],
	   u3[NA+1], u3[NB+1], u3[NC+1] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out flatness of streamwise velocity...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_FLATNESS11 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;

    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u4[NA+0],
	   coord[NB+0], coord[NB+1], u4[NB+0],
	   coord[NC+0], coord[NC+1], u4[NC+0],
	   u4[NA+0], u4[NB+0], u4[NC+0] );
  }
  fprintf( opos, "};\n\n" );


  // write out flatness of cross-stream velocity...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_FLATNESS22 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;

    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], u4[NA+1],
	   coord[NB+0], coord[NB+1], u4[NB+1],
	   coord[NC+0], coord[NC+1], u4[NC+1],
	   u4[NA+1], u4[NB+1], u4[NC+1] );
  }
  fprintf( opos, "};\n\n" );


  // write out mean pressure, <p>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_MEAN_PRESSURE );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], pr[A],
	   coord[NB+0], coord[NB+1], pr[B],
	   coord[NC+0], coord[NC+1], pr[C],
	   pr[A], pr[B], pr[C] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out mean-velocity as vector-field...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_MEANVELOCITY );
  for ( p = 0; p < npoin; p++ )
    fprintf( opos, "VP(%g,%g,0){%g,%g,0};\n",
	     coord[p*NDIMN+0], coord[p*NDIMN+1], u[p*NDIMN+0], u[p*NDIMN+1] );
  fprintf( opos, "};\n\n" );


  // write out turbulent kinetic energy...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_TKE );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;
    tke[0] = 0.5*(u2[UA+0] + u2[UA+1] + u2[UA+2]);
    tke[1] = 0.5*(u2[UB+0] + u2[UB+1] + u2[UB+2]);
    tke[2] = 0.5*(u2[UC+0] + u2[UC+1] + u2[UC+2]);
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], tke[0],
	   coord[NB+0], coord[NB+1], tke[1],
	   coord[NC+0], coord[NC+1], tke[2],
	   tke[0], tke[1], tke[2] );
  }
  fprintf( opos, "};\n\n" );


  // write out rate of dissipation of turbulent kinetic energy...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_DISSIPATION );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    UA = A*U2DOF;
    UB = B*U2DOF;
    UC = C*U2DOF;
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    tke[0] = 0.5*(u2[UA+0] + u2[UA+1] + u2[UA+2]);
    tke[1] = 0.5*(u2[UB+0] + u2[UB+1] + u2[UB+2]);
    tke[2] = 0.5*(u2[UC+0] + u2[UC+1] + u2[UC+2]);
    #ifndef WALLFUNCTIONS
    eps[0] = f[A]*(tke[0] + CT2*f[A]/RE);
    eps[1] = f[B]*(tke[1] + CT2*f[B]/RE);
    eps[2] = f[C]*(tke[2] + CT2*f[C]/RE);
    #else
    eps[0] = f[A]*tke[0];
    eps[1] = f[B]*tke[1];
    eps[2] = f[C]*tke[2];
    #endif
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], eps[0],
	   coord[NB+0], coord[NB+1], eps[1],
	   coord[NC+0], coord[NC+1], eps[2],
	   eps[0], eps[1], eps[2] );
  }
  fprintf( opos, "};\n\n" );


  // write out mean spanwise vorticity...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_VORTICITY );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    UA = A*4;
    UB = B*4;
    UC = C*4;
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;

    // spanwise vorticity: d<V>/dx - d<U>/dy
    tke[0] = 0.5*(du[UA+2] - du[UA+1]);
    tke[1] = 0.5*(du[UB+2] - du[UB+1]);
    tke[2] = 0.5*(du[UC+2] - du[UC+1]);
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], tke[0],
	   coord[NB+0], coord[NB+1], tke[1],
	   coord[NC+0], coord[NC+1], tke[2],
	   tke[0], tke[1], tke[2] );
  }
  fprintf( opos, "};\n\n" );


  /*#ifndef WALLFUNCTIONS
  // write out wiggly-p11...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P11 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+0],
	   coord[NB+0], coord[NB+1], rho[RB+0],
	   coord[NC+0], coord[NC+1], rho[RC+0],
	   rho[RA+0], rho[RB+0], rho[RC+0] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out wiggly-p22...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P22 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+4],
	   coord[NB+0], coord[NB+1], rho[RB+4],
	   coord[NC+0], coord[NC+1], rho[RC+4],
	   rho[RA+4], rho[RB+4], rho[RC+4] );
  }
  fprintf( opos, "};\n\n" );


  // write out wiggly-p33...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P33 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+8],
	   coord[NB+0], coord[NB+1], rho[RB+8],
	   coord[NC+0], coord[NC+1], rho[RC+8],
	   rho[RA+8], rho[RB+8], rho[RC+8] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out wiggly-p12...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P12 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+1],
	   coord[NB+0], coord[NB+1], rho[RB+1],
	   coord[NC+0], coord[NC+1], rho[RC+1],
	   rho[RA+1], rho[RB+1], rho[RC+1] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out wiggly-p13...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P13 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+2],
	   coord[NB+0], coord[NB+1], rho[RB+2],
	   coord[NC+0], coord[NC+1], rho[RC+2],
	   rho[RA+2], rho[RB+2], rho[RC+2] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out wiggly-p23...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P23 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+5],
	   coord[NB+0], coord[NB+1], rho[RB+5],
	   coord[NC+0], coord[NC+1], rho[RC+5],
	   rho[RA+5], rho[RB+5], rho[RC+5] );
  }
  fprintf( opos, "};\n\n" );

  
  // write out wiggly-p21...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P21 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+3],
	   coord[NB+0], coord[NB+1], rho[RB+3],
	   coord[NC+0], coord[NC+1], rho[RC+3],
	   rho[RA+3], rho[RB+3], rho[RC+3] );
  }
  fprintf( opos, "};\n\n" );


  // write out wiggly-p31...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P31 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+6],
	   coord[NB+0], coord[NB+1], rho[RB+6],
	   coord[NC+0], coord[NC+1], rho[RC+6],
	   rho[RA+6], rho[RB+6], rho[RC+6] );
  }
  fprintf( opos, "};\n\n" );


  // write out wiggly-p32...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_WIGGLY_P32 );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    RA = A*RDOF;
    RB = B*RDOF;
    RC = C*RDOF;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], rho[RA+7],
	   coord[NB+0], coord[NB+1], rho[RB+7],
	   coord[NC+0], coord[NC+1], rho[RC+7],
	   rho[RA+7], rho[RB+7], rho[RC+7] );
  }
  fprintf( opos, "};\n\n" );
  #endif // WALLFUNCTIONS*/


  #ifdef MICROMIXING
  // write out mean scalar, <C>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_MEAN_SCALAR );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], c[A],
	   coord[NB+0], coord[NB+1], c[B],
	   coord[NC+0], coord[NC+1], c[C],
	   c[A], c[B], c[C] );
  }
  fprintf( opos, "};\n\n" );


  // write out variance of scalar, <c2>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_VARIANCE_SCALAR );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], c2[A],
	   coord[NB+0], coord[NB+1], c2[B],
	   coord[NC+0], coord[NC+1], c2[C],
	   c2[A], c2[B], c2[C] );
  }
  fprintf( opos, "};\n\n" );


  // write out skewness of scalar, <c3>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_SKEWNESS_SCALAR );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], c3[A],
	   coord[NB+0], coord[NB+1], c3[B],
	   coord[NC+0], coord[NC+1], c3[C],
	   c3[A], c3[B], c3[C] );
  }
  fprintf( opos, "};\n\n" );


  // write out kurtosis of scalar, <c4>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_KURTOSIS_SCALAR );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], c4[A],
	   coord[NB+0], coord[NB+1], c4[B],
	   coord[NC+0], coord[NC+1], c4[C],
	   c4[A], c4[B], c4[C] );
  }
  fprintf( opos, "};\n\n" );


  // write out mean micrmomixing timescale, <tm>...
  fprintf( opos, "\nView \"%s\" {", VIEWNAME_MEAN_TM );
  for ( e = 0; e < nelem; e++ )
  {
    eN = e*NNODE;
    A = inpoel[eN+0];
    B = inpoel[eN+1];
    C = inpoel[eN+2];
    NA = A*NDIMN;
    NB = B*NDIMN;
    NC = C*NDIMN;
    
    fprintf( opos,
	   "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
	   coord[NA+0], coord[NA+1], tm[A],
	   coord[NB+0], coord[NB+1], tm[B],
	   coord[NC+0], coord[NC+1], tm[C],
	   tm[A], tm[B], tm[C] );
  }
  fprintf( opos, "};\n\n" );
  #endif // MICROMIXING

  fclose( opos );
}





#ifdef MICROMIXING
static void outpdf( int k, int tav, int onum, int *epdfloc, int *esupel1, int *esupel2, double *c2e,
                    int *psel1, int *psel2, double *parc, double *ce, int *sl, double *tpdf )
//
// outputs 'instantaneous' probability density function of temperature fluctuations at location k,
// depending on tav, also outputs time-averaged pdf
//
// k: location index,
// tav: 0 - output 'instantaneous' pdf only,
// 	1 - output both 'instantaneous' and time-averaged pdf.
//
{
  int i, j, l, n, e;
  char s[STRLEN];
  double stp;
  double pdf[NBI];
  FILE *fpdf, *tfpdf=0;

  
  e = epdfloc[k];

  for ( l = 0; l < NBI; l++ )
    pdf[l] = 0.0;

  // consider elements surrounding nodes of element e
  n = 0;
  for ( i = esupel2[e]+1; i <= esupel2[e+1]; i++ )// loop over elements surrounding nodes of element
    if ( c2e[esupel1[i]] > EPS )				    // only if variance != 0
      for ( j = psel2[esupel1[i]]+1; j <= psel2[esupel1[i]+1]; j++ )// loop over all particles in element esupel1[i]
      {
        n++;					// incrase number of particles considered
        stp = (parc[psel1[j]]-ce[esupel1[i]])/sqrt(c2e[esupel1[i]]);// calculate standardized temperature

        for ( l = 0; l < NBI; l++ )		// count up particles in bins
           if ( (stp > SMIN+l*BINSIZE) && (stp < SMIN+(l+1)*BINSIZE) ) pdf[l] += 1.0;
      }

  // compute pdf
  for ( l = 0; l < NBI; l++ )
    pdf[l] /= n*BINSIZE;

  
  if ( n )					// only if there is anything meaningful to save
  {
    // construct 'instantaneous' pdf filename and open file
    sprintf( s, "%s.%.4g_%.4g", PDF_BASE_FILENAME, PL[sl[k]*2+0], PL[sl[k]*2+1] );
    if ( !(fpdf = fopen(s,"w")) ) ERR("Cannot open file\n");
  
    if (tav)
    {
      // construct time-averaged pdf filename and open file
      sprintf( s, "%s.%.4g_%.4g", TPDF_BASE_FILENAME, PL[sl[k]*2+0], PL[sl[k]*2+1] );
      if ( !(tfpdf = fopen(s,"w")) ) ERR("Cannot open file\n");
    }

    for ( l = 0; l < NBI; l++ )
    {
      // pdf sample space variable for bin l
      stp = SMIN+(0.5+l)*BINSIZE;
      // output 'instantaneous' pdf
      fprintf( fpdf, "%g\t%g\n", stp, pdf[l] );
      
      if (tav)
      {
        // calculate cumulative time-averaged pdf of temperature
        if ( onum == 1 )
	{
	  tpdf[k*NBI+l] = pdf[l];
	}
	else
	{
          tpdf[k*NBI+l] = ((double)(onum-1)/onum)*(tpdf[k*NBI+l] + pdf[l]/(onum-1));
	}

        // output time-averaged pdf
        fprintf( tfpdf, "%g\t%g\n", stp, tpdf[k*NBI+l] );
      }
    }
  
    fclose(fpdf);
    if (tav)
    {
      fclose(tfpdf);
    }
  }
 
}
#endif







static void extract_element_profilepoint_values( int e, int i, int *n, int maxn, double maxy, double *yprof,
                                                 double *dete, int *inpoel, double *coord, double *u_,
                                                 double *u2_, double *u3_, double *u4_, double *f_,
                                                 double *uprof, double *vprof, double *u2prof, double *v2prof,
                                                 double *u3prof, double *v3prof, double *u4prof, double *v4prof,
                                                 double *uvprof, double *tkeprof, double *epsprof
                                                 #ifdef MICROMIXING
                                                 , double *cpeak, double *c2peak, double *c_, double *c2_,
                                                 double *c3_, double *c4_, double *tm_, double *cprof,
                                                 double *c2prof, double *c3prof, double *c4prof, double *tmprof
                                                 #endif
                                               )

//
// extract values at a profile-point along a transverse line from element e
//
// e : 2d element index,
// i : downstream location index,
// n : index of profile-point (gets increased here),
// cpeak : peak value of mean scalar,
// c2peak : peak value of scalar variance
//
{
  int p, A, B, C, AN, BN, CN, UA, UB, UC, in, needboth=0;
  double dtmp, y1, y2, NA, NB, NC;
  double x[NNODE], y[NNODE], tke[NNODE], eps[NNODE];


  A = inpoel[e*NNODE+0];		// get node information
  B = inpoel[e*NNODE+1];
  C = inpoel[e*NNODE+2];

  // order nodes according to non-decreasing x coordinate in index-order A, C, B
  if ( coord[B*NDIMN+0] < coord[A*NDIMN+0] ) { p=A; A=B; B=p; }
  if ( coord[C*NDIMN+0] < coord[A*NDIMN+0] ) { p=A; A=C; C=p; }
  if ( coord[B*NDIMN+0] < coord[C*NDIMN+0] ) { p=B; B=C; C=p; }

  // compute edge-intersections
  if ( (DL[i] > coord[A*NDIMN+0]) && (DL[i] < coord[B*NDIMN+0]) )	// if element bounding stripe
  {
     if ( DL[i] < coord[C*NDIMN+0] )	// AC edge
     {
       if ( fabs(dtmp=coord[C*NDIMN+0]-coord[A*NDIMN+0]) > EPS )	// if not parallel
         y1 = coord[A*NDIMN+1] + (DL[i]-coord[A*NDIMN+0])*(coord[C*NDIMN+1]-coord[A*NDIMN+1])/dtmp;
       else y1 = 0.5*(coord[A*NDIMN+1]+coord[C*NDIMN+1]);	// if parallel: average yA and yC
     }
     else				// CB edge
       if ( fabs(dtmp=coord[C*NDIMN+0]-coord[B*NDIMN+0]) > EPS )	// if not parallel
         y1 = coord[B*NDIMN+1] + (DL[i]-coord[B*NDIMN+0])*(coord[C*NDIMN+1]-coord[B*NDIMN+1])/dtmp;
       else y1 = 0.5*(coord[C*NDIMN+1]+coord[B*NDIMN+1]);	// if parallel: average yC and yB

     // AB edge
     y2 = coord[A*NDIMN+1] + (DL[i]-coord[A*NDIMN+0])*(coord[B*NDIMN+1]-coord[A*NDIMN+1])/
	  (coord[B*NDIMN+0]-coord[A*NDIMN+0]);

     // if an edge is aligned with the upper wall, count both intersection points
     if ( fabs(coord[A*NDIMN+1]-maxy)<EPS )
     {
       if ( (fabs(coord[B*NDIMN+1]-maxy)<EPS) || (fabs(coord[C*NDIMN+1]-maxy)<EPS) ) needboth = 1;
     }
     else
       if ( (fabs(coord[B*NDIMN+1]-maxy)<EPS) && (fabs(coord[C*NDIMN+1]-maxy)<EPS) ) needboth = 1;

     A = inpoel[e*NNODE+0];		// get node-info again
     B = inpoel[e*NNODE+1];
     C = inpoel[e*NNODE+2];
     UA = A*U2DOF;
     UB = B*U2DOF;
     UC = C*U2DOF;
     AN = A*NDIMN;
     BN = B*NDIMN;
     CN = C*NDIMN;
     x[0] = coord[A*NDIMN+0];  y[0] = coord[A*NDIMN+1];
     x[1] = coord[B*NDIMN+0];  y[1] = coord[B*NDIMN+1];
     x[2] = coord[C*NDIMN+0];  y[2] = coord[C*NDIMN+1];

     in = i*maxn+(*n);	// compute array index

     // save the lower one
     yprof[in] = fmin(y1,y2);
     // evaluate shapefunctions at particle location (Cramer's rule)
     NA = (DL[i]*(y[1]-y[2]) + x[1]*(y[2]-yprof[in]) + x[2]*(yprof[in]-y[1]))/dete[e];
     NB = (x[0]*(yprof[in]-y[2]) + DL[i]*(y[2]-y[0]) + x[2]*(y[0]-yprof[in]))/dete[e];
     NC = (x[0]*(y[1]-yprof[in]) + x[1]*(yprof[in]-y[0]) + DL[i]*(y[0]-y[1]))/dete[e];
     // interpolate nodal values
     #ifdef MICROMIXING
     cprof[in] = NA*c_[A] + NB*c_[B] + NC*c_[C];		// scalar mean
     c2prof[in] = NA*c2_[A] + NB*c2_[B] + NC*c2_[C];		// scalar variance
     c3prof[in] = NA*c3_[A] + NB*c3_[B] + NC*c3_[C];		// scalar skewness
     c4prof[in] = NA*c4_[A] + NB*c4_[B] + NC*c4_[C];		// scalar kurtosis
     if (cprof[in] > (*cpeak)) (*cpeak) = cprof[in];		// find peak of mean
     if (c2prof[in] > (*c2peak)) (*c2peak) = c2prof[in];	// find peak of variance
     tmprof[in] = NA*tm_[A] + NB*tm_[B] + NC*tm_[C];		// tm
     #endif
     tke[0] = 0.5*(u2_[UA+0] + u2_[UA+1] + u2_[UA+2]);
     tke[1] = 0.5*(u2_[UB+0] + u2_[UB+1] + u2_[UB+2]);
     tke[2] = 0.5*(u2_[UC+0] + u2_[UC+1] + u2_[UC+2]);
     tkeprof[in] = NA*tke[0] + NB*tke[1] + NC*tke[2];		// tke
     #ifndef WALLFUNCTIONS
     eps[0] = f_[A]*(tke[0] + CT2*f_[A]/RE);
     eps[1] = f_[B]*(tke[1] + CT2*f_[B]/RE);
     eps[2] = f_[C]*(tke[2] + CT2*f_[C]/RE);
     #else
     eps[0] = f_[A]*tke[0];
     eps[1] = f_[B]*tke[1];
     eps[2] = f_[C]*tke[2];
     #endif
     epsprof[in] = NA*eps[0] + NB*eps[1] + NC*eps[2];		// eps
     uprof[in] = NA*u_[AN+0] + NB*u_[BN+0] + NC*u_[CN+0];	// <U>
     vprof[in] = NA*u_[AN+1] + NB*u_[BN+1] + NC*u_[CN+1];	// <V>
     u2prof[in] = NA*u2_[UA+0] + NB*u2_[UB+0] + NC*u2_[UC+0];	// <uu>
     v2prof[in] = NA*u2_[UA+1] + NB*u2_[UB+1] + NC*u2_[UC+1];	// <vv>
     u3prof[in] = NA*u3_[AN+0] + NB*u3_[BN+0] + NC*u3_[CN+0];	// skewness u
     v3prof[in] = NA*u3_[AN+1] + NB*u3_[BN+1] + NC*u3_[CN+1];	// skewness v
     u4prof[in] = NA*u4_[AN+0] + NB*u4_[BN+0] + NC*u4_[CN+0];	// flatness u
     v4prof[in] = NA*u4_[AN+1] + NB*u4_[BN+1] + NC*u4_[CN+1];	// flatness v
     uvprof[in] = NA*u2_[UA+3] + NB*u2_[UB+3] + NC*u2_[UC+3];	// <uv>    
 
     (*n)++;							// increase number of profile-points
     in = i*maxn+(*n);	// re-compute array index

     if ( needboth )	// if edge is aligned with upper one, save the upper intersection point as well
     {
        // save the upper one
        yprof[in] = fmax(y1,y2);
	// evaluate shapefunctions at particle location (Cramer's rule)
        NA = (DL[i]*(y[1]-y[2]) + x[1]*(y[2]-yprof[in]) + x[2]*(yprof[in]-y[1]))/dete[e];
        NB = (x[0]*(yprof[in]-y[2]) + DL[i]*(y[2]-y[0]) + x[2]*(y[0]-yprof[in]))/dete[e];
        NC = (x[0]*(y[1]-yprof[in]) + x[1]*(yprof[in]-y[0]) + DL[i]*(y[0]-y[1]))/dete[e];
        // interpolate nodal values
        #ifdef MICROMIXING
        cprof[in] = NA*c_[A] + NB*c_[B] + NC*c_[C];		// scalar mean
        c2prof[in] = NA*c2_[A] + NB*c2_[B] + NC*c2_[C];		// scalar variance
        c3prof[in] = NA*c3_[A] + NB*c3_[B] + NC*c3_[C];		// scalar skewness
        c4prof[in] = NA*c4_[A] + NB*c4_[B] + NC*c4_[C];		// scalar kurtosis
        if (cprof[in] > (*cpeak)) (*cpeak) = cprof[in];		// find peak of mean
        if (c2prof[in] > (*c2peak)) (*c2peak) = c2prof[in];	// find peak of variance
        tmprof[in] = NA*tm_[A] + NB*tm_[B] + NC*tm_[C];		// tm
        #endif
        tke[0] = 0.5*(u2_[UA+0] + u2_[UA+1] + u2_[UA+2]);
        tke[1] = 0.5*(u2_[UB+0] + u2_[UB+1] + u2_[UB+2]);
        tke[2] = 0.5*(u2_[UC+0] + u2_[UC+1] + u2_[UC+2]);
        tkeprof[in] = NA*tke[0] + NB*tke[1] + NC*tke[2];	// tke
        #ifndef WALLFUNCTIONS
        eps[0] = f_[A]*(tke[0] + CT2*f_[A]/RE);
        eps[1] = f_[B]*(tke[1] + CT2*f_[B]/RE);
        eps[2] = f_[C]*(tke[2] + CT2*f_[C]/RE);
        #else
        eps[0] = f_[A]*tke[0];
        eps[1] = f_[B]*tke[1];
        eps[2] = f_[C]*tke[2];
        #endif
        epsprof[in] = NA*eps[0] + NB*eps[1] + NC*eps[2];	// eps
        uprof[in] = NA*u_[AN+0] + NB*u_[BN+0] + NC*u_[CN+0];	// <U>
        vprof[in] = NA*u_[AN+1] + NB*u_[BN+1] + NC*u_[CN+1];	// <V>
        u2prof[in] = NA*u2_[UA+0] + NB*u2_[UB+0] + NC*u2_[UC+0];// <uu>
        v2prof[in] = NA*u2_[UA+1] + NB*u2_[UB+1] + NC*u2_[UC+1];// <vv>
        u3prof[in] = NA*u3_[AN+0] + NB*u3_[BN+0] + NC*u3_[CN+0];// skewness u
        v3prof[in] = NA*u3_[AN+1] + NB*u3_[BN+1] + NC*u3_[CN+1];// skewness v
        u4prof[in] = NA*u4_[AN+0] + NB*u4_[BN+0] + NC*u4_[CN+0];// flatness u
        v4prof[in] = NA*u4_[AN+1] + NB*u4_[BN+1] + NC*u4_[CN+1];// flatness v
        uvprof[in] = NA*u2_[UA+3] + NB*u2_[UB+3] + NC*u2_[UC+3];// <uv>
     
        (*n)++;   
     }
  }
}







static void out_profiles( int ndl, int maxn, double maxy, int *es1, int *es2, int *inpoel, double *coord,
		 	  double *uprof, double *vprof, double *u2prof, double *v2prof, double *u3prof,
                          double *v3prof, double *u4prof, double *v4prof, double *uvprof, double *tkeprof,
                          double *epsprof, double *dete, double *yprof, double *u_, double *u2_,
                          double *u3_, double *u4_, double *f_,
                          #ifdef MICROMIXING
                          double *c_, double *c2_, double *c3_, double *c4_,
                          double *cprof, double *c2prof, double *c3prof, double *c4prof,
                          double *tmprof, double *tm_,
                          #endif
			  const char *cross_filename, const char *downstream_filename )
//
// outputs cross-stream profiles at selected locations
//
{
  int i, n, e, ie, p;
  char s[STRLEN];
  double cpeak, c2peak, dtmp;
  FILE *cfile, *dfile;
  
  

  // open file for downstream evolution profile
  if ( !(dfile = fopen(downstream_filename,"w")) ) ERR("Cannot open file\n");

  #ifdef _OPENMP
  #pragma omp parallel for private(i,n,e,ie,p,cpeak,c2peak,dtmp,s,cfile)
  #endif
  for ( i = 0; i < ndl; i++ )	// loop over downstream locations
  {
     n=0; cpeak=c2peak=0.0001;
     for ( e = es2[i]+1; e <= es2[i+1]; e++ )	// loop over all elements in stripe
        extract_element_profilepoint_values( es1[e], i, &n, maxn, maxy, yprof, dete, inpoel, coord, u_, u2_, u3_,
                                             u4_, f_, uprof, vprof, u2prof, v2prof, u3prof, v3prof, u4prof,
                                             v4prof, uvprof, tkeprof, epsprof
                                             #ifdef MICROMIXING
                                             , &cpeak, &c2peak, c_, c2_, c3_, c4_, tm_, cprof, c2prof, c3prof,
                                             c4prof, tmprof
                                             #endif
                                           );

     // (bubble-)sort cross-stream profile-points
     for ( p = 1; p < n; p++ )
      for ( e = 0; e < (n-p); e++ )
      {
        ie = i*maxn+e;				// compute array index
        if ( yprof[ie] > yprof[ie+1] )
        {
          SWAP(yprof[ie],yprof[ie+1],dtmp);
          SWAP(tkeprof[ie],tkeprof[ie+1],dtmp);
          SWAP(epsprof[ie],epsprof[ie+1],dtmp);
          SWAP(uprof[ie],uprof[ie+1],dtmp);
          SWAP(vprof[ie],vprof[ie+1],dtmp);
          SWAP(u2prof[ie],u2prof[ie+1],dtmp);
          SWAP(v2prof[ie],v2prof[ie+1],dtmp);
          SWAP(u3prof[ie],u3prof[ie+1],dtmp);
          SWAP(v3prof[ie],v3prof[ie+1],dtmp);
          SWAP(u4prof[ie],u4prof[ie+1],dtmp);
          SWAP(v4prof[ie],v4prof[ie+1],dtmp);
          SWAP(uvprof[ie],uvprof[ie+1],dtmp);
          #ifdef MICROMIXING
          SWAP(cprof[ie],cprof[ie+1],dtmp);
          SWAP(c2prof[ie],c2prof[ie+1],dtmp);
          SWAP(c3prof[ie],c3prof[ie+1],dtmp);
          SWAP(c4prof[ie],c4prof[ie+1],dtmp);
          SWAP(tmprof[ie],tmprof[ie+1],dtmp);
          #endif
        }
      }

     // output cross-stream profile
     // construct filename and open file
     sprintf( s, "%s.%.4g", cross_filename, DL[i] );
     if ( !(cfile = fopen(s,"w")) ) ERR("Cannot open file\n");

     // output normalized profiles:
     // y+, <U>, <V>, <uu>, <vv>, skewness u, skewness v, flatness u, flatness v, <uv>,
     // tke, eps, <C>, <c2>, <c3>, <c4>, tm
     for ( e = 0; e < n; e++ )	// loop over profile-points
     {
       ie = i*maxn+e;		// compute array index

       fprintf( cfile, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
                       #ifdef MICROMIXING
                       "\t%g\t%g\t%g\t%g\t%g"
                       #endif
                       "\n",
                       yprof[ie], uprof[ie], vprof[ie], u2prof[ie], v2prof[ie], u3prof[ie],
                       v3prof[ie], u4prof[ie], v4prof[ie], uvprof[ie], tkeprof[ie], epsprof[ie]
                       #ifdef MICROMIXING
                       , cprof[ie]/cpeak, c2prof[ie]/c2peak, c3prof[ie], c4prof[ie], tmprof[ie]
                       #endif
              );
     }

     fclose( cfile );

     // output downstream evolution point: x, peak(<C>), peak(<cc>)
     fprintf( dfile, "%g\t%g\t%g\n", DL[i], cpeak, c2peak );
  }

  fclose( dfile );
}





static void extract_element_profilepoint_values_c( int e, int *n, double maxx, double *xprof_c, double *uprof_c,
                                                   double *pprof_c,
                                                   double *dete, int *inpoel, double *coord, double *u_, double *pr_ )

//
// extract values at the centerline from element e
//
// e : 2d element index,
// n : index of profile-point (gets increased here)
//
{
  int p, eN, A, B, C, AN, BN, CN, needboth=0;
  double dtmp, x1, x2, NA, NB, NC;
  double x[NNODE], y[NNODE];

  eN = e*NNODE;

  A = inpoel[eN+0];		// get node information
  B = inpoel[eN+1];
  C = inpoel[eN+2];
  AN = A*NDIMN;
  BN = B*NDIMN;
  CN = C*NDIMN;

  // order nodes according to non-decreasing y coordinate in index-order A, C, B
  if ( coord[BN+1] < coord[AN+1] ) { p=A; A=B; B=p; }
  if ( coord[CN+1] < coord[AN+1] ) { p=A; A=C; C=p; }
  if ( coord[BN+1] < coord[CN+1] ) { p=B; B=C; C=p; }

  // compute edge-intersections
  if ( (0.0 > coord[AN+1]) && (0.0 < coord[BN+1]) )	// if element in bounding stripe
  {
     if ( 0.0 < coord[CN+1] )	// AC edge
     {
       if ( fabs(dtmp=coord[CN+1]-coord[AN+1]) > EPS )	// if not parallel
         x1 = coord[AN+0] - coord[AN+1]*(coord[CN+0]-coord[AN+0])/dtmp;
       else x1 = 0.5*(coord[AN+0]+coord[CN+0]);		// if parallel: average xA and xC
     }
     else			// CB edge
       if ( fabs(dtmp=coord[CN+1]-coord[BN+1]) > EPS )	// if not parallel
         x1 = coord[BN+0] - coord[BN+1]*(coord[CN+0]-coord[BN+0])/dtmp;
       else x1 = 0.5*(coord[BN+0]+coord[CN+0]);		// if parallel: average xB and xC

     // AB edge
     x2 = coord[AN+0] - coord[AN+1]*(coord[BN+0]-coord[AN+0]) / (coord[BN+1]-coord[AN+1]);

     // if an edge is aligned with the outflow, count both intersection points
     if ( fabs(coord[AN+0]-maxx)<EPS )
     {
       if ( (fabs(coord[BN+0]-maxx)<EPS) || (fabs(coord[CN+0]-maxx)<EPS) ) needboth = 1;
     }
     else
       if ( (fabs(coord[BN+0]-maxx)<EPS) && (fabs(coord[CN+0]-maxx)<EPS) ) needboth = 1;

     A = inpoel[eN+0];		// get node information again
     B = inpoel[eN+1];
     C = inpoel[eN+2];
     AN = A*NDIMN;
     BN = B*NDIMN;
     CN = C*NDIMN;
     x[0] = coord[AN+0];  y[0] = coord[AN+1];
     x[1] = coord[BN+0];  y[1] = coord[BN+1];
     x[2] = coord[CN+0];  y[2] = coord[CN+1];

     // save the leftmost one
     xprof_c[*n] = fmin(x1,x2);
     // evaluate shapefunctions at intersection (Cramer's rule)
     NA = (xprof_c[*n]*(y[1]-y[2]) + x[1]*y[2]-x[2]*y[1])/dete[e];
     NB = (xprof_c[*n]*(y[2]-y[0]) + x[2]*y[0]-x[0]*y[2])/dete[e];
     NC = (xprof_c[*n]*(y[0]-y[1]) + x[0]*y[1]-x[1]*y[0])/dete[e];
     // interpolate nodal values
     uprof_c[*n] = NA*u_[AN+0] + NB*u_[BN+0] + NC*u_[CN+0];	// <U>
     pprof_c[*n] = NA*pr_[A] + NB*pr_[B] + NC*pr_[C];		// 2<P>
 
     // increase number of profile-points
     (*n)++;

     if ( needboth )	// if edge is aligned with upper one, save the upper intersection point as well
     {
        // save the leftmost one
        xprof_c[*n] = fmin(x1,x2);
        // evaluate shapefunctions at intersection (Cramer's rule)
        NA = (xprof_c[*n]*(y[1]-y[2]) + x[1]*y[2]-x[2]*y[1])/dete[e];
        NB = (xprof_c[*n]*(y[2]-y[0]) + x[2]*y[0]-x[0]*y[2])/dete[e];
        NC = (xprof_c[*n]*(y[0]-y[1]) + x[0]*y[1]-x[1]*y[0])/dete[e];
        // interpolate nodal values
        uprof_c[*n] = NA*u_[AN+0] + NB*u_[BN+0] + NC*u_[CN+0];	// <U>
        pprof_c[*n] = NA*pr_[A] + NB*pr_[B] + NC*pr_[C];	// 2<P>

        (*n)++;   
     }
  }
}









static void out_centerline_profiles( int ec_size, int *ec, double maxx, int *inpoel, double *coord,
                                     double *xprof_c, double *uprof_c, double *pprof_c, double *dete,
                                     double *u_, double *pr_,
                                     const char *centerline_filename )
//
// outputs centerline streamwise profiles behind the cylinder
//
{
  int n, p, e;
  double dtmp;
  FILE *cfile;


  for ( n=e=0; e < ec_size; e++ )	// loop over all elements in centerline stripe
    extract_element_profilepoint_values_c( ec[e], &n, maxx, xprof_c, uprof_c, pprof_c, dete, inpoel, coord,
                                           u_, pr_ );


  // (bubble-)sort profile-points
  for ( p = 1; p < n; p++ )
    for ( e = 0; e < (n-p); e++ )
    {
      if ( xprof_c[e] > xprof_c[e+1] )
      {
        SWAP(xprof_c[e],xprof_c[e+1],dtmp);
        SWAP(uprof_c[e],uprof_c[e+1],dtmp);
        SWAP(pprof_c[e],pprof_c[e+1],dtmp);
      }
    }


  // open file for centerline streamwise evolution profiles
  if ( !(cfile = fopen(centerline_filename,"w")) ) ERR("Cannot open file\n");

  // output profiles:
  // x, <U>, 2<P>
  for ( e = 0; e < n; e++ )	// loop over profile-points
  {
    fprintf( cfile, "%g\t%g\t%g\n", xprof_c[e], uprof_c[e], 2.0*pprof_c[e] );
  }

  fclose( cfile );
}







static void out_wall_profiles( int angle_size, int *inprof_w, double *angle_w, double *pr_, double *du_,
                               const char *wall_filename )
//
// outputs profiles along the cylinder wall
//
{
  int p, p4;
  FILE *wfile;


  // open file for centerline streamwise evolution profiles
  if ( !(wfile = fopen(wall_filename,"w")) ) ERR("Cannot open file\n");

  // output profiles:
  // angle, 2<P>, <Wz>/sqrt(Re/2)
  for ( p = 0; p < angle_size; p++ )	// loop over profile-points
  {
    p4 = inprof_w[p]*4;
    fprintf( wfile, "%g\t%g\t%g\n", angle_w[p], 2.0*pr_[inprof_w[p]], -0.5*(du_[p4+2]-du_[p4+1])/sqrtREp2 );
  }
    
  fclose( wfile );
}









static void out_az_profiles( int *az1, int *az2, int *inpoel,
                             double *u_, double *u2_,
			     double *r_az, double *dir_az, double *N_az,
                             const char *az_base_filename )
//
// outputs profiles along azimuthal lines
//
{
  int z, zN, j, m, mN, e, eN, A, B, C, AN, BN, CN, AU, BU, CU;
  double v[NDIMN], v2[U2DOF], v2w[9], k[9], k2[9], t[9];
  char s[STRLEN];
  FILE *azfile[MAXAZ];


  // open files for each azimuthal profiles
  #ifdef _OPENMP
  #pragma omp parallel for private(z,s)
  #endif
  for ( z = 0; z < MAXAZ; z++ )				// loop over all azimuthal angles
  {
    sprintf( s, "%s.%.4g", az_base_filename, AZ[z] );	// construct filename for aizumthal direction
    if ( !(azfile[z] = fopen(s,"w")) ) ERR("Cannot open file\n");
  }


  // output azimuthal profiles
  #ifdef _OPENMP
  #pragma omp parallel for private(z,zN,m,mN,j,e,eN,A,B,C,AN,BN,CN,AU,BU,CU,v,v2,v2w,k,k2,t)
  #endif
  for ( z = 0; z < MAXAZ; z++ )				// loop over all azimuthal angles
  {
    zN = z*NDIMN;
      
    // construct transformation matrix for tensor
    t[0] = -dir_az[zN+1]; t[1] = dir_az[zN+0];  t[2] = 0.0;
    t[3] = dir_az[zN+0];  t[4] = dir_az[zN+1];  t[5] = 0.0;
    t[6] = 0.0;           t[7] = 0.0;           t[8] = 1.0;

    for ( m=0, j=az2[z]+1; j <= az2[z+1]; j++ )		// loop over all element intersected by angle z
    {
      e = az1[j];					// get element index
      eN = e*NNODE;					// get element information
      A = inpoel[eN+0];
      B = inpoel[eN+1];
      C = inpoel[eN+2];
      AN = A*NDIMN;
      BN = B*NDIMN;
      CN = C*NDIMN;
      AU = A*U2DOF;
      BU = B*U2DOF;
      CU = C*U2DOF;

      // interpolate values at intersection
      mN = (az2[z]+1+m)*NNODE;	// compute relative index for shapefunctions
      v[0] = N_az[mN+0]*u_[AN+0] + N_az[mN+1]*u_[BN+0] + N_az[mN+2]*u_[CN+0];
      v[1] = N_az[mN+0]*u_[AN+1] + N_az[mN+1]*u_[BN+1] + N_az[mN+2]*u_[CN+1];
      v2[0] = N_az[mN+0]*u2_[AU+0] + N_az[mN+1]*u2_[BU+0] + N_az[mN+2]*u2_[CU+0];
      v2[1] = N_az[mN+0]*u2_[AU+1] + N_az[mN+1]*u2_[BU+1] + N_az[mN+2]*u2_[CU+1];
      v2[2] = N_az[mN+0]*u2_[AU+2] + N_az[mN+1]*u2_[BU+2] + N_az[mN+2]*u2_[CU+2];
      v2[3] = N_az[mN+0]*u2_[AU+3] + N_az[mN+1]*u2_[BU+3] + N_az[mN+2]*u2_[CU+3];
      v2[4] = N_az[mN+0]*u2_[AU+4] + N_az[mN+1]*u2_[BU+4] + N_az[mN+2]*u2_[CU+4];
      v2[5] = N_az[mN+0]*u2_[AU+5] + N_az[mN+1]*u2_[BU+5] + N_az[mN+2]*u2_[CU+5];

      // transform Reynolds stress tensor into wall-coordinate system
      k[0] = v2[0];   k[1] = v2[3];   k[2] = v2[4];	// copy out symmetric tensor
      k[3] = v2[3];   k[4] = v2[1];   k[5] = v2[5];
      k[6] = v2[4];   k[7] = v2[5];   k[8] = v2[2];
      m3mult( t, k, k2 );				// transform Reynolds stress tensor
      m3mult( t, k2, v2w );

      // output profiles:
      // radial distance,
      // tangential mean velocity,
      // Reynolds stress in cylindrical coordinate system: <uu>, <vv>, <ww>, <uv>
      fprintf( azfile[z], "%g\t%g\t%g\t%g\t%g\t%g\n", r_az[az2[z]+1+m],
                                          v[0]*dir_az[zN+1] + v[1]*dir_az[zN+0],
                                          v2w[0], v2w[4], v2w[8], v2w[1] );

      m++;		// increase number of elements intersected (for angle z)
    }
  }


  // close all azimuthal profile files
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( z = 0; z < MAXAZ; z++ )				// loop over all azimuthal angles
    fclose(azfile[z]);
}







void out_surface_stats( double *tdrag, double *tlift )
//
// outputs time-averaged statistics on the cylinder surface
//
{
  FILE *opos;


  if ( !(opos = fopen(OUT_POSTPROCESS_SURF_FILENAME,"w")) ) ERR("Cannot open OUT_POSTPROCESS_SURF_FILENAME");

  fprintf( opos, "mean pressure drag =\t%g\n", 2.0*tdrag[0] );
  fprintf( opos, "mean viscous drag =\t%g\n", 2.0*tdrag[1] );
  fprintf( opos, "mean pressure lift =\t%g\n", 2.0*tlift[0] );
  fprintf( opos, "mean viscous lift =\t%g\n", 2.0*tlift[1] );

  fclose(opos);
}







int intersect( double *a, double *b, double *f, double *g )
//
// compute the intersection of a line (starting from (0,0) going through "f") and
// an edge given by points "a" and "b"
//
// returns 1 if the intersection is in between A and B, 0 otherwise
// returns in "g" the coordinates of the intersection
//
{
  double d, l;


  // compute denominator
  d = (b[0]-a[0])*f[1] - (b[1]-a[1])*f[0];

  if ( fabs(d) < EPS )			// if parallel
    return(0);				// exit with no intersection
  else
    l = (a[1]*f[0] - a[0]*f[1]) / d;	// compute parameter in edge equation

  //printf( "%g: (%g, %g)\t(%g, %g)\t%g\n", l, a[0],a[1],b[0],b[1],d );

  if ( (l > 0.0) && (l < 1.0) )		// in between points "a" and "b"
  {
    g[0] = a[0] + l*(b[0]-a[0]);	// compute intersection coordinates
    g[1] = a[1] + l*(b[1]-a[1]);
    return(1);
  }

  return(0);
}








static void out_spatially_averaged_fields( int npoin, int nbpoin, double t, double dt, int minnp, int pit,
                                           double *drag, double *lift, double *u, int *betags, int *binpoel,
                                           double *wenr, double *wel, double *pr, double *du, double *u2, double *f
					   #ifndef WALLFUNCTIONS
					   , int rit
					   #endif
					   #ifdef MICROMIXING
					   , double *c, double *c2
					   #endif
					 )
//
// outputs spatially-averaged fields at time t
//
{
  int i, iN, nN, p, p4, pN, pU, nwe;
  double au, av, eps, uu, vv, ww, uv, uw, vw, ap;
  double nr[NDIMN];
  FILE *opos;
  #ifdef MICROMIXING
  double ac, ac2;
  #endif


  // compute spatially averaged fields
  au = av = uu = vv = ww = uv = uw = vw = eps = ap = 0.0;
  #ifdef MICROMIXING
  ac = ac2 = 0.0;
  #endif
  
  #ifdef _OPENMP
    #ifdef MICROMIXING
      #pragma omp parallel for private(p,pU,pN) reduction(+:au,av,uu,vv,ww,uv,uw,vw,eps,ap,ac,ac2)
    #else
      #pragma omp parallel for private(p,pU,pN) reduction(+:au,av,uu,vv,ww,uv,uw,vw,eps,ap)
    #endif // MICROMIXING
  #endif // _OPENMP
  for ( p = 0; p < npoin; p++ )
  {
    pU = p*U2DOF;
    pN = p*NDIMN;
    au += u[pN+0];
    av += u[pN+1];
    uu += u2[pU+0];
    vv += u2[pU+1];
    ww += u2[pU+2];
    uv += u2[pU+3];
    uw += u2[pU+4];
    vw += u2[pU+5];
    #ifndef WALLFUNCTIONS
    eps += f[p]*(0.5*(u2[pU+0] + u2[pU+1] + u2[pU+2]) + CT*CT*f[p]/RE);
    #else
    eps += f[p]*0.5*(u2[pU+0] + u2[pU+1] + u2[pU+2]);
    #endif // WALLFUNCTIONS
    ap += pr[p];
    #ifdef MICROMIXING
    ac += c[p];
    ac2 += c2[p];
    #endif
  }
  au /= npoin;
  av /= npoin;
  uu /= npoin;
  vv /= npoin;
  ww /= npoin;
  uv /= npoin;
  uw /= npoin;
  vw /= npoin;
  eps /= npoin;
  ap /= npoin;
  #ifdef MICROMIXING
  ac /= npoin;
  ac2 /= npoin;
  #endif


  // compute drag and lift on wall-elements
  drag[0] = drag[1] = lift[0] = lift[1] = 0.0;
  for ( nwe=i=0; i < nbpoin; i++ )		// loop over all boundary elements
    if (betags[i*2+0] == 4) 			// cylinder wall
    {
       iN = i*NBNODE;				// get the index of one of the nodes of boundary element i
       p = binpoel[iN+0];
       pU = U2DOF*binpoel[iN+0];
       p4 = 4*binpoel[iN+0];

       nN = nwe*NDIMN;
       nr[0] = -wenr[nN+0];			// get inward pointing normal of wall-element nwe
       nr[1] = -wenr[nN+1];

       // compute different components of drag and lift separately
       drag[0] += -pr[p]*nr[0]*wel[nwe];				// pressure drag
       drag[1] += (du[p4+0]*nr[0] + du[p4+1]*nr[1])/RE*wel[nwe];	// viscous drag
       lift[0] += -pr[p]*nr[1]*wel[nwe];				// pressure lift
       lift[1] += (du[p4+2]*nr[0] + du[p4+3]*nr[1])/RE*wel[nwe];	// viscous lift

       nwe++;					// increase wall-element counter
    }


  // output space-averaged statistics and cylinder-wall integral statistics
  if ( !(opos = fopen(OUT_POSTPROCESS_TIME_FILENAME,"a")) ) ERR("Cannot open OUT_POSTPROCESS_TIME_FILENAME");
  // 1: time, <U>, <V>, <uu>, <vv>, <ww>, <uv>, <uw>, <vw>, eps, pr, minnpel, dt, <C>, <cc>,
  // #ifndef WALLFUNCTIONS
  // 16: rit,
  // #endif
  // 17: pit, pressure drag, viscous drag, pressure lift, viscous lift
  fprintf( opos, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%g"
                 #ifdef MICROMIXING
                 "\t%g\t%g"
		 #endif
		 #ifndef WALLFUNCTIONS
		 "\t%d"
		 #endif
		 "\t%d\t%g\t%g\t%g\t%g\n",
		 t, au, av, uu, vv, ww, uv, uw, vw, eps, ap, minnp, dt,
		 #ifdef MICROMIXING
		 ac, ac2,
		 #endif
		 #ifndef WALLFUNCTIONS
		 rit,
		 #endif
		 pit, 2*drag[0], 2*drag[1], 2*lift[0], 2*lift[1] );
  fclose( opos );
}








static void comp_time_averaged_fields( int npoin, double *u, double *pr, double *du, double *tu, double *tdu,
                                       double *tpr,
                                       double *drag, double *lift, double *tdrag, double *tlift,
                                       double *u2, double *u3, double *u4, double *f, double *tu2,
                                       double *tu3, double *tu4, double *tf,
				       #ifdef MICROMIXING
				       double *c, double *c2, double *c3, double *c4, double *tm,
				       double *uc, double *tc, double *tc2, double *tc3, double *tc4,
				       double *ttm, double *tuc,
				       #endif
				       int *onum )
//
// computes time-averaged fields
//
{
  int p, pN, pU, p4;
  
	
  (*onum)++;	// increase number of samples of time-averaging

  // calculate time-averaged fields
  #ifdef _OPENMP
    #pragma omp parallel for private(p,pN,pU,p4)
  #endif // _OPENMP
  for ( p = 0; p < npoin; p++ )
  {
    pN = p*NDIMN;
    pU = p*U2DOF;
    p4 = p*4;
    if ((*onum) == 1)
    {
      tu[pN+0] = u[pN+0];
      tu[pN+1] = u[pN+1];
      tu2[pU+0] = u2[pU+0];
      tu2[pU+1] = u2[pU+1];
      tu2[pU+2] = u2[pU+2];
      tu2[pU+3] = u2[pU+3];
      tu2[pU+4] = u2[pU+4];
      tu2[pU+5] = u2[pU+5];
      tu3[pN+0] = u3[pN+0];
      tu3[pN+1] = u3[pN+1];
      tu4[pN+0] = u4[pN+0];
      tu4[pN+1] = u4[pN+1];
      tf[p] = f[p];
      tpr[p] = pr[p];
      tdu[p4+0] = du[p4+0];
      tdu[p4+1] = du[p4+1];
      tdu[p4+2] = du[p4+2];
      tdu[p4+3] = du[p4+3];
      #ifdef MICROMIXING
      tc[p] = c[p];
      tc2[p] = c2[p];
      tc3[p] = c3[p];
      tc4[p] = c4[p];
      ttm[p] = tm[p];
      tuc[p*3+0] = uc[p*3+0];
      tuc[p*3+1] = uc[p*3+1];
      tuc[p*3+2] = uc[p*3+2];
      #endif
    }
    else
    {
      tu[pN+0] = ((double)((*onum)-1)/(*onum))*(tu[pN+0] + u[pN+0]/((*onum)-1));
      tu[pN+1] = ((double)((*onum)-1)/(*onum))*(tu[pN+1] + u[pN+1]/((*onum)-1));
      tu2[pU+0] = ((double)((*onum)-1)/(*onum))*(tu2[pU+0] + u2[pU+0]/((*onum)-1));
      tu2[pU+1] = ((double)((*onum)-1)/(*onum))*(tu2[pU+1] + u2[pU+1]/((*onum)-1));
      tu2[pU+2] = ((double)((*onum)-1)/(*onum))*(tu2[pU+2] + u2[pU+2]/((*onum)-1));
      tu2[pU+3] = ((double)((*onum)-1)/(*onum))*(tu2[pU+3] + u2[pU+3]/((*onum)-1));
      tu2[pU+4] = ((double)((*onum)-1)/(*onum))*(tu2[pU+4] + u2[pU+4]/((*onum)-1));
      tu2[pU+5] = ((double)((*onum)-1)/(*onum))*(tu2[pU+5] + u2[pU+5]/((*onum)-1));
      tu3[pN+0] = ((double)((*onum)-1)/(*onum))*(tu3[pN+0] + u3[pN+0]/((*onum)-1));
      tu3[pN+1] = ((double)((*onum)-1)/(*onum))*(tu3[pN+1] + u3[pN+1]/((*onum)-1));
      tu4[pN+0] = ((double)((*onum)-1)/(*onum))*(tu4[pN+0] + u4[pN+0]/((*onum)-1));
      tu4[pN+1] = ((double)((*onum)-1)/(*onum))*(tu4[pN+1] + u4[pN+1]/((*onum)-1));
      tf[p] = ((double)((*onum)-1)/(*onum))*(tf[p] + f[p]/((*onum)-1));
      tpr[p] = ((double)((*onum)-1)/(*onum))*(tpr[p] + pr[p]/((*onum)-1));
      tdu[p4+0] = ((double)((*onum)-1)/(*onum))*(tdu[p4+0] + du[p4+0]/((*onum)-1));
      tdu[p4+1] = ((double)((*onum)-1)/(*onum))*(tdu[p4+1] + du[p4+1]/((*onum)-1));
      tdu[p4+2] = ((double)((*onum)-1)/(*onum))*(tdu[p4+2] + du[p4+2]/((*onum)-1));
      tdu[p4+3] = ((double)((*onum)-1)/(*onum))*(tdu[p4+3] + du[p4+3]/((*onum)-1));
      #ifdef MICROMIXING
      tc[p] = ((double)((*onum)-1)/(*onum))*(tc[p] + c[p]/((*onum)-1));
      tc2[p] = ((double)((*onum)-1)/(*onum))*(tc2[p] + c2[p]/((*onum)-1));
      tc3[p] = ((double)((*onum)-1)/(*onum))*(tc3[p] + c3[p]/((*onum)-1));
      tc4[p] = ((double)((*onum)-1)/(*onum))*(tc4[p] + c4[p]/((*onum)-1));
      ttm[p] = ((double)((*onum)-1)/(*onum))*(ttm[p] + tm[p]/((*onum)-1));
      tuc[p*3+0] = ((double)((*onum)-1)/(*onum))*(tuc[p*3+0] + uc[p*3+0]/((*onum)-1));
      tuc[p*3+1] = ((double)((*onum)-1)/(*onum))*(tuc[p*3+1] + uc[p*3+1]/((*onum)-1));
      tuc[p*3+2] = ((double)((*onum)-1)/(*onum))*(tuc[p*3+2] + uc[p*3+2]/((*onum)-1));
      #endif
    }
  }

  
  // calculate time-averaged drag and lift
  if ((*onum) == 1)
  {
    tdrag[0] = drag[0];
    tdrag[1] = drag[1];
    tlift[0] = lift[0];
    tlift[1] = lift[1];
  }
  else
  {
    tdrag[0] = ((double)((*onum)-1)/(*onum))*(tdrag[0] + drag[0]/((*onum)-1));
    tdrag[1] = ((double)((*onum)-1)/(*onum))*(tdrag[1] + drag[1]/((*onum)-1));
    tlift[0] = ((double)((*onum)-1)/(*onum))*(tlift[0] + lift[0]/((*onum)-1));
    tlift[1] = ((double)((*onum)-1)/(*onum))*(tlift[1] + lift[1]/((*onum)-1));
  }
}















void out_inst_dist( int npoin, int nbpoin, int minnp, int pit, int maxn, double t, double dt, double *drag,
                    double *lift, double *tdrag, double *tlift, int *betags, int *binpoel, int *onum,
                    double *u, double *wenr, double *wel, double *du, int ndl, double maxx, double maxy,
		    int *es1, int *es2, int ec_size, int *ec, int *az1, int *az2,
		    double *xprof_c, double *uprof_c, double *pprof_c,
                    int *inpoel, double *coord, double *uprof, double *vprof, double *u2prof, double *v2prof,
                    double *u3prof, double *v3prof, double *u4prof, double *v4prof, double *uvprof,
                    double *tkeprof, double *epsprof, double *dete, double *yprof, int angle_size, int *inprof_w,
                    double *angle_w, double *r_az, double *dir_az, double *N_az,
                    #ifndef WALLFUNCTIONS
		    int rit,
		    #endif
		    #ifdef MICROMIXING
		    double *c, double *c2, double *c3, double *c4, double *tm, double *uc,
		    double *tc, double *tc2, double *tc3, double *tc4, double *ttm, double *tuc, int npl, int *epdfloc,
		    int *esupel1, int *esupel2, double *c2e, int *psel1, int *psel2, double *parc, double *ce,
                    int *sl, double *tpdf,
                    double *cprof, double *c2prof, double *c3prof, double *c4prof, double *tmprof,
		    #endif
		    double *pr, double *tu, double *u2, double *u3, double *u4, double *f, double *tu2,
                    double *tu3, double *tu4, double *tf, double *tpr, double *tdu )
//
// outputs 'instantaneous' statistics,
// computes and outputs time-averaged statistics
//
{
  #ifdef MICROMIXING
  int p;
  #endif

  // output spatially averaged fields
  out_spatially_averaged_fields( npoin, nbpoin, t, dt, minnp, pit, drag, lift, u, betags, binpoel,
                                 wenr, wel, pr, du, u2, f
                                 #ifndef WALLFUNCTIONS
				 , rit
				 #endif
                                 #ifdef MICROMIXING
				 , c, c2
				 #endif
			       );

  // output instantaneous profiles at downstream locations
  out_profiles( ndl, maxn, maxy, es1, es2, inpoel, coord, uprof, vprof, u2prof, v2prof, u3prof, v3prof,
                u4prof, v4prof, uvprof, tkeprof, epsprof, dete, yprof, u, u2, u3, u4, f,
                #ifdef MICROMIXING
                c, c2, c3, c4, cprof, c2prof, c3prof, c4prof, tmprof, tm,
                #endif
                STRIPE_BASE_FILENAME, DOWNSTREAM_EVOLUTION_FILENAME );

  // output instantaneous streamwise centerline profiles
  out_centerline_profiles( ec_size, ec, maxx, inpoel, coord, xprof_c, uprof_c, pprof_c, dete, u, pr,
                           STREAMWISE_CENTERLINE_FILENAME );

  // output instantaneous profiles along the cylinder wall
  out_wall_profiles( angle_size, inprof_w, angle_w, pr, du, WALL_FILENAME );

  // output instantaneous profiles along azimuthal lines
  out_az_profiles( az1, az2, inpoel, u, u2, r_az, dir_az, N_az, AZ_BASE_FILENAME );

  if (t > AVTIME)
  {
    // compute time-averaged fields
    comp_time_averaged_fields( npoin, u, pr, du, tu, tdu, tpr, drag, lift, tdrag, tlift, u2, u3, u4, f,
                               tu2, tu3, tu4, tf,
                               #ifdef MICROMIXING
                               c, c2, c3, c4, tm, uc, tc, tc2, tc3, tc4, ttm, tuc,
			       #endif
                               onum );				// <- modified

    #ifdef MICROMIXING
    // output instantaneous and time-averaged scalar pdfs
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( p = 0; p < npl; p++ )
      outpdf( p, 1, *onum, epdfloc, esupel1, esupel2, c2e, psel1, psel2, parc, ce, sl, tpdf );
    #endif
  
    // output time-averaged profiles at downstream locations
    out_profiles( ndl, maxn, maxy, es1, es2, inpoel, coord, uprof, vprof, u2prof, v2prof, u3prof, v3prof,
                  u4prof, v4prof, uvprof, tkeprof, epsprof, dete, yprof, tu, tu2, tu3, tu4, tf,
                  #ifdef MICROMIXING
                  tc, tc2, tc3, tc4, cprof, c2prof, c3prof, c4prof, tmprof, ttm,
                  #endif
	 	  TSTRIPE_BASE_FILENAME, TDOWNSTREAM_EVOLUTION_FILENAME );
    
    // output time-averaged streamwise centerline profiles
    out_centerline_profiles( ec_size, ec, maxx, inpoel, coord, xprof_c, uprof_c, pprof_c, dete, tu, tpr,
                             TSTREAMWISE_CENTERLINE_FILENAME );

    // output time-averaged profiles along the cylinder wall
    out_wall_profiles( angle_size, inprof_w, angle_w, tpr, tdu, TWALL_FILENAME );
    
    // output time-averaged profiles along azimuthal lines
    out_az_profiles( az1, az2, inpoel, tu, tu2, r_az, dir_az, N_az, TAZ_BASE_FILENAME );
  }
  else
  {
    #ifdef MICROMIXING
    // output instantaneous scalar pdfs
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( p = 0; p < npl; p++ )
      outpdf( p, 0, *onum, epdfloc, esupel1, esupel2, c2e, psel1, psel2, parc, ce, sl, tpdf );
    #endif
  }
}
