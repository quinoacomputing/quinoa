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
//  externally callable functions defined in pressure.cc
//


void pinit( int nelem, int nbpoin, int nthreads, int *inpoel, int *bpg, int *bpsup1, int *bpsup2,
            double maxx, double *coord,
            double *Ae, double *dNx, double *dNy, double *ofnz, sparsemat *P );

int pstep( int nwe, int npoin, int nelem, int nbpoin, int nthreads, double dt, double maxx,
           int *inpoel, int *binpoel,
           int *we, int *weo, int *bpg, int *betags,
	   double *odpn, double *oudt, double *wel, double *wenr, double *prhs, double *du,
	   double *ddu, double *pr, double *dpr, double *u2,
	   double *dNx, double *dNy, double *Ae, double *coord, sparsemat *P, double *u );

void velcorr( int npar, double dt, int *inpoel, int *elp, double *parvel, double *dpr, double *dNx, double *dNy );

void save_old_dpn( int nwe, int nbpoin, int *binpoel, int *betags, double *wenr,
                   double *ddu, double *odpn, double *oudt, int *we, int *inpoel, double *u2,
                   double *dNx, double *dNy, double *u, double *du );
