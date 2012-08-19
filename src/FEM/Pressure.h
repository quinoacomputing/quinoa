// -----------------------------------------------------------------------------
// \file    src/FEM/Pressure.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Poisson solution for mean pressure
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

void pinit(int nelem, int nbpoin, int nthreads, int *inpoel, int *bpg,
           int *bpsup1, int *bpsup2,
           double maxx, double *coord,
           double *Ae, double *dNx, double *dNy, double *ofnz, sparsemat *P );

int pstep(int nwe, int npoin, int nelem, int nbpoin, int nthreads, double dt,
          double maxx, int *inpoel, int *binpoel,
          int *we, int *weo, int *bpg, int *betags,
	  double *odpn, double *oudt, double *wel, double *wenr, double *prhs,
          double *du,
	  double *ddu, double *pr, double *dpr, double *u2,
	  double *dNx, double *dNy, double *Ae, double *coord, sparsemat *P,
          double *u );

void velcorr(int npar, double dt, int *inpoel, int *elp, double *parvel,
             double *dpr, double *dNx, double *dNy );

void save_old_dpn(int nwe, int nbpoin, int *binpoel, int *betags, double *wenr,
                  double *ddu, double *odpn, double *oudt, int *we, int *inpoel,
                  double *u2,
                  double *dNx, double *dNy, double *u, double *du );
