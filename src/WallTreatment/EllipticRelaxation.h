// -----------------------------------------------------------------------------
// \file    src/WallTreatment/EllipticRelaxation.h
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Elliptic relaxation in 2D
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

int rstep(int nelem, int nbpoin, int nthreads, int *bpg, int *inpoel, int *wlc,
          int *bpsup1, int *bpsup2, int *bptags, double *Ae, double *wpnr,
          double *rrhs, double *rho, double *dNx, double *dNy, double *u2,
          double *f, double *wnz, double *du, sparsemat *R );
