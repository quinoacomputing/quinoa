// -----------------------------------------------------------------------------
// \file    src/SDE/Langevin.h
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Langevin equation
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

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
             );
