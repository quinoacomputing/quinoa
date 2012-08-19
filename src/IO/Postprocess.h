// -----------------------------------------------------------------------------
// \file    src/IO/Postprocess.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Postprocess functions
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

void outparpos(int npar, double t, double *parcoord);

void out_3dstep(int nelem, int npoin, int *inpoel, double *coord, double *u,
                double *pr, double *u2, double *u3, double *u4, double *f,
                double *du,
		/*#ifndef WALLFUNCTIONS
		  double *rho,
		#endif*/
		#ifdef MICROMIXING
		double *c, double *c2, double *c3, double *c4, double *tm,
		#endif
		const char *filename, int av);

void out_tav_dist(int nbpoin, double maxx, int *bpg, double *coord,
		  double *tu, double *tu2, double *tf);

void out_surface_stats(double *tdrag, double *tlift);

int intersect(double *a, double *b, double *f, double *g);

void out_inst_dist(int npoin, int nbpoin, int minnp, int pit, int maxn,
                   double t, double dt, double *drag,
                   double *lift, double *tdrag, double *tlift, int *betags,
                   int *binpoel, int *onum,
                   double *u, double *wenr, double *wel, double *du, int ndl,
                   double maxx, double maxy,
		   int *es1, int *es2, int ec_size, int *ec, int *az1, int *az2,
		   double *xprof_c, double *uprof_c, double *pprof_c,
                   int *inpoel, double *coord, double *uprof, double *vprof,
                   double *u2prof, double *v2prof,
                   double *u3prof, double *v3prof, double *u4prof,
                   double *v4prof, double *uvprof,
                   double *tkeprof, double *epsprof, double *dete,
                   double *yprof, int angle_size, int *inprof_w,
                   double *angle_w, double *r_az, double *dir_az, double *N_az,
                   #ifndef WALLFUNCTIONS
		   int rit,
		   #endif
		   #ifdef MICROMIXING
		   double *c, double *c2, double *c3, double *c4, double *tm,
                   double *uc,
		   double *tc, double *tc2, double *tc3, double *tc4,
                   double *ttm, double *tuc, int npl, int *epdfloc,
		   int *esupel1, int *esupel2, double *c2e, int *psel1,
                   int *psel2, double *parc, double *ce,
                   int *sl, double *tpdf,
                   double *cprof, double *c2prof, double *c3prof,
                   double *c4prof, double *tmprof,
		   #endif
		   double *pr, double *tu, double *u2, double *u3, double *u4,
                   double *f, double *tu2,
                   double *tu3, double *tu4, double *tf, double *tpr,
                   double *tdu );
