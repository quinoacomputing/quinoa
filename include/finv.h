 
/* finv.h for ANSI C */

#ifndef FINV_H
#define FINV_H
 
#include "gdef.h"     /* From the library mylib */
#include "fmass.h"
#include "fdist.h"
#include "wdist.h"


double finv_Expon (double u);


double finv_Weibull (double alpha, double u);



double finv_ExtremeValue (double u);


double finv_Logistic (double u);


double finv_Pareto (double c, double u);



double finv_Normal1 (double u);



double finv_Normal2 (double u);



double finv_Normal3 (double u);



double finv_LogNormal (double mu, double sigma, double u);



double finv_JohnsonSB (double alpha, double beta, double a, double b,
                       double u);



double finv_JohnsonSU (double alpha, double beta, double u);



double finv_ChiSquare1 (long k, double u);



double finv_ChiSquare2 (long k, double u);



double finv_Student (long n, double u);



double finv_BetaSymmetric (double p, double u);



double finv_GenericC (wdist_CFUNC F, double par[], double u, int d,
                      int detail);


long finv_GenericD1 (fmass_INFO W, double u);


#if 0
long finv_GenericD2 (wdist_DFUNC F, fmass_INFO W, double u);
#endif



long finv_Geometric (double p, double u);


 
#if 0
long finv_Poisson2 (fmass_INFO W, double u);



long finv_Binomial2 (fmass_INFO W, double u);


#endif
#endif
 

