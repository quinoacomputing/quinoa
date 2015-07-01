 
/* fdist.h for ANSI C */
#ifndef FDIST_H
#define FDIST_H
 
#include "gdef.h"
#include "fmass.h"


double fdist_Unif (double x);



double fdist_Expon (double x);



double fdist_Weibull (double alpha, double x);



double fdist_ExtremeValue (double x);



double fdist_Logistic (double x);



double fdist_Pareto (double c, double x);



double fdist_Normal1 (double x);



double fdist_Normal2 (double x);



#ifdef HAVE_ERF
   double fdist_Normal3 (double x);
#endif



double fdist_Normal4 (double x);



double fdist_BiNormal1 (double x, double y, double rho, int ndig);



double fdist_BiNormal2 (double x, double y, double rho);



double fdist_LogNormal (double mu, double sigma, double x);



double fdist_JohnsonSB (double alpha, double beta, double a, double b,
                        double x);



double fdist_JohnsonSU (double alpha, double beta, double x);



double fdist_ChiSquare1 (long k, double x);



double fdist_ChiSquare2 (long k, int d, double x);



double fdist_Student1 (long n, double x);



double fdist_Student2 (long n, int d, double x);



double fdist_Gamma (double a, int d, double x);



double fdist_Beta (double p, double q, int d, double x);



double fdist_BetaSymmetric (double p, double x);



double fdist_KSPlus (long n, double x);



double fdist_KS1 (long n, double x);



double fdist_KS2 (long n, double x);



double fdist_KSPlusJumpOne (long n, double a, double x);

#if 0

void fdist_FindJumps (fdist_FUNC_JUMPS *H, int Detail);



void fdist_FreeJumps (fdist_FUNC_JUMPS *H);



double fdist_KSMinusJumpsMany (fdist_FUNC_JUMPS *H, double x);



double fdist_KSPlusJumpsMany (fdist_FUNC_JUMPS *H, double x);
#endif



double fdist_CramerMises (long n, double x);



double fdist_WatsonG (long n, double x);



double fdist_WatsonU (long n, double x);



double fdist_AndersonDarling (long n, double x);



double fdist_AndersonDarling2 (long n, double x);


double fdist_Geometric (double p, long s);



double fdist_Poisson1 (double lambda, long s);



double fdist_Poisson2 (fmass_INFO W, long s);



double fdist_Binomial1 (long n, double p, long s);



double fdist_Binomial2 (fmass_INFO W, long s);



double fdist_NegaBin1 (long n, double p, long s);



double fdist_NegaBin2 (fmass_INFO W, long s);



double fdist_Scan (long N, double d, long m);

 
#endif
 

