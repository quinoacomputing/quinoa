 
/* fbar.h for ANSI C */
#ifndef FBAR_H
#define FBAR_H
 
#include "gdef.h"
#include "fmass.h"


double fbar_Unif (double x);



double fbar_Expon (double x);



double fbar_Weibull (double alpha, double x);



double fbar_Logistic (double x);



double fbar_Pareto (double c, double x);



double fbar_Normal1 (double x);



double fbar_Normal2 (double x);



#ifdef HAVE_ERF
   double fbar_Normal3 (double x);
#endif



double fbar_BiNormal1 (double x, double y, double rho, int ndig);



double fbar_BiNormal2 (double x, double y, double rho);



double fbar_ChiSquare1 (long N, double x);



double fbar_ChiSquare2 (long N, int d, double x);



double fbar_Gamma (double a, int d, double x);



double fbar_KS1 (long n, double x);



double fbar_KSPlus (long n, double x);



double fbar_LogNormal (double mu, double sigma, double x);

double fbar_JohnsonSB (double alpha, double beta, double a, double b,
                       double x);

double fbar_JohnsonSU (double alpha, double beta, double x);

double fbar_CramerMises (long n, double x);

double fbar_WatsonU (long n, double x);

double fbar_WatsonG (long n, double x);

double fbar_AndersonDarling (long n, double x);


double fbar_Geometric (double p, long s);



double fbar_Poisson1 (double lambda, long s);



double fbar_Poisson2 (fmass_INFO W, long s);



double fbar_Binomial2 (fmass_INFO W, long s);



double fbar_NegaBin2 (fmass_INFO W, long s);



double fbar_Scan (long N, double d, long m);

 
#endif
 

