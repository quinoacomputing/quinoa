 
/* fmass.h for ANSI C */
#ifndef FMASS_H
#define FMASS_H
 

struct fmass_INFO_T;
 
/*
   For better precision in the tails, we keep the cumulative probabilities
   (F) in cdf[s] for s <= smed (i.e. cdf[s] is the sum off all the probabi-
   lities pdf[i] for i <= s),
   and the complementary cumulative probabilities (1 - F) in cdf[s] for
   s > smed (i.e. cdf[s] is the sum off all the probabilities pdf[i]
   for i >= s).
*/ 
struct fmass_INFO_T {
   double *cdf;                    /* cumulative probabilities */
   double *pdf;                    /* probability terms or mass distribution */
   double *paramR;                 /* real parameters of the distribution */
   long *paramI;                   /* integer parameters of the distribution */
   long smin;                      /* pdf[s] = 0 for s < smin */
   long smax;                      /* pdf[s] = 0 for s > smax */  
   long smed;                      /* cdf[s] = F(s) for s <= smed, and 
                                      cdf[s] = bar_F(s) for s > smed */
};
 
typedef struct fmass_INFO_T *fmass_INFO;


extern double fmass_Epsilon;


extern double fmass_MaxLambdaPoisson;  /* = 100000  */

extern double fmass_MaxnBinomial;      /* = 100000  */

extern double fmass_MaxnNegaBin;       /* = 100000  */


double fmass_PoissonTerm1 (double lambda, long s);



fmass_INFO fmass_CreatePoisson (double lambda);



void fmass_DeletePoisson (fmass_INFO W);


  
double fmass_PoissonTerm2 (fmass_INFO W, long s);



double fmass_BinomialTerm3 (long n, double p, long s);



double fmass_BinomialTerm1 (long n, double p, double q, long s);



double fmass_BinomialTerm4 (long n, double p, double p2, long s);



fmass_INFO fmass_CreateBinomial (long n, double p, double q);



void fmass_DeleteBinomial (fmass_INFO W);



double fmass_BinomialTerm2 (fmass_INFO W, long s);


double fmass_NegaBinTerm1 (long n, double p, long s);



fmass_INFO fmass_CreateNegaBin (long n, double p);



void fmass_DeleteNegaBin (fmass_INFO W);


  
double fmass_NegaBinTerm2 (fmass_INFO W, long s);

 
#endif
 

