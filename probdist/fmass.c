/*************************************************************************\
 *
 * Package:        ProbDist
 * File:           fmass.c
 * Environment:    ANSI C
 *
 * Copyright (c) 2002 Pierre L'Ecuyer, DIRO, Université de Montréal.
 * e-mail: lecuyer@iro.umontreal.ca
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted without a fee for private, research,
 * academic, or other non-commercial purposes.
 * Any use of this software in a commercial environment requires a
 * written licence from the copyright owner.
 *
 * Any changes made to this package must be clearly identified as such.
 *
 * In scientific publications which used this software, a reference to it
 * would be appreciated.
 *
 * Redistributions of source code must retain this copyright notice
 * and the following disclaimer.
 *
 * THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
\*************************************************************************/


#include "fmass.h"

#include "util.h"
#include "num.h"
#include "num2.h"

#include <stddef.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>


#define TRACE1(N) printf ("*********   " #N " = %d\n", N);
#define TRACE2(x) printf ("*********   " #x " = %g\n", x);


double fmass_Epsilon = 1.0e-16;

/* When we precompute probability terms until terms are smaller than
   fmass_Epsilon, the last few terms will not be very precise. Instead we
   add terms as small as fmass_Epsilon * EPS_EXTRA to get a few correct digits 
   at the tails of the precomputed distributions. */
static const double EPS_EXTRA = 1 / 100.0;

double fmass_MaxLambdaPoisson = 100000.0;

double fmass_MaxnBinomial = 100000.0;

double fmass_MaxnNegaBin = 100000.0;




/*=========================================================================*/

double fmass_PoissonTerm1 (double lam, long s)
{
   const double lamlim = 20.0;
   double y;
   double x = s;
   double Res;

   if (s < 0)
      return 0.0;

   if ((lam < lamlim) && (x < 2.0 * lamlim)) {
      Res = exp (-lam) * pow (lam, x) / num2_Factorial (s);

   } else {
      y = x * log (lam) - num2_LnGamma (x + 1.0) - lam;
      Res = exp (y);
   }

   return Res;
}

/*=========================================================================*/

fmass_INFO fmass_CreatePoisson (double lam)
{
   double epsilon;
   long i, mid, Nmax;
   long imin, imax;
   double sum;
   fmass_INFO W;
   double *P;                     /* Poisson probability terms */
   double *F;                     /* Poisson cumulative probabilities */

   util_Assert (lam >= 0.0, "fmass_CreatePoisson:   lambda < 0");
   W = (fmass_INFO) util_Malloc (sizeof (struct fmass_INFO_T));
   W->paramI = NULL;
   W->paramR = (double *) util_Malloc (sizeof (double));
   W->paramR[0] = lam;

   /* For lam > fmass_MaxLambdaPoisson, we do not use pre-computed arrays */
   if (lam > fmass_MaxLambdaPoisson) {
      W->pdf = NULL;
      W->cdf = NULL;
      return W;
   }

   /* In theory, the Poisson distribution has an infinite range. But */
   /* for i > Nmax, probabilities should be extremely small. */
   Nmax = (long) (lam + 16 * (2 + sqrt (lam)));
   P = (double *) util_Calloc ((size_t) (1 + Nmax), sizeof (double));
   F = (double *) util_Calloc ((size_t) (1 + Nmax), sizeof (double));

   mid = (long) lam;
   epsilon = EPS_EXTRA * fmass_Epsilon / fmass_PoissonTerm1 (lam, mid);
   /* For large lam, fmass_PoissonTerm1 will lose a few digits of precision */
   /* We shall normalize by explicitly summing all terms >= epsilon */
   sum = P[mid] = 1.0;

   /* Start from the maximum and compute terms > epsilon on each side. */
   i = mid;
   while (i > 0 && P[i] > epsilon) {
      P[i - 1] = P[i] * i / lam;
      i--;
      sum += P[i];
   }
   W->smin = imin = i;

   i = mid;
   while (P[i] > epsilon) {
      P[i + 1] = P[i] * lam / (i + 1);
      i++;
      sum += P[i];
      if (i >= Nmax - 1) {
         Nmax *= 2;
         P = (double *) util_Realloc (P, (1 + Nmax) * sizeof (double));
         F = (double *) util_Realloc (F, (1 + Nmax) * sizeof (double));
         /* util_Warning (TRUE, "fmass_CreatePoisson: Calling Realloc"); */
      }
   }
   W->smax = imax = i;

   /* Renormalize the sum of probabilities to 1 */
   for (i = imin; i <= imax; i++) {
      P[i] /= sum;
   }

   /* Compute the cumulative probabilities until F >= 0.5, and keep them in
      the lower part of array, i.e. F[s] contains all P[i] for i <= s */
   F[imin] = P[imin];
   i = imin;
   while (i < imax && F[i] < 0.5) {
      i++;
      F[i] = P[i] + F[i - 1];
   }
   /* This is the boundary between F and 1 - F in the CDF */
   W->smed = i;
 
   /* Compute the cumulative probabilities of the complementary distribution
      and keep them in the upper part of the array. i.e. F[s] contains all
      P[i] for i >= s */
   F[imax] = P[imax];
   i = imax - 1;
   while (i > W->smed) {
      F[i] = P[i] + F[i + 1];
      i--;
   };

   /* Reset imin because we lose too much precision for a few terms near
      imin when we stop adding terms < epsilon. */
   i = imin;
   while (i < W->smed && F[i] < fmass_Epsilon)
      i++; 
   W->smin = imin = i;

   /* Same thing with imax */
   i = imax;
   while (i > W->smed && F[i] < fmass_Epsilon)
      i--; 
   W->smax = imax = i;

   W->pdf = (double *) util_Calloc ((size_t) (imax + 1 - imin), sizeof (double));
   W->cdf = (double *) util_Calloc ((size_t) (imax + 1 - imin), sizeof (double));
   for (i = imin; i <= imax; i++) {
      W->pdf[i - imin] = P[i];
      W->cdf[i - imin] = F[i];
   }
   util_Free (P);
   util_Free (F);
   return W;
}

/*-------------------------------------------------------------------------*/

double fmass_PoissonTerm2 (fmass_INFO W, long s)
{
   double lam;

   util_Assert (W != NULL,
      "fmass_PoissonTerm2:  fmass_INFO is NULL pointer");
   lam = W->paramR[0];
   if (s < 0)
      return 0.0;
   if (W->pdf == NULL)
      return fmass_PoissonTerm1 (lam, s);
   if (s > W->smax || s < W->smin)
      return fmass_PoissonTerm1 (lam, s);
   return W->pdf[s - W->smin];
}

/*-------------------------------------------------------------------------*/

void fmass_DeletePoisson (fmass_INFO W)
{
   if (W == NULL)
      return;
   util_Free (W->paramR);
   util_Free (W->pdf);
   util_Free (W->cdf);
   util_Free (W);
}


/*=========================================================================*/

double fmass_BinomialTerm1 (long n, double p, double q, long s)
{
   const long slim = 30;          /* To avoid overflow */
   const double maxexp = (DBL_MAX_EXP - 1) * num_Ln2; /* To avoid overflow */
   const double minexp = (DBL_MIN_EXP - 1) * num_Ln2; /* To avoid underflow */
   int signe = 1;
   double Res;

   util_Assert (n >= 0, "fmass_BinomialTerm1:   n < 0");
   if (0 == n)
      return 1.0;
   if (s < 0 || s > n)
      return 0.0;

   /* Combination(n, s) are symmetric between s and n-s */
   if (s > n / 2) {
      s = n - s;
      Res = p;
      p = q;
      q = Res;
   }

   if (p < 0.0) {
      p = -p;
      if (s & 1)
         signe *= -1;             /* odd s */
   }
   if (q < 0.0) {
      q = -q;
      if ((n - s) & 1)
         signe *= -1;             /* odd n - s */
   }

   if (n <= slim) {
      Res = pow (p, (double) s) * num2_Combination (n, s) * pow (q,
         (double) (n - s));
      return signe * Res;

   } else {
      /* This could be calculated with more precision as there is some
         cancellation because of subtraction of the large LnFactorial: the
         last few digits can be lost. But we need the function lgammal in
         long double precision. Another possibility would be to use an
         asymptotic expansion for the binomial coefficient. */
      Res = s * log (p) + (n - s) * log (q) + num2_LnFactorial (n)
         - num2_LnFactorial (n - s) - num2_LnFactorial (s);
      util_Assert (Res < maxexp, "fmass_BinomialTerm1:   term overflow");

      if (Res < minexp)
         return 0.0;

      return signe * exp (Res);
   }
}


/*=========================================================================*/

double fmass_BinomialTerm4 (long n, double p, double p2, long s)
{
   const long slim = 30;          /* To avoid overflow */
   const double maxexp = (DBL_MAX_EXP - 1) * num_Ln2; /* To avoid overflow */
   const double minexp = (DBL_MIN_EXP - 1) * num_Ln2; /* To avoid underflow */
   double Res;

   util_Assert (p >= 0.0 && p <= 1.0, "fmass_BinomialTerm4:   p not in [0, 1]");
   util_Assert (p2 >= 0.0 && p2 <= 1.0, "fmass_BinomialTerm4:   p2 not in [0, 1]");
   util_Assert (n >= 0, "fmass_BinomialTerm4:   n < 0");
   if (0 == n)
      return 1.0;
   if (s < 0 || s > n)
      return 0.0;

   if (n <= slim) {
      if (p2 > 1.0e-1) {
         Res = pow (p, (double) s) * num2_Combination (n, s) * pow (1.0 - p2,
               (double) (n - s));
      } else {
         double temp = (n - s)*num2_log1p (-p2);
         Res = pow (p, (double) s) * num2_Combination (n, s) * exp(temp);
      }
      return Res;

   } else {
      /* This could be calculated with more precision as there is some
         cancellation because of subtraction of the large LnFactorial: the
         last few digits can be lost. But we need the function lgammal in
         long double precision. Another possibility would be to use an
         asymptotic expansion for the binomial coefficient. */
      Res = s * log (p) + (n - s) * num2_log1p(-p2) + num2_LnFactorial (n)
         - num2_LnFactorial (n - s) - num2_LnFactorial (s);
      util_Assert (Res < maxexp, "fmass_BinomialTerm4:   term overflow");

      if (Res < minexp)
         return 0.0;

      return exp (Res);
   }
}


/*=========================================================================*/

double fmass_BinomialTerm3 (long n, double p, long s)
{
   const long slim = 50;          /* To avoid overflow */
   const double maxexp = (DBL_MAX_EXP - 1) * num_Ln2; /* To avoid overflow */
   const double minexp = (DBL_MIN_EXP - 1) * num_Ln2; /* To avoid underflow */
   int signe = 1;
   double Res;
   double q = 1.0 - p;

   /* util_Assert (p >= 0.0 && p <= 1.0, "fmass_BinomialTerm3: p not in [0,
      1]"); */
   util_Assert (n >= 0, "fmass_BinomialTerm3:   n < 0");
   if (0 == n)
      return 1.0;
   if (s < 0 || s > n)
      return 0.0;

   /* Combination(n, s) are symmetric between s and n-s */
   if (s > n / 2) {
      s = n - s;
      Res = p;
      p = q;
      q = Res;
   }

   if (p < 0.0) {
      p = -p;
      if (s & 1)
         signe *= -1;             /* odd s */
   }
   if (q < 0.0) {
      q = -q;
      if ((n - s) & 1)
         signe *= -1;             /* odd n - s */
   }

   if (n <= slim) {
      if (p > 1.0e-1) {
         Res = pow (p, (double) s) * num2_Combination (n, s) * pow (q,
               (double) (n - s));
      } else {
         double temp = (n - s)*num2_log1p (-p);
         Res = pow (p, (double) s) * num2_Combination (n, s) * exp(temp);
      }
      return signe * Res;

   } else {
      /* This could be calculated with more precision as there is some
         cancellation because of subtraction of the large LnFactorial: the
         last few digits can be lost. But we need the function lgammal in
         long double precision. Another possibility would be to use an
         asymptotic expansion for the binomial coefficient. */
      Res = s * log (p) + (n - s) * num2_log1p (-p) + num2_LnFactorial (n)
         - num2_LnFactorial (n - s) - num2_LnFactorial (s);
      util_Assert (Res < maxexp, "fmass_BinomialTerm3:   term overflow");

      if (Res < minexp)
         return 0.0;

      return signe * exp (Res);
   }
}


/*=========================================================================*/

fmass_INFO fmass_CreateBinomial (long n, double p, double q)
{
/* 
 * Compute all probability terms of the binomial distribution; start near
 * the mean, and calculate probabilities on each side until they become
 * smaller than epsilon, then stop there.
 * However, this is more general than the binomial probability distribu-
 * tion as this will compute the binomial terms when p + q != 1, and
 * even when p or q are negative. However in this case, the cumulative
 * terms are meaningless and are not computed.
 */
   const double epsilon = fmass_Epsilon * EPS_EXTRA;
   long i, mid;
   long imin, imax;
   double z = 0;
   fmass_INFO W;
   double *P;                     /* Binomial "probability" terms */
   double *F;                     /* Binomial cumulative "probabilities" */

   util_Assert (n > 0, "fmass_CreateBinomial:  n <= 0");

   W = (fmass_INFO) util_Malloc (sizeof (struct fmass_INFO_T));
   W->paramI = (long *) util_Malloc (sizeof (long));
   W->paramR = (double *) util_Calloc ((size_t) 2, sizeof (double));
   W->paramI[0] = n;
   W->paramR[0] = p;
   W->paramR[1] = q;

   /* For n > fmass_MaxnBinomial, we shall not use pre-computed arrays */
   if (n > fmass_MaxnBinomial) {
      W->pdf = NULL;
      W->cdf = NULL;
      return W;
   }

   P = (double *) util_Calloc ((size_t) (1 + n), sizeof (double));
   F = (double *) util_Calloc ((size_t) (1 + n), sizeof (double));

   /* the maximum term in absolute value */
   mid = (long) ((n + 1) * fabs (p) / (fabs (p) + fabs (q)));
   if (mid > n)
      mid = n;
   P[mid] = fmass_BinomialTerm1 (n, p, q, mid);

   if (fabs(p) > 0.0) {
      z = q / p;
   } else {
      z = 0.0;
      util_Warning (1, "fmass_CreateBinomial:   q / p = infinite");
   }
   i = mid;
   while (i > 0 && fabs (P[i]) > epsilon) {
      P[i - 1] = P[i] * z * i / (n - i + 1);
      i--;
   }
   imin = i;

   if (fabs(q) > 0.0) {
      z = p / q;
   } else {
      z = 0.0;
      util_Warning (1, "fmass_CreateBinomial:   p / q = infinite");
   }
   i = mid;
   while (i < n && fabs (P[i]) > epsilon) {
      P[i + 1] = P[i] * z * (n - i) / (i + 1);
      i++;
   }
   imax = i;

   /* Here, we assume that we are dealing with a probability distribution. */
   /* Compute the cumulative probabilities for F and keep them in the */
   /* lower part of CDF. */
   F[imin] = P[imin];
   i = imin;
   while (i < n && F[i] < 0.5) {
      i++;
      F[i] = F[i - 1] + P[i];
   }

   /* This is the boundary between F (i <= smed) and 1 - F (i > smed) in */
   /* the array CDF */
   W->smed = i;

   /* Compute the cumulative probabilities of the complementary */
   /* distribution and keep them in the upper part of the array */
   F[imax] = P[imax];
   i = imax - 1;
   while (i > W->smed) {
      F[i] = P[i] + F[i + 1];
      i--;
   }

   /* Reset imin because we lose too much precision for a few terms near
      imin when we stop adding terms < epsilon. */
   i = imin;
   while (i < W->smed && F[i] < fmass_Epsilon)
      i++; 
   W->smin = imin = i;

   /* Same thing with imax */
   i = imax;
   while (i > W->smed && F[i] < fmass_Epsilon)
      i--; 
   W->smax = imax = i;

   W->pdf = (double *) util_Calloc ((size_t) (imax + 1 - imin), sizeof (double));
   W->cdf = (double *) util_Calloc ((size_t) (imax + 1 - imin), sizeof (double));
   for (i = imin; i <= imax; i++) {
      W->pdf[i - imin] = P[i];
      W->cdf[i - imin] = F[i];
   }
   util_Free (P);
   util_Free (F);

   return W;
}

/*-------------------------------------------------------------------------*/

double fmass_BinomialTerm2 (fmass_INFO W, long s)
{
   long n;
   double p, q;

   util_Assert (W != NULL,
      "fmass_BinomialTerm2: fmass_INFO is NULL pointer");
   n = W->paramI[0];
   if (0 == n)
      return 1.0;
   if (s < 0 || s > n)
      return 0.0;
   p = W->paramR[0];
   if (p == 0.0) {
      if (s > 0)
         return 0.0;
      else
         return 1.0;
   }
   q = W->paramR[1];
   if (q == 0.0) {
      if (s < n)
         return 0.0;
      else
         return 1.0;
   }
   if (W->pdf == NULL)
      return fmass_BinomialTerm1 (n, p, q, s);

   if (s > W->smax || s < W->smin)
      return fmass_BinomialTerm1 (n, p, q, s);

   return W->pdf[s - W->smin];
}

/*-------------------------------------------------------------------------*/

void fmass_DeleteBinomial (fmass_INFO W)
{
   if (W == NULL)
      return;
   util_Free (W->paramI);
   util_Free (W->paramR);
   util_Free (W->pdf);
   util_Free (W->cdf);
   util_Free (W);
}


/*=========================================================================*/

double fmass_NegaBinTerm1 (long n, double p, long s)
{
   const long slim = 15;          /* To avoid overflow */
   const double maxexp = (DBL_MAX_EXP - 1) * num_Ln2; /* To avoid overflow */
   const double minexp = (DBL_MIN_EXP - 1) * num_Ln2; /* To avoid underflow */
   double y;

   util_Assert (p >= 0.0 && p <= 1.0,
      "fmass_NegaBinTerm1:   p not in [0, 1]");
   util_Assert (n > 0, "fmass_NegaBinTerm1:   n < 1");
   if (s < 0)
      return 0.0;
   if (p >= 1.0) {                /* In fact, p == 1 */
      if (0 == s)
         return 1.0;
      else
         return 0.0;
   }
   if (p <= 0.0)                  /* In fact, p == 0 */
      return 0.0;

   if (s <= slim || n <= slim) {
      y = pow (p, (double) n) * num2_Combination (n + s - 1, s) *
         pow (1.0 - p, (double) s);
      return y;

   } else {
      y = s * num2_log1p (-p) + n * log (p) + num2_LnFactorial (n + s - 1)
         - num2_LnFactorial (n - 1) - num2_LnFactorial (s);
      util_Assert (y < maxexp, "fmass_NegaBinTerm1:   term overflow");
      if (y <= minexp)
         return 0.0;
      else
         return exp (y);
   }
}


/*=========================================================================*/

fmass_INFO fmass_CreateNegaBin (long n, double p)
/* 
 * Compute all probability terms of the negative binomial distribution;
 * start at the mode, and calculate probabilities on each side until they
 * become smaller than epsilon. Set all others to 0.
 */
{
   double epsilon;
   long i, mode, Nmax;
   long imin, imax;
   double sum;
   fmass_INFO W;
   double *P;                     /* Negative Binomial mass probabilities */
   double *F;                     /* Negative Binomial cumulative
                                     probabilities */

   util_Assert (p >= 0.0 && p <= 1.0,
      "fmass_CreateNegaBin:   p not in [0, 1]");
   util_Assert (n > 0, "fmass_CreateNegaBin:  n < 1");

   W = (fmass_INFO) util_Malloc (sizeof (struct fmass_INFO_T));
   W->paramI = (long *) util_Malloc (sizeof (long));
   W->paramR = (double *) util_Malloc (sizeof (double));
   W->paramI[0] = n;
   W->paramR[0] = p;

   /* Compute the mode (at the maximum term) */
   mode = (long) (1 + (n * (1.0 - p) - 1.0) / p);

   /* For mode > fmass_MaxnNegaBin, we shall not use pre-computed arrays.
      mode < 0 should be impossible, unless overflow of long occur, in
      which case mode will be = LONG_MIN. */
   if (mode < 0 || mode > fmass_MaxnNegaBin) {
      W->pdf = NULL;
      W->cdf = NULL;
      return W;
   }

   /* In theory, the negative binomial distribution has an infinite range. */
   /* But for i > Nmax, probabilities should be extremely small. */
   /* Nmax = Mean + 16 * Standard deviation. */
   Nmax = (long) (n * (1.0 - p) / p + 16 * sqrt (n * (1.0 - p) / (p * p)));
   if (Nmax < 32)
      Nmax = 32;
   P = (double *) util_Calloc ((size_t) (1 + Nmax), sizeof (double));
   F = (double *) util_Calloc ((size_t) (1 + Nmax), sizeof (double));

   epsilon = fmass_Epsilon * EPS_EXTRA / fmass_NegaBinTerm1 (n, p, mode);

   /* We shall normalize by explicitly summing all terms >= epsilon */
   sum = P[mode] = 1.0;

   /* Start from the maximum and compute terms > epsilon on each side. */
   i = mode;
   while (i > 0 && P[i] >= epsilon) {
      P[i - 1] = P[i] * i / ((1.0 - p) * (n + i - 1));
      i--;
      sum += P[i];
   }
   imin = i;

   i = mode;
   while (P[i] >= epsilon) {
      P[i + 1] = P[i] * (1.0 - p) * (n + i) / (i + 1);
      i++;
      sum += P[i];
      if (i == Nmax - 1) {
         Nmax *= 2;
         P = (double *) util_Realloc (P, (1 + Nmax) * sizeof (double));
         F = (double *) util_Realloc (F, (1 + Nmax) * sizeof (double));
         /* util_Warning (TRUE, "fmass_CreateNegaBin: Calling Realloc"); */
      }
   }
   imax = i;

   /* Renormalize the sum of probabilities to 1 */
   for (i = imin; i <= imax; i++) {
      P[i] /= sum;
   }

   /* Compute the cumulative probabilities for F and keep them in the */
   /* lower part of CDF. */
   F[imin] = P[imin];
   i = imin;
   while (i < imax && F[i] < 0.5) {
      i++;
      F[i] = F[i - 1] + P[i];
   }

   /* This is the boundary between F (i <= smed) and 1 - F (i > smed) in */
   /* the array CDF */
   W->smed = i;

   /* Compute the cumulative probabilities of the complementary */
   /* distribution 1 - F and keep them in the upper part of the array */
   F[imax] = P[imax];
   i = imax - 1;
   while (i > W->smed) {
      F[i] = P[i] + F[i + 1];
      i--;
   }

   /* Reset imin because we lose too much precision for a few terms near
      imin when we stop adding terms < epsilon. */
   i = imin;
   while (i < W->smed && F[i] < fmass_Epsilon)
      i++; 
   W->smin = imin = i;

   /* Same thing with imax */
   i = imax;
   while (i > W->smed && F[i] < fmass_Epsilon)
      i--; 
   W->smax = imax = i;

   W->pdf = (double *) util_Calloc ((size_t) (imax + 1 - imin), sizeof (double));
   W->cdf = (double *) util_Calloc ((size_t) (imax + 1 - imin), sizeof (double));
   for (i = imin; i <= imax; i++) {
      W->pdf[i - imin] = P[i];
      W->cdf[i - imin] = F[i];
   }
   util_Free (P);
   util_Free (F);

   return W;
}

/*-------------------------------------------------------------------------*/

double fmass_NegaBinTerm2 (fmass_INFO W, long s)
{
   double p;
   long n;

   util_Assert (W != NULL,
      "fmass_NegaBinTerm2:  fmass_INFO is NULL pointer");
   if (s < 0)
      return 0.0;
   n = W->paramI[0];
   p = W->paramR[0];
   if (p == 0.0)
      return 0.0;
   if (p == 1.0) {
      if (s > 0)
         return 0.0;
      else
         return 1.0;
   }

   if (W->pdf == NULL)
      return fmass_NegaBinTerm1 (n, p, s);

   if (s > W->smax || s < W->smin)
      return fmass_NegaBinTerm1 (n, p, s);

   return W->pdf[s - W->smin];
}

/*-------------------------------------------------------------------------*/

void fmass_DeleteNegaBin (fmass_INFO W)
{
   if (W == NULL)
      return;
   util_Free (W->paramI);
   util_Free (W->paramR);
   util_Free (W->pdf);
   util_Free (W->cdf);
   util_Free (W);
}
