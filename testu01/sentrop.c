/*************************************************************************\
 *
 * Package:        TestU01
 * File:           sentrop.c
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


#include "gdef.h"
#include "util.h"
#include "tables.h"
#include "chrono.h"
#include "num.h"

#include "sentrop.h"
#include "smultin.h"
#include "unif01.h"
#include "swrite.h"
#include "wdist.h"

#include "statcoll.h"
#include "gofw.h"
#include "fbar.h"

#include <math.h>
#include <stddef.h>
#include <string.h>


/*------------------------------ Constants --------------------------------*/

/* We use a numerically stable algorithm to compute the variance and the
   correlation in sentrop_EntropyDiscOver and sentrop_EntropyDiscOver2.
   If it is undefined, we use the naive unstable method. */
#define STABLE

 /* Euler constant */
static const double Euler = 0.57721566490153286;

/* Necessary so that the multiplication of a large number of U01 will
  not underflow in sentrop_EntropyDM and sentrop_EntropyDMCirc;
  the precise value is not so important: it could be bigger. */
static const double Epsilon = 1.0e-50;


/* We may choose n >>> NLIM in the tests where xLgx is used, because
 * the probabilities will be concentrated near i = 0, and large values of
 * -i/n * Lg (i/n) will never occur if n/2^L is not too large */
#define NLIM 16384




/*------------------------------ Functions --------------------------------*/


static double FoncMNEntropie (double junk, double n, long j)
/*
 * One term of entropy when there are j balls in an urn
 */
{
   double x;
   if (j == 0)
      return 0.0;
   x = j / n;
   return -x * num_Log2 (x);
}


/*=========================================================================*/

static void CalcLgx (double xLgx[], long n)
{
   long i, nLim;
   double nLR = n;
   double temp;

   if (n < NLIM)
      nLim = n;
   else
      nLim = NLIM;

   for (i = 1; i <= nLim; i++) {
      temp = i / nLR;
      xLgx[i] = -temp * num_Log2 (temp);
   }
   xLgx[0] = 0.0;
}


/*=========================================================================*/

static void InitRes (
   sentrop_Res *res,          /* Results holder */
   long N,                    /* Number of replications */
   int jmax,                  /* Maximal index for Count */
   char *nam
)
/* 
 * Initializes res
 */
{
   sres_InitBasic (res->Bas, N, nam);
   if (jmax > res->jmax)
      res->Count = util_Realloc (res->Count, (jmax + 1) * sizeof (long));
   res->jmax = jmax;
}


/*-------------------------------------------------------------------------*/

sentrop_Res * sentrop_CreateRes (void)
{
   sentrop_Res *res;
   res = util_Malloc (sizeof (sentrop_Res));
   memset (res, 0, sizeof (sentrop_Res));
   res->Bas = sres_CreateBasic ();
   res->Count = util_Calloc (1, sizeof (long));
   res->jmax = 0;
   res->jmin = 0;
   return res;
}


/*-------------------------------------------------------------------------*/

void sentrop_DeleteRes (sentrop_Res *res)
{
   if (res == NULL)
      return;
   sres_DeleteBasic (res->Bas);
   util_Free (res->Count);
   util_Free (res);
}


/*=========================================================================*/

static void InitExactOver (
   long n,                    /* Sample size */
   long L,                    /* Blocks of L bits */
   double *Mu,                /* Mean */
   double *Sigma              /* Standard deviation */
   )
/* 
 * Initializes the exact values of the mean and the standard deviation
 * for the entropy, if these are known (precomputed). Otherwise, sets
 * the standard deviation to -1.
 */
{
   long j;
   long i;
   double ExactVar[23][3];
   double ExactMean[23][3];

   *Mu = 0.0;
   *Sigma = -1.0;
   if (n <= 30 && n >= 8 && L <= 5 && L >= 3) {
      for (i = 8; i <= 30; i++) {
         for (j = 3; j <= 5; j++) {
            ExactVar[i - 8][j - 3] = -1.0;
         }
      }
      ExactMean[0][0]  = 2.29977;     ExactVar[0][0]   = 1.86729e-1;
      ExactMean[8][1]  = 3.23872;     ExactVar[8][1]   = 1.00739e-1;
      ExactMean[12][2] = 3.817;       ExactVar[12][2]  = 8.1539e-2;
      ExactMean[17][2] = 4.01429;     ExactVar[17][2]  = 6.9463e-2;
      ExactMean[22][2] = 4.160005;    ExactVar[22][2]  = 5.91489e-2;
      if (ExactVar[n - 8][L - 3] > 0.0) {
         *Mu = ExactMean[n - 8][L - 3];
         *Sigma = sqrt (ExactVar[n - 8][L - 3]);
      }
   }
}


/*-------------------------------------------------------------------------*/

static void WriteDataDisc (unif01_Gen * gen, char *TestName,
   long N, long n, int r, int s, int L, double Mu, double Sigma)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   s = %1d,   L = %1d\n", s, L);

   if (Sigma > 0.0) {
      printf ("   Mu    = ");
      num_WriteD (Mu, 15, 7, 7);
      printf ("\n   Sigma = ");
      num_WriteD (Sigma, 15, 7, 7);
      printf ("\n");
   } else {
      printf ("   Mu and Sigma are unknown.\n");
   }
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

static void WriteResultsDiscOver (sentrop_Res *res, double NLR, double Sum2,
   double SumSq, double Mu, double Sigma, double Mean, double Var,
   double Corr)
{
   printf ("\nEmpirical mean                        :");
   num_WriteD (Mean, 12, 8, 2);
   printf ("\nEmpirical variance                    :");
   num_WriteD (Var, 12, 8, 2);
   printf ("\n");

#ifndef STABLE
   if (swrite_Basic) {
      printf ("\nInformation on the numerical error in the computation\n"
              "   of the variance and the correlation:\n"
              "  Mean * Mean                         :");
      num_WriteD (Mean * Mean, 20, 16, 2);
      printf ("\n  Sum2 / N                            :");
      num_WriteD (Sum2 / NLR, 20, 16, 2);
      printf ("\n  SumSq / (N-1)                       :");
      num_WriteD (SumSq / (NLR - 1.0), 20, 16, 2);
      printf ("\n\n");
   }
#endif

   if (Sigma > 0.0) {
      printf ("\nDeviation from the theoretical mean   :");
      num_WriteD (Mean - Mu, 12, 8, 2);
      printf ("\nStandardized standard deviation       :");
      gofw_Writep2 (res->Bas->sVal2[gofw_Mean], res->Bas->pVal2[gofw_Mean]);
      printf ("\n");
   } else {
      printf ("\n\n");
   }
   printf ("Empirical correlation                 :");
   num_WriteD (Corr, 12, 8, 2);
   printf ("\nEmpirical correlation * sqrt(N)       :");
   gofw_Writep2 (res->Bas->sVal2[gofw_Cor], res->Bas->pVal2[gofw_Cor]);
   printf ("\n");
}


/*-------------------------------------------------------------------------*/

void sentrop_EntropyDiscOver (unif01_Gen * gen, sentrop_Res * res,
   long N, long n, int r, int s, int L)
{

   long i;                        /* Index */
   unsigned long Block1, Block0;  /* Blocks of bits */
   long Seq;                      /* Replication number */
   double Entropy;                /* Value of the entropy S */
   double tempPrev;               /* Previous value of entropy */
   double SumSq;                  /* To compute the covariance */
   double Corr;                   /* Empirical correlation */
   double Var;                    /* Empirical variance */
   double Mean;                   /* Empirical mean */
   double Sigma, Mu;              /* Parameters of the normal law */
   double Sum2, Sum;              /* Temporary variables */
   long d;                        /* 2^s */
   long C;                        /* 2^L */
   long nSurs;                    /* n / s */
   double xLgx[NLIM + 1];         /* = -i/n * Lg (i/n) */
   double NLR = N;
   double temp, E1;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sentrop_EntropyDiscOver test";

   Timer = chrono_Create ();
   InitExactOver (n, L, &Mu, &Sigma);
   if (swrite_Basic)
      WriteDataDisc (gen, TestName, N, n, r, s, L, Mu, Sigma);

   util_Assert (L <= n - L, "sentrop_EntropyDiscOver:   L > n-L");
   util_Assert (n <= 31, "sentrop_EntropyDiscOver:   n > 31");
   util_Assert (r <= 31, "sentrop_EntropyDiscOver:   r > 31");
   util_Assert (s <= 31, "sentrop_EntropyDiscOver:   s > 31");
   util_Assert (n % s == 0, "sentrop_EntropyDiscOver:   n % s != 0");
   util_Assert (N > 1, "sentrop_EntropyDiscOver:   N <= 1");

   d = num_TwoExp[s];
   C = num_TwoExp[L];
   nSurs = n / s;

   if (res == NULL) {
      localRes = TRUE;
      res = sentrop_CreateRes ();
   }
   InitRes (res, N, C - 1, "sentrop_EntropyDiscOver");
   CalcLgx (xLgx, n);
   tempPrev = SumSq = Sum2 = Sum = 0.0;

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i < C; i++)
         res->Count[i] = 0;

      Block0 = unif01_StripB (gen, r, s);
      for (i = 2; i <= nSurs; i++)
         Block0 = Block0 * d + unif01_StripB (gen, r, s);

      /* Compute entropy of the block of n bits = Block0. */
      /* This block has less than 31 bits. */
      Block1 = Block0;
      for (i = 0; i <= n - L - 1; i++) {
         ++res->Count[Block1 % C];
         Block1 >>= 1;
      }
      Block1 = (Block1 % C) + C * (Block0 % C);
      for (i = n - L; i < n; i++) {
         ++res->Count[Block1 % C];
         Block1 >>= 1;
      }

      Entropy = 0.0;
      for (i = 0; i < C; i++) {
         util_Assert (res->Count[i] <= NLIM,
            "sentrop_EntropyDiscOver:   NLIM is too small");
         Entropy += xLgx[res->Count[i]];
      }

#ifdef STABLE
      /* Ideally, we should use the moving average for numerical stability.
         But we shall use the first observed value instead; it should be
         typical and will prevent loss of precision (unless it is 0). */
      if (1 == Seq)
         E1 = Entropy;
      temp = Entropy - E1;
      Sum += temp;
      Sum2 += temp * temp;
      SumSq += temp * tempPrev;
      tempPrev = temp;

#else
      /* The naive method: it is simple but numerically unstable. It can be
         used for debugging and testing the more stable calculation in the
         case of small samples. */
      Sum += Entropy;
      Sum2 += Entropy * Entropy;
      SumSq += Entropy * tempPrev;
      tempPrev = Entropy;
#endif

      if (swrite_Counters)
         tables_WriteTabL (res->Count, 0, C - 1, 5, 10, "Counters:");

      if (swrite_Collectors) {
         printf ("Entropy = ");
         num_WriteD (Entropy, 15, 6, 1);
         printf ("\n");
      }
   }

   /* We now test the correlation between successive values of the entropy.
      Corr should have mean 0 and variance 1. We use a numerically stable
      calculation. */

#ifdef STABLE
   Mean = Sum / NLR + E1;
   Var = Sum2 / NLR - (E1 - Mean) * (E1 - Mean);
   Var *= NLR / (NLR - 1.0);
   temp = (Entropy + E1 * NLR - 2.0 * NLR * Mean) * E1 / (NLR - 1.0);
   Corr = SumSq / (NLR - 1.0) - temp - Mean * Mean;
   if (Var <= 0.0) {
      Corr = 1.0e100;
      util_Warning (TRUE,
         "Empirical variance <= 0.   Correlation set to 1e100.");
   } else
      Corr /= Var;

#else
   /* Naive calculations. Here, there could be huge losses of precision
      because Mean*Mean, Sum2/NLR, and SumSq/(NLR - 1.0) may be very close. */
   Mean = Sum / NLR;
   Var = (Sum2 / NLR - Mean * Mean) * NLR / (NLR - 1.0);
   Corr = (SumSq / (NLR - 1.0) - Mean * Mean) / Var;

#endif


   if (Sigma > 0.0) {
      /* We know the true values of Mu and Sigma */
      res->Bas->sVal2[gofw_Mean] = (Mean - Mu) * sqrt (NLR) / Sigma;
      res->Bas->pVal2[gofw_Mean] = fbar_Normal1 (res->Bas->sVal2[gofw_Mean]);
   } else
      res->Bas->pVal2[gofw_Mean] = -1.0;

   res->Bas->sVal2[gofw_Cor] = Corr * sqrt (NLR);
   res->Bas->pVal2[gofw_Cor] = fbar_Normal1 (res->Bas->sVal2[gofw_Cor]);

   if (swrite_Basic) {
      WriteResultsDiscOver (res, NLR, Sum2, SumSq, Mu, Sigma, Mean, Var,
         Corr);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sentrop_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sentrop_EntropyDiscOver2 (unif01_Gen * gen, sentrop_Res * res,
   long N, long n, int r, int s, int L)
{
   long i, j;                     /* Indices */
   unsigned long B2, B1, B0;      /* Blocks of bits */
   long Seq;                      /* Replication number */
   double Entropy;                /* Value of the entropy S */
   double tempPrev;               /* Previous value of the entropy */
   double SumSq;                  /* To compute the covariance */
   double Corr;                   /* Empirical correlation */
   double Var;                    /* Empirical variance */
   double Mean;                   /* Empirical mean */
   double Sigma, Mu;              /* Parameters of the normal law */
   double Sum2, Sum;              /* Temporary variables */
   unsigned long d;               /* 2^s */
   long C;                        /* 2^L */
   unsigned long CLC;             /* 2^L */
   long m0;                       /* m0 = ceil (L/s) */
   long m;                        /* m = n/s */
   double xLgx[NLIM + 1];         /* = -i/n * Lg (i/n) */
   double NLR = N;
   double temp, E1;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sentrop_EntropyDiscOver2 test";

   Timer = chrono_Create ();
   InitExactOver (n, L, &Mu, &Sigma);
   if (swrite_Basic)
      WriteDataDisc (gen, TestName, N, n, r, s, L, Mu, Sigma);

   util_Assert (L <= n, "sentrop_EntropyDiscOver2:   L > n");
   util_Assert (L <= 15, "sentrop_EntropyDiscOver2:   L > 15");
   util_Assert (r <= 31, "sentrop_EntropyDiscOver2:   r > 31");
   util_Assert (s <= 31, "sentrop_EntropyDiscOver2:   s > 31");
   util_Assert (L + s <= 31, "sentrop_EntropyDiscOver2:   L+s > 31");
   util_Assert (n % s == 0, "sentrop_EntropyDiscOver2:   n % s != 0");

   d = num_TwoExp[s];
   m = n / s;
   m0 = L / s;
   if (m0 * s < L)
      ++m0;
   /* B0 must not be larger than LONG_MAX (31 bits) */
   util_Assert (m0 * s <= 31, "sentrop_EntropyDiscOver2:   m0 * s > 31");
   C = num_TwoExp[L];
   CLC = num_TwoExp[L];

   if (res == NULL) {
      localRes = TRUE;
      res = sentrop_CreateRes ();
   }
   InitRes (res, N, C - 1, "sentrop_EntropyDiscOver2");
   tempPrev = SumSq = Sum2 = Sum = 0.0;
   CalcLgx (xLgx, n);

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i < C; i++)
         res->Count[i] = 0;

      B0 = unif01_StripB (gen, r, s);
      for (j = 2; j <= m0; j++)
         B0 = B0 * d + unif01_StripB (gen, r, s);

      /* B0 now contains the bits 0,...,0,b_1,...,b_{m0*s} */
      B2 = B0;

      /* Count the blocks of L bits in b_1,...,b_{m0*s} */
      for (i = 0; i <= m0 * s - L; i++) {
         ++res->Count[B2 % CLC];
         B2 >>= 1;
      }
      B1 = B0 % CLC;
      B0 = B2 % CLC;
      /* B1 contains 0,...,0,b_{m0*s-L+1},...,b_{m0*s} */
      /* B0 contains 0,...,0,b_1,...,b_{L-1} */
      for (j = 1; j <= m - m0; j++) {
         B1 = B1 * d + unif01_StripB (gen, r, s);
         B2 = B1;
         B1 %= CLC;
         /* B1 and B2 contain L bits and L+s bits, resp. */
         for (i = 1; i <= s; i++) {
            ++res->Count[B2 % CLC];
            B2 >>= 1;
         }
      }

      /* B1 contains 0,...,0,b_{m*s-L+1},...,b_{m*s}. */
      /* Her we must have 2 * L <= 31. */
      B2 = B0 + B1 * (CLC / 2);
      /* B2 contains 0,..,0,b_{m*s-L+1},..,b_{m*s},b_1,...,b_{L-1}. */
      /* Now count blocks with overlap. */
      for (i = 1; i < L; i++) {
         ++res->Count[B2 % CLC];
         B2 >>= 1;
      }

      /* Compute entropy */
      Entropy = 0.0;
      for (i = 0; i < C; i++) {
         util_Assert (res->Count[i] <= NLIM,
            "sentrop_EntropyDiscOver2:   NLIM is too small");
         Entropy += xLgx[res->Count[i]];
      }

#ifdef STABLE
      /* Ideally, we should use the moving average for numerical stability.
         But we shall use the first observed value of instead; it should be
         typical and will prevent loss of precision (unless it is 0). */

      if (1 == Seq)
         E1 = Entropy;
      temp = Entropy - E1;
      Sum += temp;
      Sum2 += temp * temp;
      SumSq += temp * tempPrev;
      tempPrev = temp;

#else
      /* The naive unstable method */
      Sum += Entropy;
      Sum2 += Entropy * Entropy;
      SumSq += Entropy * tempPrev;
      tempPrev = Entropy;
#endif

      if (swrite_Counters)
         tables_WriteTabL (res->Count, 0, C - 1, 5, 10, "Counters:");

      if (swrite_Collectors) {
         printf ("Entropy = ");
         num_WriteD (Entropy, 15, 6, 1);
         printf ("\n");
      }
   }

   /* We now test the correlation between successive values of the */
   /* entropy. Corr should have mean 0 and variance 1. */

#ifdef STABLE
   Mean = Sum / NLR + E1;
   Var = Sum2 / NLR - (E1 - Mean) * (E1 - Mean);
   Var *= NLR / (NLR - 1.0);
   temp = (Entropy + E1 * NLR - 2.0 * NLR * Mean) * E1 / (NLR - 1.0);
   Corr = SumSq / (NLR - 1.0) - temp - Mean * Mean;
   if (Var <= 0.0) {
      Corr = 1.0e100;
      util_Warning (TRUE,
         "Empirical variance <= 0.   Correlation set to 1e100.");
   } else
      Corr /= Var;

#else
   /* Naive calculations. Here, there could be huge losses of precision
      because Mean*Mean, Sum2/NLR, and SumSq/(NLR - 1.0) may be very close. */
   Mean = Sum / NLR;
   Var = (Sum2 / NLR - Mean * Mean) * NLR / (NLR - 1.0);
   Corr = (SumSq / (NLR - 1.0) - Mean * Mean) / Var;
#endif

   if (Sigma > 0.0) {
      /* We know the true values of Mu and Sigma */
      res->Bas->sVal2[gofw_Mean] = (Mean - Mu) * sqrt (NLR) / Sigma;
      res->Bas->pVal2[gofw_Mean] = fbar_Normal1 (res->Bas->sVal2[gofw_Mean]);
   } else
      res->Bas->sVal2[gofw_Mean] = -1.0;

   res->Bas->sVal2[gofw_Cor] = Corr * sqrt (NLR);
   res->Bas->pVal2[gofw_Cor] = fbar_Normal1 (res->Bas->sVal2[gofw_Cor]);

   if (swrite_Basic) {
      WriteResultsDiscOver (res, NLR, Sum2, SumSq, Mu, Sigma, Mean, Var,
         Corr);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sentrop_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static void WriteResultsDisc (long N, gofw_TestArray sVal2, gofw_TestArray pVal2,
    sres_Basic *res)
{
   if (N > 1) {
      gofw_WriteActiveTests0 (N, sVal2, pVal2);
      swrite_NormalSumTest (N, res);
      printf ("Standardized empirical correlation    :");
      gofw_Writep2 (sVal2[gofw_Cor], pVal2[gofw_Cor]);
   } else {
      printf ("Standardized statistic value          :");
      gofw_Writep2 (sVal2[gofw_Mean], pVal2[gofw_Mean]);
   }
}


/*-------------------------------------------------------------------------*/

static void EntropyDisc00 (unif01_Gen * gen, sentrop_Res * res,
   long N, long n, int r, int s, int L)
/*
 * Test based on the discrete entropy, proposed by Compagner and L'Ecuyer
 */
{
   long Seq;
   long j;
   long i;
   double EntropyNorm;            /* Normalized entropy */
   double Entropy;                /* Value of entropy S */
   double EntropyPrev;            /* Previous value of entropy */
   double SumSq;                  /* To compute the covariance */
   double Sigma, Mu;              /* Parameters of the normal law */
   double tem;
   double nLR = n;
   long d;                        /* 2^s */
   long C;                        /* 2^L */
   long LSurs;                    /* L / s */
   long sSurL;                    /* s / L */
   long nLSurs;                   /* nL / s */
   double xLgx[NLIM + 1];         /* = -i/n * Lg (i/n) */
   unsigned long Block;
   unsigned long Number;
   unsigned int LL = L;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sentrop_EntropyDisc test";

   Timer = chrono_Create ();
   if (s <= L && L % s) {
      util_Error ("EntropyDisc00:   s <= L and L % s != 0");
   }
   if (s > L && s % L) {
      util_Error ("EntropyDisc00:   s > L and s % L != 0");
   }

   d = num_TwoExp[s];
   C = num_TwoExp[L];

   if (s <= L)
      LSurs = L / s;
   else {
      sSurL = s / L;
      nLSurs = n / sSurL;
      if (n % sSurL)
         ++nLSurs;
   }

   util_Assert (n / num_TwoExp[L] < NLIM,
      "sentrop_EntropyDisc:    n/2^L is too large");
   smultin_MultinomMuSigma (n, num_TwoExp[L], 0.0, (double) n,
        FoncMNEntropie, &Mu, &Sigma);

   if (swrite_Basic)
      WriteDataDisc (gen, TestName, N, n, r, s, L, Mu, Sigma);

   if (res == NULL) {
      localRes = TRUE;
      res = sentrop_CreateRes ();
   }
   InitRes (res, N, C - 1, "sentrop_EntropyDisc");
   CalcLgx (xLgx, n);

   statcoll_SetDesc (res->Bas->sVal1, "EntropyDisc sVal1");
   statcoll_SetDesc (res->Bas->pVal1, "EntropyDisc pVal1");
   SumSq = EntropyPrev = 0.0;

   for (Seq = 1; Seq <= N; Seq++) {
      for (i = 0; i < C; i++)
         res->Count[i] = 0;

      if (s <= L) {
         for (i = 1; i <= n; i++) {
            Block = unif01_StripB (gen, r, s);
            for (j = 2; j <= LSurs; j++)
               Block = Block * d + unif01_StripB (gen, r, s);
            ++res->Count[Block];
         }

      } else {                    /* s > L */
         for (i = 1; i <= nLSurs; i++) {
            Number = unif01_StripB (gen, r, s);
            for (j = 1; j <= sSurL; j++) {
               Block = Number % C;
               ++res->Count[Block];
               Number >>= LL;
            }
         }
      }

      /* Compute entropy */
      Entropy = 0.0;
      for (i = 0; i < C; i++) {
         if (res->Count[i] > NLIM) {
            tem = res->Count[i] / nLR;
            tem *= -num_Log2 (tem);
            Entropy += tem;
         } else if (res->Count[i] > 0) {
            Entropy += xLgx[res->Count[i]];
         }
      }
      EntropyNorm = (Entropy - Mu) / Sigma;
      statcoll_AddObs (res->Bas->sVal1, EntropyNorm);
      SumSq += EntropyNorm * EntropyPrev;
      EntropyPrev = EntropyNorm;

      if (swrite_Counters)
         tables_WriteTabL (res->Count, 0, C - 1, 5, 10, "Counters:");

      if (swrite_Collectors) {
         printf ("Entropy = ");
         num_WriteD (Entropy, 15, 6, 1);
         printf ("\n");
      }
   }

   gofw_ActiveTests2 (res->Bas->sVal1->V, res->Bas->pVal1->V, N, wdist_Normal,
      (double *) NULL, res->Bas->sVal2, res->Bas->pVal2);
   res->Bas->pVal1->NObs = N;
   sres_GetNormalSumStat (res->Bas);

   /* We now test the correlation between successive values of the entropy.
      The next SumSq should have mean 0 and variance 1. */
   if (N > 1) {
      res->Bas->sVal2[gofw_Cor] = SumSq / sqrt ((double) N);
      res->Bas->pVal2[gofw_Cor] = fbar_Normal1 (res->Bas->sVal2[gofw_Cor]);
   }
   if (swrite_Collectors) {
      statcoll_Write (res->Bas->sVal1, 5, 14, 4, 3);
   }
   if (swrite_Basic) {
      WriteResultsDisc (N, res->Bas->sVal2, res->Bas->pVal2, res->Bas);
      swrite_Final (gen, Timer);
   }

   if (localRes)
      sentrop_DeleteRes (res);
   chrono_Delete (Timer);
}


/*-------------------------------------------------------------------------*/

static void WriteDataEnt (long N, long n, int r, int s, int L)
{
   printf ("***********************************************************\n"
           "Test sentrop_EntropyDisc calling smultin_Multinomial\n\n");
   printf ("   N = %2ld,  n = %8ld,  r = %2d", N, n, r);
   printf (",   s = %1d,   L = %1d\n\n", s, L);
}


/*-------------------------------------------------------------------------*/

void sentrop_EntropyDisc (unif01_Gen * gen, sentrop_Res * res,
   long N, long n, int r, int s, int L)
/*
 * Use smultin_Multinomial to replace sentrop_EntropyDisc
 */
{
   int i;
   int t;
   long d;
   double x;
   smultin_Param *par;
   double ValDelta[] = { 0.0 };

   if (s > L) {
      EntropyDisc00 (gen, res, N, n, r, s, L);
      return;                     /* <--- Case "s > L" end up here.  */
   }

   if (swrite_Basic)
      WriteDataEnt (N, n, r, s, L);
   util_Assert (L % s == 0, "sentrop_EntropyDisc:   L % s != 0");

   d = num_TwoExp[s];
   t = L / s;

   x = d;
   for (i = 2; i <= t; i++)
      x *= d;

   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, 3);
   if (NULL == res) {
      if (n / x <= 8.0)
         smultin_Multinomial (gen, par, NULL, N, n, r, d, t, TRUE);
      else
         smultin_Multinomial (gen, par, NULL, N, n, r, d, t, FALSE);
   } else {
      smultin_Res *resm;
      resm = smultin_CreateRes (par);
      if (n / x <= 8.0)
         smultin_Multinomial (gen, par, resm, N, n, r, d, t, TRUE);
      else
         smultin_Multinomial (gen, par, resm, N, n, r, d, t, FALSE);
      InitRes (res, N, 0, "sentrop_EntropyDisc");
      statcoll_SetDesc (res->Bas->sVal1, "EntropyDisc sVal1");
      res->Bas->sVal1->NObs = resm->Collector[0]->NObs;
      tables_CopyTabD (resm->Collector[0]->V, res->Bas->sVal1->V, 1, N);
      tables_CopyTabD (resm->sVal2[0], res->Bas->sVal2, 0, gofw_NTestTypes - 1);
      tables_CopyTabD (resm->pVal2[0], res->Bas->pVal2, 0, gofw_NTestTypes - 1);
      smultin_DeleteRes (resm);
   }
   smultin_DeleteParam (par);
}


/*=========================================================================*/

static void WriteDataDM (unif01_Gen * gen, char *TestName,
   long N, long n, int r, long m)
{
   swrite_Head (gen, TestName, N, n, r);
   printf (",   m = %1ld\n\n", m);
}


/*-------------------------------------------------------------------------*/

void sentrop_EntropyDM (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r, long m)
{
   long i;                        /* Index */
   double I0, x;
   long Seq;                      /* Replication number */
   double Entropy;                /* Value of the statistic H(m, n) */
   double LnEntropy;
   double SumR;                   /* 1 + 1/2 + ... + 1/(2m-1) */
   double nLR = n;
   double Twom;                   /* 2m */
   double *AU;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sentrop_EntropyDM test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataDM (gen, TestName, N, n, r, m);

   Twom = 2 * m;
   I0 = Twom - 1.0;
   SumR = 0.0;
   for (i = 2 * m - 1; i >= 1; i--) {
      SumR += 1.0 / I0;
      I0 -= 1.0;
   }

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "sentrop_EntropyDM");
   AU = util_Calloc ((size_t) (n + 1), sizeof (double));
   statcoll_SetDesc (res->sVal1,
      "The N statistic values (a standard normal)");

   for (Seq = 1; Seq <= N; Seq++) {

      /* Generate the sample and sort */
      for (i = 1; i <= n; i++)
         AU[i] = unif01_StripD (gen, r);
      tables_QuickSortD (AU, 1, n);

      /* Compute empirical entropy */
      Entropy = 1.0;
      LnEntropy = 0.0;
      for (i = 1; i <= n; i++) {
         if (i - m < 1) {
            Entropy *= (AU[i + m] - AU[1]);
         } else if (i + m > n) {
            Entropy *= (AU[n] - AU[i - m]);
         } else
            Entropy *= (AU[i + m] - AU[i - m]);

         /* Avoid underflow */
         if (Entropy < Epsilon) {
            LnEntropy += log (Entropy);
            Entropy = 1.0;
         }
      }

      Entropy = (LnEntropy + log (Entropy)) / nLR + log (nLR / Twom);

      /* Compute standardized statistic */
      x  =  sqrt (3.0 * Twom * nLR) * (Entropy + log (Twom) + Euler - SumR);
      statcoll_AddObs (res->sVal1, x);
   }

   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_Normal,
      (double *) NULL, res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetNormalSumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2,
         "Normal statistic                      :");
      swrite_NormalSumTest (N, res);
      swrite_Final (gen, Timer);
   }

   util_Free (AU);
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sentrop_EntropyDMCirc (unif01_Gen * gen, sres_Basic * res,
   long N, long n, int r, long m)
{
   long i;                        /* Index */
   double I0, x;
   long Seq;                      /* Replication number */
   double Entropy;                /* Value of the statistic H(m, n) */
   double LnEntropy;
   double SumR;                   /* 1 + 1/2 + ... + 1/(2m-1) */
   double nLR = n;
   double Twom;                   /* 2m */
   double *AU;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sentrop_EntropyDMCirc test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataDM (gen, TestName, N, n, r, m);

   Twom = 2 * m;
   I0 = Twom - 1.0;
   SumR = 0.0;
   for (i = 2 * m - 1; i >= 1; i--) {
      SumR += 1.0 / I0;
      I0 -= 1.0;
   }

   if (res == NULL) {
      localRes = TRUE;
      res = sres_CreateBasic ();
   }
   sres_InitBasic (res, N, "sentrop_EntropyDMCirc");
   AU = util_Calloc ((size_t) (n + 1), sizeof (double));
   statcoll_SetDesc (res->sVal1,
      "The N statistic values (a standard normal)");

   for (Seq = 1; Seq <= N; Seq++) {
      /* Generate the sample and sort */
      for (i = 1; i <= n; i++)
         AU[i] = unif01_StripD (gen, r);
      tables_QuickSortD (AU, 1, n);

      /* Compute empirical entropy */
      Entropy = 1.0;
      LnEntropy = 0.0;
      for (i = 1; i <= n; i++) {
         if (i - m < 1) {
            Entropy *= (AU[i + m] - AU[n + i - m] + 1.0);
         } else if (i + m > n) {
            Entropy *= (AU[i + m - n] - AU[i - m] + 1.0);
         } else
            Entropy *= (AU[i + m] - AU[i - m]);

         if (Entropy < Epsilon) {
            LnEntropy += log (Entropy);
            Entropy = 1.0;
         }
      }

      Entropy = (LnEntropy + log (Entropy)) / nLR + log (nLR / Twom);

      /* Compute standardized statistic */
      x  = sqrt (3.0 * Twom * nLR) * ((Entropy + log (Twom) + Euler) - SumR);
      statcoll_AddObs (res->sVal1, x);
   }

   gofw_ActiveTests2 (res->sVal1->V, res->pVal1->V, N, wdist_Normal,
                     (double *) NULL, res->sVal2, res->pVal2);
   res->pVal1->NObs = N;
   sres_GetNormalSumStat (res);

   if (swrite_Collectors)
      statcoll_Write (res->sVal1, 5, 14, 4, 3);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->sVal2, res->pVal2,
         "Normal statistic                      :");
      swrite_NormalSumTest (N, res);
      swrite_Final (gen, Timer);
   }

   util_Free (AU);
   if (localRes)
      sres_DeleteBasic (res);
   chrono_Delete (Timer);
}
