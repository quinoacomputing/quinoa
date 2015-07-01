/*************************************************************************\
 *
 * Package:        TestU01
 * File:           sspectral.c
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

#include "util.h"
#include "num.h"
#include "chrono.h"
#include "tables.h"

#include "sspectral.h"
#include "swrite.h"
#include "wdist.h"
#include "unif01.h"

#include "gofw.h"
#include "statcoll.h"

#include "fftc.c"

#include <math.h>
#include <stdlib.h>
#include <string.h>






/*-------------------------------- Functions ------------------------------*/


static void InitRes (
   sspectral_Res *res,
   long N,
   long jmin,
   long jmax,
   char *nam
)
/* 
 * Initializes the sspectral_Res structure
 */
{
   long j;
   sres_InitBasic (res->Bas, N, nam);
   if (jmax > res->jmax)
      res->Coef = util_Realloc (res->Coef, (jmax + 200) * sizeof (double));
   for (j = 0; j <= jmax; j++)
      res->Coef[j] = 0.0;
   res->jmin = jmin;
   res->jmax = jmax;
   res->Bas->name = util_Realloc (res->Bas->name,
                                  1 + strlen (nam) * sizeof (char));
   strcpy (res->Bas->name, nam);
}


/*-------------------------------------------------------------------------*/

sspectral_Res * sspectral_CreateRes (void)
{
   sspectral_Res *res;
   res = util_Malloc (sizeof (sspectral_Res));
   res->Bas = sres_CreateBasic ();
   res->Coef = util_Calloc (1, sizeof (double));
   res->jmax = 0;
   return res;
}


/*-------------------------------------------------------------------------*/

void sspectral_DeleteRes (sspectral_Res *res)
{
   if (res == NULL)
      return;
   sres_DeleteBasic (res->Bas);
   util_Free (res->Coef);
   util_Free (res);
}


/*=========================================================================*/

static void WriteDataFour (
   unif01_Gen *gen,      /* generator */
   char *Test,           /* Test name */
   long N,               /* Number of replications */
   int k,                /* Sample size n = 2^k */
   int r,                /* r first bits of each random number dropped */
   int s                 /* s bits of each random number used */
)
{
   long n;
   n = num_TwoExp[k];
   swrite_Head (gen, Test, N, n, r);
   printf (",   s = %4d,   k = %4d\n\n", s, k);
}


/*-------------------------------------------------------------------------*/

void sspectral_Fourier1 (unif01_Gen *gen, sspectral_Res *res,
   long N, int t, int r, int s)
{
   const unsigned long SBIT = 1UL << (s - 1);
   unsigned long jBit;
   unsigned long Z;
   long k, KALL, Seq, n, i;
   double x, NbExp, h, per;
   long co;
   double *A;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sspectral_Fourier1 test";

   Timer = chrono_Create ();
   util_Assert (t <= 20, "sspectral_Fourier1:   k > 20");
   util_Assert (t > 1, "sspectral_Fourier1:   k < 2");
   if (swrite_Basic)
      WriteDataFour (gen, TestName, N, t, r, s);
   if (res == NULL) {
      localRes = TRUE;
      res = sspectral_CreateRes ();
   }
   n = num_TwoExp[t];
   KALL = n / s;
   if (n % s > 0)
      KALL++;
   per = 0.95;
   NbExp = per * (n / 2 + 1);
/*   h = 3.0 * n; */
   h = 2.995732274 * n;
   InitRes (res, N, 0, n, "sspectral_Fourier1");
   statcoll_SetDesc (res->Bas->sVal1, "sVal1:   a standard normal");
   A = res->Coef;

   for (Seq = 1; Seq <= N; Seq++) {
      /* Fill array A: 1 for bit 1, -1 for bit 0 */
      i = 0;
      for (k = 0; k < KALL; k++) {
         Z = unif01_StripB (gen, r, s);
         jBit = SBIT;
         while (jBit) {
            if (jBit & Z)
               A[i] = 1.0;
            else
               A[i] = -1.0;
            jBit >>= 1;
            i++;
         }
      }
      /* 
       * Compute the Fourier transform of A and return the result in A. The
       * first half of the array, (from 0 to n/2) is filled with the real
       * components of the FFT. The second half of the array (from n/2+1 to
       * n-1) is filled with the imaginary components of the FFT.
       * The n new elements of A are thus:
       *      [Re(0), Re(1), ...., Re(n/2), Im(n/2-1), ..., Im(1)]
       * The procedure is due to H.V. Sorensen, University of Pennsylvania 
       * and is found in file fftc.c.
       */
      rsrfft (A, t);

      /* Count the number of Fourier coefficients smaller than h */
      co = 0;
      for (i = 1; i < n / 2; i++) {
         x = A[i] * A[i] + A[n - i] * A[n - i];
         if (x < h)
            co++;
      }
      if (A[0] * A[0] < h)
         co++;

      /* Compute the NIST statistic */
      x = (co - NbExp) / sqrt (NbExp * (1.0 - per));
      statcoll_AddObs (res->Bas->sVal1, x);

      if (swrite_Counters) {
         tables_WriteTabD (res->Coef, 0, n - 1, 5, 14, 5, 5,
            "Fourier coefficients");
      }
   }

   gofw_ActiveTests2 (res->Bas->sVal1->V, res->Bas->pVal1->V, N, wdist_Normal,
      (double *) NULL, res->Bas->sVal2, res->Bas->pVal2);
   res->Bas->pVal1->NObs = N;
   sres_GetNormalSumStat (res->Bas);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->Bas->sVal2, res->Bas->pVal2,
         "Normal statistic                      :");
      swrite_NormalSumTest (N, res->Bas);
      if (swrite_Collectors)
         statcoll_Write (res->Bas->sVal1, 5, 14, 4, 3);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sspectral_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sspectral_Fourier3 (unif01_Gen *gen, sspectral_Res *res,
   long N, int t, int r, int s)
{
   const unsigned long SBIT = 1UL << (s - 1);
   unsigned long jBit;
   unsigned long Z;
   long k, KALL, Seq, n, i;
   double *A, *B;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sspectral_Fourier3 test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataFour (gen, TestName, N, t, r, s);
   util_Assert (r + s <= 32, "sspectral_Fourier3:   r + s > 32");
   util_Assert (t <= 26, "sspectral_Fourier3:   k > 26");
   util_Assert (t >= 2, "sspectral_Fourier3:   k < 2");
   if (res == NULL) {
      localRes = TRUE;
      res = sspectral_CreateRes ();
   }
   n = num_TwoExp[t];
   KALL = n / s + 1;
   InitRes (res, n/4 + 1, 0, n, "sspectral_Fourier3");
   statcoll_SetDesc (res->Bas->sVal1, "sVal1:   a standard normal");
   B = res->Bas->sVal1->V;
   A = res->Coef;
   for (i = 0; i <= n / 4; i++)
      B[i] = 0.0;

   for (Seq = 1; Seq <= N; Seq++) {
      /* Fill array A: 1 for bit 1, -1 for bit 0 */
      i = 0;
      for (k = 0; k < KALL; k++) {
         Z = unif01_StripB (gen, r, s);
         jBit = SBIT;
         while (jBit) {
            if (jBit & Z)
               A[i] = 1.0;
            else
               A[i] = -1.0;
            jBit >>= 1;
            i++;
         }
      }
      /* 
       * Compute the Fourier transform of A and return the result in A. The
       * first half of the array, (from 0 to n/2) is filled with the real
       * components of the FFT. The second half of the array (from n/2+1 to
       * n-1) is filled with the imaginary components of the FFT.
       * The n new elements of A are thus:
       *      [Re(0), Re(1), ...., Re(n/2), Im(n/2-1), ..., Im(1)]
       * The procedure is due to H.V. Sorensen, University of Pennsylvania 
       * and is found in file fftc.c.
       */
      rsrfft (A, t);

      /* Add the squares of the Fourier coefficients over the N replications
         for each i = [1, ..., n/4], and keep them in B[i] */
      for (i = 1; i <= n / 4; i++)
         B[i] += A[i] * A[i] + A[n - i] * A[n - i];

      if (0 && swrite_Counters)
	 tables_WriteTabD (B, 1, n / 4, 5, 14, 5, 5,
	     "Sums of square of Fourier coefficients");
   }

   /* There is an extra sqrt (n) factor between the Fourier coefficients
      of Sorensen and those of Erdmann */
   for (i = 1; i <= n / 4; i++)
      B[i] /= n;

   /* The N random variables have been added for each i and kept in B[i].
      Their mean (1) and variance (~1) is known from Diane Erdmann. Now
      consider the B[i] as n/4 normal random variables. */
   for (i = 1; i <= n / 4; i++) {
      B[i] = (B[i] - N) / sqrt (N * (1.0 - 2.0 / n));
      statcoll_AddObs (res->Bas->sVal1, B[i]);
   }

   gofw_ActiveTests2 (res->Bas->sVal1->V, res->Bas->pVal1->V, n/4, wdist_Normal,
      (double *) NULL, res->Bas->sVal2, res->Bas->pVal2);
   res->Bas->pVal1->NObs = n/4;

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (n/4, res->Bas->sVal2, res->Bas->pVal2,
         "Normal statistic                      :");
      if (swrite_Collectors)
         statcoll_Write (res->Bas->sVal1, 5, 14, 4, 3);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sspectral_DeleteRes (res);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void sspectral_Fourier2 (unif01_Gen *gen, sspectral_Res *res,
   long N, int t, int r, int s)
{
   const unsigned long SBIT = 1UL << (s - 1);
   unsigned long jBit;
   unsigned long Z;
   long k, KALL, Seq, n, i;
   double *A;
   double x, sum;
   lebool localRes = FALSE;
   chrono_Chrono *Timer;
   char *TestName = "sspectral_Fourier2 test";

   Timer = chrono_Create ();
   if (swrite_Basic)
      WriteDataFour (gen, TestName, N, t, r, s);
   util_Assert (r + s <= 32, "sspectral_Fourier2:   r + s > 32");
   util_Assert (t <= 26, "sspectral_Fourier2:   k > 26");
   util_Assert (t >= 2, "sspectral_Fourier2:   k < 2");
   if (res == NULL) {
      localRes = TRUE;
      res = sspectral_CreateRes ();
   }
   n = num_TwoExp[t];
   KALL = n / s + 1;
   InitRes (res, N, 0, n, "sspectral_Fourier2");
   statcoll_SetDesc (res->Bas->sVal1, "sVal1:   a standard normal");
   A = res->Coef;

   for (Seq = 1; Seq <= N; Seq++) {
      /* Fill array A: 1 for bit 1, -1 for bit 0 */
      i = 0;
      for (k = 0; k < KALL; k++) {
         Z = unif01_StripB (gen, r, s);
         jBit = SBIT;
         while (jBit) {
            if (jBit & Z)
               A[i] = 1.0;
            else
               A[i] = -1.0;
            jBit >>= 1;
            i++;
         }
      }
      /* 
       * Compute the Fourier transform of A and return the result in A. The
       * first half of the array, (from 0 to n/2) is filled with the real
       * components of the FFT. The second half of the array (from n/2+1 to
       * n-1) is filled with the imaginary components of the FFT.
       * The n new elements of A are thus:
       *      [Re(0), Re(1), ...., Re(n/2), Im(n/2-1), ..., Im(1)]
       * The procedure is due to H.V. Sorensen, University of Pennsylvania 
       * and is found in file fftc.c.
       */
      rsrfft (A, t);

      /* Sum the square of the Fourier coefficients (only half of them) */
      sum = 0.0;
      for (i = 1; i <= n / 4; i++)
         sum += A[i] * A[i] + A[n - i] * A[n - i];

      /* There is an extra sqrt (n) factor between the Fourier coefficients
         of Sorensen and those of Erdmann */
      sum /= n;

      /* Standardize the statistic */
      x = 2.0*(sum - n / 4.0) / sqrt (n - 2.0);
      statcoll_AddObs (res->Bas->sVal1, x);

      if (swrite_Counters) {
         tables_WriteTabD (res->Coef, 0, n - 1, 5, 14, 5, 5,
            "Fourier coefficients");
      }
   }

   gofw_ActiveTests2 (res->Bas->sVal1->V, res->Bas->pVal1->V, N, wdist_Normal,
      (double *) NULL, res->Bas->sVal2, res->Bas->pVal2);
   res->Bas->pVal1->NObs = N;
   sres_GetNormalSumStat (res->Bas);

   if (swrite_Basic) {
      gofw_WriteActiveTests2 (N, res->Bas->sVal2, res->Bas->pVal2,
         "Normal statistic                      :");
      swrite_NormalSumTest (N, res->Bas);
      if (swrite_Collectors)
         statcoll_Write (res->Bas->sVal1, 5, 14, 4, 3);
      swrite_Final (gen, Timer);
   }
   if (localRes)
      sspectral_DeleteRes (res);
   chrono_Delete (Timer);
}

