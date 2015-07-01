/*************************************************************************\
 *
 * Package:        ProbDist
 * File:           gofs.c
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
#include "tables.h"
#include "num.h"
#include "num2.h"

#include "gofs.h"
#include "fdist.h"
#include "wdist.h"

#include <float.h>
#include <math.h>
#include <stdio.h>


#define TRACE0(x) printf ("***   " #x "     ");
#define TRACE1(N, form) printf ("***   " #N " = %"#form "     ", N);



/*---------------------------- extern variables ---------------------------*/

double gofs_MinExpected = 10.0;

double gofs_EpsilonAD = DBL_EPSILON / 2.0;



/*---------------------------- module variables ---------------------------*/

/* Used in discontinuous distributions */
static double EpsilonD = 1.0E-15;






/*-------------------------------- functions ------------------------------*/


void gofs_ContUnifTransform (double V[], long N, wdist_CFUNC F,
                             double par[], double U[])
{
   long i;
   for (i = 1; i <= N; i++)
      U[i] = F (par, V[i]);
}

/*-------------------------------------------------------------------------*/

void gofs_DiscUnifTransform (double V[], long N, wdist_DFUNC F,
                             fmass_INFO W, double U[])
{
   long i;
   for (i = 1; i <= N; i++)
      U[i] = F (W, (long) V[i]);
}

/*-------------------------------------------------------------------------*/

void gofs_DiffD (double U[], double D[], long N1, long N2, 
                 double a, double b)
{
   long i;
   D[N1 - 1] = U[N1] - a;
   for (i = N1; i < N2; i++)
      D[i] = U[i + 1] - U[i];
   D[N2] = b - U[N2];
}

/*-------------------------------------------------------------------------*/
#ifdef USE_LONGLONG

void gofs_DiffLL (longlong U[], longlong D[], long N1, long N2,
                  longlong a, longlong b)
{
   long i;
   D[N1 - 1] = U[N1] - a;
   for (i = N1; i < N2; i++)
      D[i] = U[i + 1] - U[i];
   D[N2] = b - U[N2];
}

/*-------------------------------------------------------------------------*/

void gofs_DiffULL (ulonglong U[], ulonglong D[], long N1, long N2,
                   ulonglong a, ulonglong b)
{
   long i;
   D[N1 - 1] = U[N1] - a;
   for (i = N1; i < N2; i++)
      D[i] = U[i + 1] - U[i];
   D[N2] = b - U[N2];
}

#endif
/*-------------------------------------------------------------------------*/

void gofs_DiffL (long U[], long D[], long N1, long N2, long a, long b)
{
   long i;
   D[N1 - 1] = U[N1] - a;
   for (i = N1; i < N2; i++)
      D[i] = U[i + 1] - U[i];
   D[N2] = b - U[N2];
}

/*-------------------------------------------------------------------------*/

void gofs_IterateSpacings (double V[], double S[], long N)
{
   long i;
   tables_QuickSortD (S, 0, N);
   for (i = 0; i < N; i++)
      S[N - i] = (i + 1) * (S[N - i] - S[N - i - 1]);
   S[0] = (N + 1) * S[0];
   V[1] = S[0];
   for (i = 2; i <= N; i++)
      V[i] = V[i - 1] + S[i - 1];
}

/*-------------------------------------------------------------------------*/

void gofs_PowerRatios (double U[], long N)
{
   long i;
   /* Assumes that the U[i] are already sorted in increasing order. */
   for (i = 1; i < N; i++) {
      if (U[i + 1] == 0.0 || U[i + 1] == -0.0) {
         /* util_Warning (1, "gofs_PowerRatios: 0 divisor"); */
         U[i] = 1.0;
      } else
         U[i] = pow (U[i] / U[i + 1], (double) i);
   }
   U[N] = pow (U[N], (double) N);
   tables_QuickSortD (U, 1, N);
}

/*-------------------------------------------------------------------------*/

void gofs_MergeClasses (double NbExp[], long Loc[],
                        long *smin, long *smax, long *NbClasses)
{
   long s0, j, s;
   double somme;

   *NbClasses = 0;
   s = *smin;
   while (s <= *smax) {
      /* Merge classes to ensure that the number expected in each class is
         >= gofs_MinExpected. */
      if (NbExp[s] < gofs_MinExpected) {
         s0 = s;
         somme = NbExp[s];
         while (somme < gofs_MinExpected && s < *smax) {
            NbExp[s] = 0.0;
            ++s;
            somme += NbExp[s];
         }
         NbExp[s] = somme;
         for (j = s0; j <= s; j++)
            Loc[j] = s;
      } else {
         Loc[s] = s;
      }
      ++*NbClasses;
      ++s;
   }
   *smin = Loc[*smin];

   /* Special case: the last class, if NbExp < MinExpected */
   if (NbExp[*smax] < gofs_MinExpected) {
      if (s0 > *smin)
         --s0;
      NbExp[s0] += NbExp[*smax];
      NbExp[*smax] = 0.0;
      --*NbClasses;
      for (j = s0 + 1; j <= *smax; j++)
         Loc[j] = s0;
      *smax = s0;
   }
   util_Warning (*NbClasses < 2, "gofs_MergeClasses:   NumClasses < 2.\n"
                                 "   The chi-square test is not done.");
   /*
   util_Assert (*NbClasses > 1, "gofs_MergeClasses:   NumClasses < 2");
   */
}


/*-------------------------------------------------------------------------*/

void gofs_WriteClasses (double NbExp[], long Loc[], 
                        long smin, long smax, long NbClasses)
{
   /* Writes the groupings of cells before or after a merging that has */
   /* been done by a previous call to gofs_MergeClasses.  */
   long s, s0;
   double somme;
   const double epsilon = 5.0E-16;

   /* Before merging classes or cells */
   if (NbClasses <= 0) {
      somme = 0.0;
      printf ("-----------------------------------------------\n"
              "Expected numbers per class before merging:\n\n"
              "Class s        NumExpected[s]\n");

      /* Don't print classes for which the expected number < epsilon */
      /* Instead reset smin */
      s = smin;
      while (NbExp[s] < epsilon)
         s++;
      if (s > smin) {
         smin = s;
         s--;
         printf ("<= %3ld", s);
         num_WriteD (NbExp[s], 18, 4, 4);
         printf ("\n");
      }
      /* Reset smax also */
      s0 = s = smax;
      while (NbExp[s] < epsilon)
         s--;
      if (s < smax)
         smax = s;

      /* Now print the classes with their expected numbers */
      for (s = smin; s <= smax; s++) {
         somme += NbExp[s];
         printf ("%6ld", s);
         num_WriteD (NbExp[s], 20, 4, 4);
         printf ("\n");
      }

      if (s0 > smax) {
         s = smax + 1;
         printf (">= %3ld", s);
         num_WriteD (NbExp[s], 18, 4, 4);
         printf ("\n");
      }

      printf ("\n");
      printf ("Total No. Expected = %18.2f\n\n", somme);
      return;
   }

   /* NbClasses > 0: After merging classes */
   printf ("-----------------------------------------------\n"
           "Expected numbers per class after merging:\n"
           "Number of classes: %4ld\n\n", NbClasses);
   printf ("Class s     NumExpected[s]\n");

   somme = 0.0;
   for (s = smin; s <= smax; s++) {
      if (Loc[s] == s) {
         somme += NbExp[s];
         printf ("%4ld %18.4f\n", s, NbExp[s]);
      }
   }
   printf ("\nTotal NumExpected = %18.2f\n\n", somme);
   printf ("The groupings :\n Class s        Loc[s]\n");
   for (s = smin; s <= smax; s++) {
      if (s == smin)
         printf ("<= ");
      else if (s == smax)
         printf (">= ");
      else
         printf ("   ");
      printf ("%4ld  %12ld\n", s, Loc[s]);
   }
   printf ("\n\n");
}


/*-------------------------------------------------------------------------*/

/*******************************\

  Computing EDF test statistics 

\*******************************/


double gofs_Chi2 (double NbExp[], long Count[], long smin, long smax)
{
   double Diff, Khi;
   long s;

   Khi = 0.0;
   for (s = smin; s <= smax; s++) {
      if (NbExp[s] <= 0.0) {
         util_Assert (Count[s] == 0,
            "gofs_Chi2:   NbExp[s] = 0 and Count[s] > 0");
      } else {
         Diff = Count[s] - NbExp[s];
         Khi += Diff * Diff / NbExp[s];
      }
   }
   return Khi;
}

/*-------------------------------------------------------------------------*/

double gofs_Chi2Equal (double NbExp, long Count[], long smin, long smax)
{
   double Diff, Khi;
   long s;
   Khi = 0.0;
   for (s = smin; s <= smax; s++) {
      Diff = Count[s] - NbExp;
      Khi += Diff * Diff;
   }
   return Khi / NbExp;
}

/*-------------------------------------------------------------------------*/

long gofs_Scan (double U[], long N, double d)
{
   long m, j = 1, i = 0;
   double High;

   High = 0.0;
   m = 1;
   while (j < N && High < 1.0) {
      ++i;
      /* Low = U[i]; */
      High = U[i] + d;
      while (j <= N && U[j] < High)
         ++j;
      /* j is now the index of the first obs. to the right of High. */
      if (j - i > m)
         m = j - i;
   }
   /* p-value = fbar_Scan (N, d, m); */
   return m;
}

/*-------------------------------------------------------------------------*/

double gofs_CramerMises (double U[], long N)
{
   long i;
   double W, W2;

   if (N <= 0) {
      util_Warning (TRUE, "gofs_CramerMises:   N <= 0");
      return 0.0;
   }

   W2 = 1.0 / (12 * N);
   for (i = 1; i <= N; i++) {
      W = U[i] - (i - 0.5) / N;
      W2 += W * W;
   }
   return W2;
   /* p-value = fbar_CramerMises (N, W2); */
}

/*-------------------------------------------------------------------------*/

double gofs_WatsonG (double U[], long N)
{
   long i;
   double SumZ;
   double D2;
   double DP, G;
   double UnSurN = 1.0 / N;

   if (N <= 0) {
      util_Warning (TRUE, "gofs_WatsonG:   N <= 0");
      return 0.0;
   }

   /* degenerate case N = 1 */
   if (N == 1)
      return 0.0;

   /* We assume that U is already sorted.  */
   DP = SumZ = 0.0;
   for (i = 1; i <= N; i++) {
      D2 = i * UnSurN - U[i];
      if (D2 > DP)
         DP = D2;
      SumZ += U[i];
   }
   SumZ = SumZ * UnSurN - 0.5;
   G = sqrt ((double) N) * (DP + SumZ);
   return G;
   /* p-value = fbar_WatsonG (N, G); */
}

/*-------------------------------------------------------------------------*/

double gofs_WatsonU (double U[], long N)
{
   long i;
   double SumZ, W, W2, U2;

   if (N <= 0) {
      util_Warning (TRUE, "gofs_WatsonU:   N <= 0");
      return 0.0;
   }

   /* degenerate case N = 1 */
   if (N == 1) {
      return 1.0 / 12.0;
   }

   SumZ = 0.0;
   W2 = 1.0 / (12 * N);
   for (i = 1; i <= N; i++) {
      SumZ += U[i];
      W = U[i] - (i - 0.5) / N;
      W2 += W * W;
   }
   SumZ = SumZ / N - 0.5;
   U2 = W2 - SumZ * SumZ * N;
   return U2;
   /* p-value = fbar_WatsonU (N, U2); */
}

/*-------------------------------------------------------------------------*/

double gofs_AndersonDarling (double V[], long N)
{
   long i;
   double U1;
   double U, A2;

   if (N <= 0) {
      util_Warning (TRUE, "gofs_AndersonDarling:   N <= 0");
      return 0.0;
   }

   A2 = 0.0;
   for (i = 1; i <= N; i++) {
      U1 = U = V[i];
      if (U <= gofs_EpsilonAD) {
         U1 = U = gofs_EpsilonAD;
      } else if (U >= 1 - gofs_EpsilonAD)
         U1 = 1.0 - gofs_EpsilonAD;
      A2 += (2 * i - 1) * log (U) + (1 + 2 * (N - i)) * num2_log1p (-U1);
   }
   A2 = -N - A2 / N;
   return A2;
   /* p-value = fbar_AndersonDarling (N, A2); */
}

/*-------------------------------------------------------------------------*/

void gofs_KSJumpOne (double U[], long N, double a, double *DP, double *DM)
   /* Statistics KS+ and KS-. Case with 1 jump at a, near the lower tail of
      the distribution. */
{
   long j, i;
   double D2, D1, UnSurN;

   if (N <= 0) {
      *DP = *DM = 0.0;
      util_Warning (TRUE, "gofs_KSJumpOne:   N <= 0");
      return;
   }

   *DP = 0.0;
   *DM = 0.0;
   UnSurN = 1.0 / N;
   j = 1;
   while (j < N && U[j] <= a + EpsilonD)
      ++j;
   for (i = j - 1; i <= N; i++) {
      if (i >= 1) {
         D1 = i * UnSurN - U[i];
         if (D1 > *DP)
            *DP = D1;
      }
      if (i >= j) {
         D2 = U[i] - (i - 1) * UnSurN;
         if (D2 > *DM)
            *DM = D2;
      }
   }
}

/*-------------------------------------------------------------------------*/

void gofs_KS (double U[], long N, double *DP, double *DM, double *D)
{
   if (N <= 0) {
      *DP = *DM = *D = 0.0;
      util_Warning (TRUE, "gofs_KS:   N <= 0");
      return;
   }

   gofs_KSJumpOne (U, N, 0.0, DP, DM);
   if (*DM > *DP)
      *D = *DM;
   else
      *D = *DP;
   /*   pp = fbar_KSPlus (N, *DP);
        pm = fbar_KSPlus (N, *DM);
        p  = fbar_KS (N, *D);      */
}

/*-------------------------------------------------------------------------*/
#if 0

void gofs_KSJumpsMany (double X[], int N, wdist_CFUNC F, double W[],
                       double *DP, double *DM, int Detail)
{
   int i;
   double y, UnSurN, D;

   if (N <= 0) {
      *DP = *DM = 0.0;
      util_Warning (TRUE, "gofs_KSJumpsMany:   N <= 0");
      return;
   }

   util_Assert (N > 0, "gofs_KSJumpsMany:   N <= 0");
   UnSurN = 1.0 / N;
   *DP = 0.0;
   *DM = 0.0;

   if (Detail > 0) {
      printf ("-----------------------------------------------\n"
              "Values of the distribution F(x+0) :\n\n");
   }
   /* Assume that the X[i] are already sorted */
   for (i = 1; i <= N; i++) {
      /* Compute KS+ */
      y = F (W, X[i]);
      D = i * UnSurN - y;
      if (D > *DP)
         *DP = D;
      if (Detail > 0) {
         printf ("%14.6f  %14.6f\n", X[i], y);
      }
      /* Compute KS- */
      y = F (W, X[i] - EpsilonD);
      D = y - (i - 1) * UnSurN;
      if (D > *DM)
         *DM = D;
   }
   if (Detail > 0)
      printf ("\n\n");
}
#endif
