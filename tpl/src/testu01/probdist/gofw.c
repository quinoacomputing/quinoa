/*************************************************************************\
 *
 * Package:        ProbDist
 * File:           gofw.c
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


#include "gofw.h"
#include "gofs.h"
#include "fdist.h"
#include "wdist.h"
#include "fbar.h"

#include "util.h"
#include "num.h"
#include "tables.h"
#include "mystr.h"
#include "bitset.h"

#include <float.h>
#include <string.h>
#include <stdio.h>
#include <math.h>




/*---------------------------- extern variables ----------------------------*/

gofw_GraphType gofw_GraphSoft = gofw_Gnuplot;

double gofw_Suspectp = 0.001;

double gofw_Epsilonp = 1.0E-300;
double gofw_Epsilonp1 = 1.0E-15;

char *gofw_TestNames[gofw_NTestTypes] = {
   "KSPlus", "KSMinus", "KS", "Anderson-Darling",
   "Cramer-vonMises", "Watson G", "Watson U",
   "Mean", "Variance", "Correlation", "Sum"
};

bitset_BitSet gofw_ActiveTests = (bitset_BitSet) 0 |
   (1U << gofw_KSP) | (1U << gofw_KSM) | (1U << gofw_AD);



/*---------------------------- module variables ----------------------------*/

#define LEN1 100
#define LEN2 10

static char desc[LEN1];
static char str[LEN2];




/*--------------------------------------------------------------------------*/

static void printMath2 (FILE * f, double x, double y)
{
   /* Writes the pair (x, y) in file f, in a format understood */
   /* by Mathematica */
   char S[41];

   fprintf (f, "   { ");
   if ((x != 0.0) && (x < 0.1 || x > 1.0)) {
      sprintf (S, "%16.7E", x);
      mystr_Subst (S, "E", "*10^(");
      strcat (S, ")");
   } else {
      sprintf (S, "%16.8g", x);
   }
   fprintf (f, S);
   fprintf (f, ",     ");

   if (y != 0.0 && (y < 0.1 || y > 1.0)) {
      sprintf (S, "%16.7E", y);
      mystr_Subst (S, "E", "*10^(");
      strcat (S, ")");
   } else {
      sprintf (S, "%16.8g", y);
   }
   fprintf (f, S);
   fprintf (f, " }");
}

/*--------------------------------------------------------------------------*/

void gofw_GraphDistUnif (FILE * f, double U[], long N, char Desc[])
{
   long i;
   double UnSurN = 1.0 / N;
   if (f == NULL)
      f = stdout;

   switch (gofw_GraphSoft) {

   case gofw_Gnuplot:
      fprintf (f, "#----------------------------------\n");
      fprintf (f, "# %-70s\n\n", Desc);
      fprintf (f, "%16.8g  %16.8g\n", 0.0, 0.0);
      for (i = 1; i <= N; i++)
         fprintf (f, "%16.8g  %16.8g\n", U[i], i * UnSurN);
      fprintf (f, "%16.8g  %16.8g\n\n", 1.0, 1.0);
      break;

   case gofw_Mathematica:
      fprintf (f, "(*----------------------------------*)\n");
      fprintf (f, "(* %-70s\n *)\n\npoints = { \n", Desc);
      printMath2 (f, 0.0, 0.0);
      fprintf (f, ",\n");
      for (i = 1; i <= N; i++) {
         printMath2 (f, U[i], i * UnSurN);
         fprintf (f, ",\n");
      }
      printMath2 (f, 1.0, 1.0);
      fprintf (f, "\n}\n\n");
      break;

   default:
      util_Error ("gofw_GraphDistUnif:   gofw_GraphSoft unknown");
      break;
   }
}


/*--------------------------------------------------------------------------*/

void gofw_GraphFunc (FILE *f, wdist_CFUNC F, double par[], double a,
                     double b, int M, int mono, char Desc[])
{
   int i;
   double yprec, y, x, h;
   if (f == NULL)
      f = stdout;

   switch (gofw_GraphSoft) {

      /* Il y a trop de repetition de code ici.  */
   case gofw_Gnuplot:
      fprintf (f, "#----------------------------------\n");
      fprintf (f, "# %-70s\n\n", Desc);
      h = (b - a) / M;
      if (mono == 1)
         yprec = -DBL_MAX;
      else if (mono == -1)
         yprec = DBL_MAX;
      else
         yprec = 0.0;
      for (i = 0; i <= M; i++) {
         x = a + i * h;
         y = F (par, x);
         fprintf (f, "%16.8g      %16.8g", x, y);
         switch (mono) {
         case 1:
            if (y < yprec)
               fprintf (f, "    #  DECREASING");
            break;
         case -1:
            if (y > yprec)
               fprintf (f, "    #  INCREASING");
            break;
         default:
            break;
         }
         fprintf (f, "\n");
         yprec = y;
      }
      fprintf (f, "\n");
      break;

   case gofw_Mathematica:
      fprintf (f, "(*----------------------------------*)\n");
      fprintf (f, "(* %-70s\n *)\n\npoints = { \n", Desc);
      h = (b - a) / M;
      if (mono == 1)
         yprec = -DBL_MAX;
      else if (mono == -1)
         yprec = DBL_MAX;
      else
         yprec = 0.0;
      for (i = 0; i <= M; i++) {
         x = a + i * h;
         y = F (par, x);
         printMath2 (f, x, y);
         if (i < M)
            fprintf (f, ",");

         switch (mono) {
         case 1:
            if (y < yprec)
               fprintf (f, "   (* DECREASING *)");
            break;
         case -1:
            if (y > yprec)
               fprintf (f, "   (* INCREASING *)");
            break;
         default:
            break;
         }
         fprintf (f, "\n");
         yprec = y;
      }
      fprintf (f, "}\n\n");
      break;

   default:
      util_Error ("gofw_GraphFunc:   gofw_GraphSoft unknown");
      break;
   }
}


/*--------------------------------------------------------------------------*/

double gofw_pDisc (double pLeft, double pRight)
{
   double p;

   if (pRight < pLeft)
      p = pRight;
   else if (pLeft > 0.5)
      p = 0.5;
   else
      p = 1.0 - pLeft;
   /* Note: si p est tres proche de 1, on perd toute la precision ici! */
   /* Note2: je ne pense pas que cela puisse se produire a cause des if (RS) 
    */
   return p;
}


/*--------------------------------------------------------------------------*/

void gofw_Writep0 (double p)
   /* Prints the significance level of a test, without a descriptor */
{
   if ((p >= 0.01) && (p <= 0.99))
      num_WriteD (p, 8, 2, 1);
   else if (p < gofw_Epsilonp)
      printf ("   eps  ");
   else if (p < 0.01)
      num_WriteD (p, 8, 2, 2);
   else if (p >= 1.0 - gofw_Epsilonp1)
      printf (" 1 - eps1");
   else if (p < 1.0 - 1.0e-4)
      printf ("    %.4f", p);
   else {
      printf (" 1 - ");
      num_WriteD (1.0 - p, 7, 2, 2);
   }
}


/*--------------------------------------------------------------------------*/

void gofw_Writep1 (double p)
   /* Prints the significance level of a test, with a descriptor. */
{
/* printf ("Significance level of test            :"); */
   printf ("p-value of test                       :");
   gofw_Writep0 (p);
   if (p < gofw_Suspectp || p > 1.0 - gofw_Suspectp) {
      printf ("    *****");
   }
   printf ("\n\n");
}


/*--------------------------------------------------------------------------*/

void gofw_Writep2 (double x, double p)
   /* Prints the statistic x and its significance level p. */
{
   if ((x < 1.0e5 && x >= 0.1) || (x > -1.0e4 && x <= -0.1))
      num_WriteD (x, 8, 2, 1);
   else if ((x < 0.1 && x >= 0.01) || (x > -0.1 && x <= -0.01))
      num_WriteD (x, 8, 3, 2);
   else
      num_WriteD (x, 8, 3, 3);
   printf ("\n");
   gofw_Writep1 (p);
}


/*--------------------------------------------------------------------------*/

void gofw_WriteKS0 (long N, double DP, double DM, double D)
   /* Prints the results of a Kolmogorov-Smirnov test */
{
   printf ("\n\nKolmogorov-Smirnov+ statistic = D+    :");
   gofw_Writep2 (DP, fbar_KSPlus (N, DP));
   printf ("Kolmogorov-Smirnov- statistic = D-    :");
   gofw_Writep2 (DM, fbar_KSPlus (N, DM));
   printf ("Kolmogorov-Smirnov statistic = D      :");
   gofw_Writep2 (D, fbar_KS1 (N, D));
   printf ("\n\n");
}


/*--------------------------------------------------------------------------*/

void gofw_WriteKS1 (double V[], long N, wdist_CFUNC F, double par[])
{
   double *U;
   double D, DM, DP;

   U = (double *) util_Calloc ((size_t) N + 1, sizeof (double));
   gofs_ContUnifTransform (V, N, F, par, U);
   tables_QuickSortD (U, 1, N);
   gofs_KS (U, N, &DP, &DM, &D);
   gofw_WriteKS0 (N, DP, DM, D);
   util_Free (U);
}


/*--------------------------------------------------------------------------*/

void gofw_WriteKSJumpOne0 (long N, double a, double DP)
{
   double d;

   printf ("\nKolmogorov-Smirnov+ statistic = D+    :%8.2g\n", DP);
   d = 1.0 - fdist_KSPlusJumpOne (N, a, DP);
   gofw_Writep1 (d);
   printf ("\n");
}


/*--------------------------------------------------------------------------*/

void gofw_WriteKSJumpOne1 (double V[], long N, wdist_CFUNC F, double par[],
                           double a)
{
   double *U;
   double DP, DM;

   U = (double *)util_Calloc ((size_t) N + 1, sizeof (double));
   gofs_ContUnifTransform (V, N, F, par, U);
   tables_QuickSortD (U, 1, N);
   gofs_KSJumpOne (U, N, a, &DP, &DM);
   gofw_WriteKSJumpOne0 (N, a, DP);
   util_Free (U);
}


/*--------------------------------------------------------------------------*/

#if 0

void gofw_KSJumpsMany0 (double DP, double DM, fdist_FUNC_JUMPS * H)
{
   double d;

   printf ("\nKolmogorov-Smirnov+ statistic = D+    :%8.2g\n", DP);
   d = 1.0 - fdist_KSPlusJumpsMany (H, DP);
   gofw_Desc1 (d);

   printf ("\nKolmogorov-Smirnov- statistic = D-    :%8.2g\n", DM);
   d = 1.0 - fdist_KSMinusJumpsMany (H, DM);
   gofw_Desc1 (d);
   printf ("\n");
}


/*--------------------------------------------------------------------------*/

void gofw_KSJumpsMany2 (statcoll_Collector *S, fdist_FUNC_JUMPS *H,
                        int Detail)
{
   double DM, DP;
   double *X;
   wdist_CFUNC F = H->F;
   double *W = H->par;

   /* The implementation of fdist_KSPlusJumpsMany and fdist_KSMinusJumpsMany
      works only for NObs <= 64: instability for larger NObs.  */
   if (S->NObs > 64) {
      printf ("\nKolmogorov-Smirnov, sample too large\n\n\n"
	      "------------------------------------------\n");
      return;
   }

   X = (double *) util_Calloc (1 + (size_t) S->NObs, sizeof (double));
   tables_CopyTabD (S->St, X, 1, S->NObs);
   tables_QuickSortD (X, 1, S->NObs);
   statcalc_KSJumpsMany (X, S->NObs, F, W, &DP, &DM, Detail);
   gofw_KSJumpsMany0 (DP, DM, H);
   util_Free (X);
   printf ("\n");
}

#endif
/*--------------------------------------------------------------------------*/

void gofw_InitTestArray (gofw_TestArray A, double x)
{
   int i;
   for (i = 0; i < gofw_NTestTypes; i++)
      A[i] = x;
}


/*--------------------------------------------------------------------------*/

void gofw_Tests0 (double U[], long N, gofw_TestArray sVal)
{
   long i;
   double A2 = 0.0, W2, DM = 0.0, DP = 0.0, W;
   double U1, Ui, D2, D1;
   double SumZ;
   double UnSurN;

   util_Assert (N > 0, "gofw_Tests0:   N <= 0");

   /* We assume that U is already sorted. */
   if (N == 1) {
      sVal[gofw_KSP] = 1.0 - U[1];
      sVal[gofw_Mean] = U[1];
      return;
   }
   UnSurN = 1.0 / N;
   W2 = UnSurN / 12.0;
   SumZ = 0.0;
   for (i = 1; i <= N; i++) {
      /* Statistics KS */
      D1 = U[i] - (i - 1) * UnSurN;
      D2 = i * UnSurN - U[i];
      if (D1 > DM)
         DM = D1;
      if (D2 > DP)
         DP = D2;
      /* Watson U and G */
      SumZ += U[i];
      W = U[i] - (i - 0.5) * UnSurN;
      W2 += W * W;
      /* Anderson-Darling */
      Ui = U[i];
      U1 = 1.0 - Ui;
      if (Ui < gofs_EpsilonAD)
         Ui = gofs_EpsilonAD;
      else if (U1 < gofs_EpsilonAD)
         U1 = gofs_EpsilonAD;
      A2 += (2 * i - 1) * log (Ui) + (1 + 2 * (N - i)) * log (U1);
   }
   if (DM > DP)
      sVal[gofw_KS] = DM;
   else
      sVal[gofw_KS] = DP;
   sVal[gofw_KSM] = DM;
   sVal[gofw_KSP] = DP;
   SumZ = SumZ * UnSurN - 0.5;
   sVal[gofw_CM] = W2;
   sVal[gofw_WG] = sqrt ((double) N) * (DP + SumZ);
   sVal[gofw_WU] = W2 - SumZ * SumZ * N;
   sVal[gofw_AD] = -N - A2 * UnSurN;
/*   sVal[gofw_Mean] = SumZ + 0.5; */ /* Nouveau ... */
}

/*-------------------------------------------------------------------------*/

void gofw_Tests1 (double V[], long N, wdist_CFUNC F, double par[],
                  gofw_TestArray sVal)
{
   double *U;
   util_Assert (N > 0, "gofw_Tests1:   N <= 0");
   U = (double *) util_Calloc ((size_t) N + 1, sizeof (double));
   gofs_ContUnifTransform (V, N, F, par, U);
   tables_QuickSortD (U, 1, N);
   gofw_Tests0 (U, N, sVal);
   if (N == 1)
      sVal[gofw_Mean] = V[1];   /* On veut V[1], pas U[1] */
   util_Free (U);
}

/*-------------------------------------------------------------------------*/

void gofw_ActiveTests0 (double U[], long N, 
                        gofw_TestArray sVal, gofw_TestArray pVal)
{
   util_Assert (N > 0, "gofw_ActiveTests0:   N <= 0");
   if (N == 1) {
      sVal[gofw_Mean] = U[1];
      pVal[gofw_Mean] = 1.0 - U[1];
      sVal[gofw_KSP] = 1.0 - U[1];
      pVal[gofw_KSP] = 1.0 - U[1];
      pVal[gofw_AD] = -1.0;        /* My bug detector */
      return;
   }
   /* We assume that U is already sorted.  */
   gofw_Tests0 (U, N, sVal);

   if (bitset_TestBit (gofw_ActiveTests, gofw_KSP))
      pVal[gofw_KSP] = fbar_KSPlus (N, sVal[gofw_KSP]);

   if (bitset_TestBit (gofw_ActiveTests, gofw_KSM))
      pVal[gofw_KSM] = fbar_KSPlus (N, sVal[gofw_KSM]);

   if (bitset_TestBit (gofw_ActiveTests, gofw_KS))
      pVal[gofw_KS] = fbar_KS1 (N, sVal[gofw_KS]);

   if (bitset_TestBit (gofw_ActiveTests, gofw_AD))
      pVal[gofw_AD] = fbar_AndersonDarling (N, sVal[gofw_AD]);

   if (bitset_TestBit (gofw_ActiveTests, gofw_CM))
      pVal[gofw_CM] = fbar_CramerMises (N, sVal[gofw_CM]);

   if (bitset_TestBit (gofw_ActiveTests, gofw_WG))
      pVal[gofw_WG] = fbar_WatsonG (N, sVal[gofw_WG]);

   if (bitset_TestBit (gofw_ActiveTests, gofw_WU))
      pVal[gofw_WU] = fbar_WatsonU (N, sVal[gofw_WU]);
}

/*-------------------------------------------------------------------------*/

void gofw_ActiveTests1 (double V[], long N, wdist_CFUNC F, double par[],
                        gofw_TestArray sVal, gofw_TestArray pVal)
{
   double *U;
   util_Assert (N > 0, "gofw_ActiveTests1:   N <= 0");
   U = (double *) util_Calloc ((size_t) N + 1, sizeof (double));
   gofs_ContUnifTransform (V, N, F, par, U);
   tables_QuickSortD (U, 1, N);
   gofw_ActiveTests0 (U, N, sVal, pVal);
   if (N == 1)
      sVal[gofw_Mean] = V[1];
   util_Free (U);
}

/*-------------------------------------------------------------------------*/

void gofw_ActiveTests2 (double V[], double U[], long N, wdist_CFUNC F,
   double par[], gofw_TestArray sVal, gofw_TestArray pVal)
{
   util_Assert (N > 0, "gofw_ActiveTests1:   N <= 0");
   tables_QuickSortD (V, 1, N);
   gofs_ContUnifTransform (V, N, F, par, U);
   gofw_ActiveTests0 (U, N, sVal, pVal);
   if (N == 1)
      sVal[gofw_Mean] = V[1];
}

/*-------------------------------------------------------------------------*/

void gofw_WriteActiveTests0 (long N, gofw_TestArray sVal,
                             gofw_TestArray pVal)
{
   if (N == 1) {
      gofw_Writep1 (pVal[gofw_KSP]);
      return;
   }
   printf ("\n");
   if (bitset_TestBit (gofw_ActiveTests, gofw_KSP)) {
      printf ("Kolmogorov-Smirnov+ statistic = D+    :");
      gofw_Writep2 (sVal[gofw_KSP], pVal[gofw_KSP]);
   }
   if (bitset_TestBit (gofw_ActiveTests, gofw_KSM)) {
      printf ("Kolmogorov-Smirnov- statistic = D-    :");
      gofw_Writep2 (sVal[gofw_KSM], pVal[gofw_KSM]);
   }
   if (bitset_TestBit (gofw_ActiveTests, gofw_KS)) {
      printf ("Kolmogorov-Smirnov statistic  = D     :");
      gofw_Writep2 (sVal[gofw_KS], pVal[gofw_KS]);
   }
   if (bitset_TestBit (gofw_ActiveTests, gofw_AD)) {
      printf ("Anderson-Darling statistic = A2       :");
      gofw_Writep2 (sVal[gofw_AD], pVal[gofw_AD]);
   }
   if (bitset_TestBit (gofw_ActiveTests, gofw_CM)) {
      printf ("Cramer-von Mises statistic = W2       :");
      gofw_Writep2 (sVal[gofw_CM], pVal[gofw_CM]);
   }
   if (bitset_TestBit (gofw_ActiveTests, gofw_WG)) {
      printf ("Watson statistic = G                  :");
      gofw_Writep2 (sVal[gofw_WG], pVal[gofw_WG]);
   }
   if (bitset_TestBit (gofw_ActiveTests, gofw_WU)) {
      printf ("Watson statistic = U2                 :");
      gofw_Writep2 (sVal[gofw_WU], pVal[gofw_WU]);
   }
}


/*--------------------------------------------------------------------------*/

void gofw_WriteActiveTests1 (double V[], long N, wdist_CFUNC F, double par[])
{
   gofw_TestArray sv, pv;

   gofw_ActiveTests1 (V, N, F, par, sv, pv);
   gofw_WriteActiveTests0 (N, sv, pv);
}


/*--------------------------------------------------------------------------*/

void gofw_WriteActiveTests2 (long N, gofw_TestArray sVal,
   gofw_TestArray pVal, char S[])
{
   printf ("\n-----------------------------------------------\n");
   if (N == 1) {
      printf (S);
      gofw_Writep2 (sVal[gofw_Mean], pVal[gofw_Mean]);
   } else {
      gofw_WriteActiveTests0 (N, sVal, pVal);
   }
}


/*--------------------------------------------------------------------------*/

void gofw_IterSpacingsTests0 (double U[], long N, int k,
   lebool printval, lebool graph, FILE * f)
   /* Assumes that U is sorted.  */
{
   int j;
   long i;
   double *S, *UU;
   gofw_TestArray sVal, pVal;

   UU = (double *) util_Calloc (1 + (size_t) N, sizeof (double));
   S = (double *) util_Calloc (1 + (size_t) N, sizeof (double));
   printf ("\n");
   for (i = 1; i <= N; i++)
      UU[i] = U[i];               /* UU is a copy of U */
   for (j = 1; j <= k; j++) {
      printf ("-----------------------------------\n"
         "EDF Tests after \"gofw_IterateSpacings\", level :%2d\n", j);
      gofs_DiffD (UU, S, 1, N, 0.0, 1.0);
      gofs_IterateSpacings (UU, S, N);
      tables_QuickSortD (UU, 1, N);
      gofw_ActiveTests0 (UU, N, sVal, pVal);
      gofw_WriteActiveTests0 (N, sVal, pVal);
      strncpy (desc, "Values of Uniforms after IterateSpacings, level ",
         (size_t) LEN1);
      sprintf (str, "%2d", j);
      strncat (desc, str, (size_t) LEN2);
      if (printval > 0)
         tables_WriteTabD (UU, 1, N, 5, 15, 6, 6, desc);
      if (graph > 0)
         gofw_GraphDistUnif (f, UU, N, desc);
   }
   util_Free (UU);
   util_Free (S);
}


/*--------------------------------------------------------------------------*/


void gofw_IterPowRatioTests0 (double U[], long N, int k,
   lebool printval, lebool graph, FILE * f)
{
   int i;
   long j;
   double *UU;
   gofw_TestArray sVal, pVal;

   UU = (double *) util_Calloc (1 + (size_t) N, sizeof (double));
   printf ("\n");
   for (j = 1; j <= N; j++)
      UU[j] = U[j];
   for (i = 1; i <= k; i++) {
      gofs_PowerRatios (UU, N);
      printf ("-----------------------------------\n"
              "EDF Tests after \"gofw_PowerRatios\", level :%2d\n", i);
      tables_QuickSortD (UU, 1, N);
      gofw_ActiveTests0 (UU, N, sVal, pVal);
      gofw_WriteActiveTests0 (N, sVal, pVal);
      strncpy (desc, "Values of Uniforms after PowerRatios, level ",
               (size_t) LEN1);
      sprintf (str, "%2d", i);
      strncat (desc, str, (size_t) LEN2);
      if (printval > 0)
         tables_WriteTabD (UU, 1, N, 5, 15, 6, 6, desc);
      if (graph > 0)
         gofw_GraphDistUnif (f, UU, N, desc);
   }
   util_Free (UU);
}

/*--------------------------------------------------------------------------*/
