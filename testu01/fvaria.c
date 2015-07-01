/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fvaria.c
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
#include "gofs.h"
#include "num.h"

#include "fvaria.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"

#include "svaria.h"
#include "sres.h"



long fvaria_MaxN = 1048576 * 4;
long fvaria_Maxn = 1048576 * 4;
long fvaria_Maxk = 1048576 * 4;
long fvaria_MaxK = 1048576 * 4;


enum {
   A_SAMPLEMEAN,
   A_SAMPLECORR,
   A_SAMPLEPROD,
   A_SUMLOGS,
   A_SUMCOLLECTOR,
   A_APPEARANCE,
   A_WEIGHTDISTRIB,
   A_N
};




/*----------------------------- Functions --------------------------------*/


static void PrintHead (char *name, ffam_Fam *fam, int test, void *par1,
   int Nr, int j1, int j2, int jstep)
{
   long *Par = par1;
   double *ParD = par1;

   printf
   ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", name);

   switch (test) {
   case A_SAMPLEMEAN:
      printf ("   n  = %ld,   r = %d", Par[0], (int) Par[1]);
      break;
   case A_SAMPLECORR:
      printf ("   N  = %ld,   r = %d,   k = %d",
              Par[0], (int) Par[1], (int) Par[2]);
      break;
   case A_SAMPLEPROD:
      printf ("   N  = %ld,   r = %d,   t = %d",
              Par[0], (int) Par[1], (int) Par[2]);
      break;
   case A_SUMLOGS:
      printf ("   N  = %ld,   r = %d", Par[0], (int) Par[1]);
      break;
   case A_SUMCOLLECTOR:
      printf ("   N  = %ld,   r = %d,   g = %f", (long) ParD[0], (int) ParD[1],
        ParD[2]);
      break;
   case A_APPEARANCE:
      printf ("   N  = %ld,   r = %d,   s = %d,   L = %d",
              Par[0], (int) Par[1], (int) Par[2], (int) Par[3]);
      break;
   case A_WEIGHTDISTRIB:
      printf ("   N  = %ld,   n  = %ld,   r = %d,   k = %ld,\n   alpha  = %6.4g,   beta = %6.4g",
              (long) ParD[0], (long) ParD[1], (int) ParD[2], (long) ParD[3],
               ParD[4], ParD[5]);
      break;
   default:
      util_Error ("in fknuth, PrintHead:  no such case");
   }

   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);
}


/*=========================================================================*/

static void TabSampleMean (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r;
   long N, n;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Basic *sres;

   n = Par[0];
   r = Par[1];

   N = fcho_ChooseParamL (cho, 0, fvaria_MaxN, i, j);
   if (N <= 0)
      return;

   sres = sres_CreateBasic ();
   svaria_SampleMean (fam->Gen[irow], sres, N, n, r);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteBasic (sres);
}


/*------------------------------------------------------------------------*/

void fvaria_SampleMean1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long n, int r,
                         int Nr, int j1, int j2, int jstep)
{
   long Par[2];
   lebool localRes;

   Par[0] = n;
   Par[1] = r;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fvaria_SampleMean1",
              fam, A_SAMPLEMEAN, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, 2, Nr, j1, j2, jstep, "fvaria_SampleMean1");
   ftab_MakeTables (fam, res, cho, Par, TabSampleMean, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabSampleCorr (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, k;
   long N, n;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Basic *sres;

   N = Par[0];
   r = Par[1];
   k = Par[2];

   n = fcho_ChooseParamL (cho, 2, fvaria_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateBasic ();
   svaria_SampleCorr (fam->Gen[irow], sres, N, n, r, k);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteBasic (sres);
}


/*------------------------------------------------------------------------*/

void fvaria_SampleCorr1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long N, int r, int k,
                         int Nr, int j1, int j2, int jstep)
{
   long Par[3];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = k;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fvaria_SampleCorr1",
              fam, A_SAMPLECORR, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fvaria_SampleCorr1");
   ftab_MakeTables (fam, res, cho, Par, TabSampleCorr, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabSampleProd (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, t;
   long N, n;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Basic *sres;

   N = Par[0];
   r = Par[1];
   t = Par[2];

   n = fcho_ChooseParamL (cho, 1, fvaria_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateBasic ();
   svaria_SampleProd (fam->Gen[irow], sres, N, n, r, t);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteBasic (sres);
}


/*------------------------------------------------------------------------*/

void fvaria_SampleProd1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long N, int r, int t,
                         int Nr, int j1, int j2, int jstep)
{
   long Par[3];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = t;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fvaria_SampleProd1",
              fam, A_SAMPLEPROD, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fvaria_SampleProd1");
   ftab_MakeTables (fam, res, cho, Par, TabSampleProd, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabSumLogs (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r;
   long N, n;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   r = Par[1];

   n = fcho_ChooseParamL (cho, 1, fvaria_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateChi2 ();
   svaria_SumLogs (fam->Gen[irow], sres, N, n, r);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*------------------------------------------------------------------------*/

void fvaria_SumLogs1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                      long N, int r,
                      int Nr, int j1, int j2, int jstep)
{
   long Par[2];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fvaria_SumLogs1",
              fam, A_SUMLOGS, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fvaria_SumLogs1");
   ftab_MakeTables (fam, res, cho, Par, TabSumLogs, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabSumCollector (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r;
   long N, n;
   double g;
   double *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   r = Par[1];
   g = Par[2];

   n = fcho_ChooseParamL (cho, (long) (3.0 * gofs_MinExpected),
        fvaria_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateChi2 ();
   svaria_SumCollector (fam->Gen[irow], sres, N, n, r, g);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*------------------------------------------------------------------------*/

void fvaria_SumCollector1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                           long N, int r, double g,
                           int Nr, int j1, int j2, int jstep)
{
   double Par[3];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = g;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fvaria_SumCollector1",
              fam, A_SUMCOLLECTOR, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fvaria_SumCollector1");
   ftab_MakeTables (fam, res, cho, Par, TabSumCollector, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabAppearance (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, s, L;
   long N, Q, K;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Basic *sres;

   N = Par[0];
   r = Par[1];
   s = Par[2];
   L = Par[3];

   s = fcho_Chooses (r, s, fam->Resol[irow]);
   if (s <= 0)
      return ;
   if ((s > L) && (s % L))
      return ;
   Q = num_TwoExp[L + 4];
   if (Q > fvaria_MaxK) {
      printf ("Q > %ld\n\n", fvaria_MaxK);
      return;
   }
   K = fcho_ChooseParamL (cho, 1, fvaria_MaxK, i, j);
   if (K <= 0)
      return;

   sres = sres_CreateBasic ();
   svaria_AppearanceSpacings (fam->Gen[irow], sres, N, Q, K, r, s, L);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteBasic (sres);
}


/*------------------------------------------------------------------------*/

void fvaria_Appearance1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                         long N, int r, int s, int L,
                         int Nr, int j1, int j2, int jstep)
{
   long Par[4];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = s;
   Par[3] = L;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fvaria_Appearance1",
              fam, A_APPEARANCE, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fvaria_Appearance1");
   ftab_MakeTables (fam, res, cho, Par, TabAppearance, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabWeightDistrib (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r;
   long N, n, k;
   double alpha, beta;
   double *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;
   fcho_Cho2 *cho2 = cho;
   fcho_Cho *chon;
   fcho_Cho *chok;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   k = Par[3];
   alpha = Par[4];
   beta = Par[5];

   util_Assert (cho, "fvaria:   cho is NULL");
   chon = cho2->Chon;
   chok = cho2->Chop2;
   util_Assert (n < 0 || k < 0,
         "fvaria_WeightDistrib1:   Either n or k must be < 0");

   if (n < 0) {
      util_Assert (chon, "fvaria_WeightDistrib1:   n < 0 and chon is NULL");
      n = fcho_ChooseParamL (chon, (long) (3.0 * gofs_MinExpected),
               fvaria_Maxn, i, j);
      if (n <= 0)
         return;
   }

   if (k < 0) {
      util_Assert (chok, "fvaria_WeightDistrib1:   k < 0 and chop2 is NULL");
      k = fcho_ChooseParamL (chok, 1, fvaria_Maxk, i, j);
      if (k <= 0)
         return;
   }

   sres = sres_CreateChi2 ();
   svaria_WeightDistrib (fam->Gen[irow], sres, N, n, r, k, alpha, beta);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*------------------------------------------------------------------------*/

void fvaria_WeightDistrib1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                            long N, long n, int r, long k,
                            double alpha, double beta,
                            int Nr, int j1, int j2, int jstep)
{
   double Par[6];
   lebool localRes;

   Par[0] = N;
   Par[1] = n;
   Par[2] = r;
   Par[3] = k;
   Par[4] = alpha;
   Par[5] = beta;

   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fvaria_WeightDistrib1",
              fam, A_WEIGHTDISTRIB, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fvaria_WeightDistrib1");
   ftab_MakeTables (fam, res, cho, Par, TabWeightDistrib, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/


