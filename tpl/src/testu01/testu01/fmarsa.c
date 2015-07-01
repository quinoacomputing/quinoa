/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fmarsa.c
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

#include "fmarsa.h"
#include "fcho.h"
#include "ffam.h"
#include "fres.h"
#include "ftab.h"
#include "smarsa.h"
#include "unif01.h"

#include <string.h>
#include <limits.h>
#include <math.h>

long fmarsa_Maxn = 1024 * 1024 * 32;
long fmarsa_MaxL = 1024 * 4;




/*------------------------------ Functions --------------------------------*/


static void InitRes2 (
   ffam_Fam *fam,
   fmarsa_Res2 *res,          /* Results holder */
   int N,                     /* Number of replications */
   int Nr,
   int j1, int j2, int jstep,
   char *name1,
   char *name2
)
/* 
 * Initializes the fmarsa_Res2 structure
 */
{
   fres_InitCont (fam, res->GCD, N, Nr, j1, j2, jstep, name1);
   fres_InitCont (fam, res->NumIter, N, Nr, j1, j2, jstep, name2);
}


/*-------------------------------------------------------------------------*/

fmarsa_Res2 * fmarsa_CreateRes2 (void)
{
   fmarsa_Res2 *res;
   res = util_Malloc (sizeof (fmarsa_Res2));
   res->NumIter = fres_CreateCont ();
   res->GCD = fres_CreateCont ();
   return res;
}


/*-------------------------------------------------------------------------*/

void fmarsa_DeleteRes2 (fmarsa_Res2 *res)
{
   if (res == NULL)
      return;
   fres_DeleteCont (res->GCD);
   fres_DeleteCont (res->NumIter);
   util_Free (res);
}


/*=========================================================================*/

static void PrintHead (char *test, ffam_Fam * fam,
   long N, long n, int r, int s, int L, int t, int p,
   int Nr, int j1, int j2, int jstep)
{
   printf
   ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", test);
   printf ("   N  = %ld,", N);
   if (n)
      printf ("   n = %ld,", n);
   printf ("   r = %d,", r);
   if (s)
      printf ("   s = %d,", s);
   if (L)
      printf ("   L = %d", L);
   if (t)
      printf ("   t = %d,", t);
   if (p)
      printf ("   p = %d", p);
   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);
}


/*=========================================================================*/

static int CheckParamMat (int prec, void *cho,
   long *pn, int *pr, int *ps, long *pL, long LMin, int i, int j)
/*
 * Set the values of the parameters for the test. If a parameter is < 0,
 * will call a choose function to set it. Otherwise, will accept it as is.
 * Returns 0 if parameters are ok for the test, returns -1 if the test
 * should not be done for these parameters.
 */
{
   fcho_Cho2 *cho2 = cho;
   fcho_Cho *chon;
   fcho_Cho *choL;

   util_Assert (cho, "fmarsa:   cho is NULL");
   chon = cho2->Chon;
   choL = cho2->Chop2;
   if (*pn < 0) {
      util_Assert (chon, "fmarsa:   n < 0 and chon is NULL");
      *pn = chon->Choose (chon->param, i, j);

      if (*pn <= 3.0 * gofs_MinExpected) {
         printf ("n is too small\n\n");
         return -1;
      }
      if (*pn > fmarsa_Maxn) {
         printf ("n > %2ld\n\n", fmarsa_Maxn);
         return -1;
      }
   }

   *ps = fcho_Chooses (*pr, *ps, prec);
   if (*ps <= 0)
      return -1;

   if (*pL < 0) {
      util_Assert (choL, "fmarsa:   L < 0 and chop2 is NULL");
      *pL = choL->Choose (choL->param, i, j);

      if (*pL <= LMin) {
         printf ("L is too small\n\n");
         return -1;
      }
      if (*pL > fmarsa_MaxL) {
         printf ("L > %2ld\n\n", fmarsa_MaxL);
         return -1;
      }
   }
   return 0;
}


/*=========================================================================*/


static void TabMatrixR (ffam_Fam * fam, void *res1, void *cho,
   void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n, L;
   const long *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];
   L = Par[4];

   if (CheckParamMat (fam->Resol[irow], cho, &n, &r, &s, &L, 1, i, j))
      return;

   sres = sres_CreateChi2 ();
   smarsa_MatrixRank (fam->Gen[irow], sres, N, n, r, s, L, L);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*-------------------------------------------------------------------------*/

void fmarsa_MatrixR1 (ffam_Fam * fam, fres_Cont * res, fcho_Cho2 * cho,
   long N, long n, int r, int s, int L, int Nr, int j1, int j2, int jstep)
{
   long Par[5] = { 0 };
   lebool localRes;

   Par[0] = N;
   Par[1] = n;
   Par[2] = r;
   Par[3] = s;
   Par[4] = L;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   util_Assert (n < 0 || L < 0, 
      "fmarsa_MatrixR1:   Either n or L must be < 0" );
   PrintHead ("fmarsa_MatrixR1", fam, N, n, r, s, L, 0, 0, Nr, j1, j2,
      jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fmarsa_MatrixR1");
   ftab_MakeTables (fam, res, cho, Par, TabMatrixR, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*=========================================================================*/

static void WriteBirthEC (void *vpar, long junk1, long junk2)
{
   double *Par = vpar;
   double EC = Par[2];
   printf ("Choose d such that EC = %f\n\n", EC);
}

/*-------------------------------------------------------------------------*/

static double ChooseBirthEC (void *vpar, long n, long junk)
{
   double *Par = vpar;
   long N = Par[0];
   int t = Par[1];
   double EC = Par[2];
   long d;
   double k, dr;
   double Mu;

   WriteBirthEC (vpar, 0, 0);
   k = (N * (double) n * n * n) / (4.0 * EC);
   if (k >= smarsa_Maxk) {
      printf ("k >= %2.0f\n\n", smarsa_Maxk);
      return -1.0;
   }
   d = dr = pow (k, 1.0 / t);
   if (dr > LONG_MAX) {
      printf ("d > LONG_MAX\n\n");
      return -1.0;
   }

   k = pow ((double) d, (double) t);
   Mu = N * (double) n * n * n / (4.0 * k);
   if (8.0 * Mu > sqrt (sqrt (k))) {
      printf ("8 EC > k^(1/4)\n\n");
      return -1.0;
   }
   return (double) d;
}

/*-------------------------------------------------------------------------*/

fcho_Cho *fmarsa_CreateBirthEC (long N, int t, double EC)
{
   fcho_Cho *cho;
   double *Par;

   cho = util_Malloc (sizeof (fcho_Cho));
   Par = util_Calloc (3, sizeof (double));
   Par[0] = N;
   Par[1] = t;
   Par[2] = EC;
   cho->param = Par;
   cho->Write = WriteBirthEC;
   cho->Choose = ChooseBirthEC;
   cho->name = util_Calloc (2, sizeof (char));
   strcpy (cho->name, "d");
   return cho;
}

/*-------------------------------------------------------------------------*/

void fmarsa_DeleteBirthEC (fcho_Cho * cho)
{
   if (NULL == cho)
      return;
   cho->name = util_Free (cho->name);
   cho->param = util_Free (cho->param);
   util_Free (cho);
}


/*=========================================================================*/

static int CheckParamBirth (int prec, void *cho,
   long *pn, int *pr, long *pd, int i, int j)
/*
 * Set the values of the parameters for the test.
 * Returns 0 if parameters are ok for the test, returns -1 if the test
 * should not be done for these parameters.
 */
{
   fcho_Cho2 *cho2 = cho;
   fcho_Cho *chon;
   fcho_Cho *chod;
   int s;

   util_Assert (cho, "fmarsa:   cho is NULL");
   chon = cho2->Chon;
   chod = cho2->Chop2;
   util_Assert (chon, "fmarsa:   chon is NULL");
   *pn = chon->Choose (chon->param, i, j);
   if (*pn > fmarsa_Maxn) {
      printf ("n > %2ld\n\n", fmarsa_Maxn);
      return -1;
   }

   util_Assert (chod, "fmarsa:   chop2 is NULL");
   *pd = chod->Choose (chod->param, *pn, 0);
   if (*pd <= 1.0)
      return -1;

   s = num_Log2 ((double) *pd);
   if (*pr + s > prec) {
      printf ("r + Lg(d) > Resolution of generator\n\n");
      return -1;
   }

   return 0;
}


/*=========================================================================*/

static void TabBirthdayS (ffam_Fam * fam, void *vres, void *cho,
   void *vpar, int i, int j, int irow, int icol)
{
   int r, t, p;
   long N, n, d;
   const long *Par = vpar;
   fres_Poisson *fres = vres;
   sres_Poisson *sres;

   N = Par[0];
   r = Par[1];
   t = Par[2];
   p = Par[3];

   if (CheckParamBirth (fam->Resol[irow], cho, &n, &r, &d, i, j))
      return;

   sres = sres_CreatePoisson ();
   smarsa_BirthdaySpacings (fam->Gen[irow], sres, N, n, r, d, t, p);
   fres_FillTableEntryPoisson (fres, sres->Mu, sres->sVal2, sres->pLeft,
      sres->pRight, sres->pVal2, irow, icol);
   sres_DeletePoisson (sres);
}


/*-------------------------------------------------------------------------*/

void fmarsa_BirthdayS1 (ffam_Fam * fam, fres_Poisson * res, fcho_Cho2 * cho,
   long N, int r, int t, int p, int Nr, int j1, int j2, int jstep)
{
   long Par[4] = { 0 };
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = t;
   Par[3] = p;

   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreatePoisson ();
   } else
      localRes = FALSE;

   PrintHead ("fmarsa_BirthdayS1",
      fam, N, 0, r, 0, 0, t, p, Nr, j1, j2, jstep);
   fres_InitPoisson (fam, res, Nr, j1, j2, jstep, "fmarsa_BirthdayS1");
   ftab_MakeTables (fam, res, cho, Par, TabBirthdayS, Nr, j1, j2, jstep);
   ftab_PrintTable2 (res->Exp, res->Obs, FALSE);
   ftab_PrintTable (res->PVal2);
   if (localRes)
      fres_DeletePoisson (res);
}


/*========================================================================*/

void fmarsa_SerialOver1 (void)
{
   util_Error ("fmarsa_SerialOver1:   use fmultin_SerialOver1 instead");
}

/*========================================================================*/

void fmarsa_CollisionOver1 (void)
{
   util_Error ("fmarsa_CollisionOver1:   use fmultin_SerialOver1 instead");
}


/*=========================================================================*/

static void TabGCD (ffam_Fam * fam, void *res1, void *cho,
   void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n;
   const long *Par = par1;
   fmarsa_Res2 *fres = res1;
   smarsa_Res2 *sres;

   N = Par[0];
   r = Par[1];
   s = Par[2];

   n = fcho_ChooseParamL (cho, (long) (3.0 * gofs_MinExpected),
          fmarsa_Maxn, i, j);
   if (n <= 0)
      return;
   s = fcho_Chooses (r, s, fam->Resol[irow]);
   if (s <= 0)
      return;

   sres = smarsa_CreateRes2 ();
   smarsa_GCD (fam->Gen[irow], sres, N, n, r, s);
   fres_FillTableEntryC (fres->GCD, sres->GCD->pVal2, N, irow, icol);
   fres_FillTableEntryC (fres->NumIter, sres->NumIter->pVal2, N, irow, icol);
   smarsa_DeleteRes2 (sres);
}


/*-------------------------------------------------------------------------*/

void fmarsa_GCD1 (ffam_Fam *fam, fmarsa_Res2 *res, fcho_Cho *cho,
   long N, int r, int s, int Nr, int j1, int j2, int jstep)
{
   long Par[3] = { 0 };
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = s;
   if (res == NULL) {
      localRes = TRUE;
      res = fmarsa_CreateRes2 ();
   } else
      localRes = FALSE;

   PrintHead ("fmarsa_GCD1", fam, N, 0, r, s, 0, 0, 0, Nr, j1, j2, jstep);
   InitRes2 (fam, res, N, Nr, j1, j2, jstep,
             "fmarsa_GCD1, GCD", "fmarsa_GCD1, NumIter");
   ftab_MakeTables (fam, res, cho, Par, TabGCD, Nr, j1, j2, jstep);
   fres_PrintCont (res->GCD);
   /*   fres_PrintCont (res->NumIter); */
   if (localRes)
      fmarsa_DeleteRes2 (res);
}


/*=========================================================================*/
