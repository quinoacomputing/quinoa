/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fnpair.c
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
#include "fnpair.h"
#include <string.h>
#include <math.h>

long fnpair_Maxn = 4194304;


typedef enum {
   A_CLOSEPAIRS,
   A_BICKEL,
   A_BITMATCH,
   A_N
} TestType;


/*------------------------------- Functions ------------------------------*/


static void PrintRes1 (fnpair_Res1 * res, TestType test, long N, int m)
{
   switch (test) {
   case A_BITMATCH:
      ftab_PrintTable (res->PVal[snpair_BM]);
      break;
   case A_BICKEL:
      ftab_PrintTable (res->PVal[snpair_BB]);
      break;
   default:
      ftab_PrintTable (res->PVal[snpair_NP]);
      if (m > 1) {
         ftab_PrintTable (res->PVal[snpair_mNP]);
         if (N > 1) {
            ftab_PrintTable (res->PVal[snpair_mNP1]);
            ftab_PrintTable (res->PVal[snpair_mNP2]);
         }
      }
   }
}


/*-------------------------------------------------------------------------*/

static ftab_Table *InitTable (ffam_Fam * fam,
   char *name, int Nr, int j1, int j2, int jstep)
{
   int j;
   ftab_Table *T;

   Nr = util_Min (Nr, fam->Ng);
   T = ftab_CreateTable (Nr, j1, j2, jstep, name, ftab_pVal2, 0);
   ftab_InitMatrix (T, -1.0);
   for (j = 0; j < Nr; j++)
      T->LSize[j] = fam->LSize[j];
   return T;
}


/*-------------------------------------------------------------------------*/

static void InitRes1 (
   ffam_Fam *fam,
   TestType test,
   fnpair_Res1 *res,
   int N,                     /* Number of replications */
   int Nr,
   int j1, int j2, int jstep
)
/* 
 * Initializes the fnpair_Res1 structure
 */
{
   int j;
   Nr = util_Min (Nr, fam->Ng);

   for (j = 0; j < snpair_StatType_N; j++)
      ftab_DeleteTable (res->PVal[j]);
   memset (res, 0, sizeof (fnpair_Res1));

   if (test == A_BITMATCH) {
      res->PVal[snpair_BM] = InitTable (fam, "ClosePairsBitMatch",
         Nr, j1, j2, jstep);
      return;
   }
   if (test == A_BICKEL) {
      res->PVal[snpair_BB] = InitTable (fam,
         "The pVal of AD for BickelBreiman", Nr, j1, j2, jstep);
      return;
   }

   res->PVal[snpair_NP] = InitTable (fam, "", Nr, j1, j2, jstep);
   if (N == 1)
      ftab_SetDesc (res->PVal[snpair_NP],
         "ClosePairs: The closest distance");
   else
      ftab_SetDesc (res->PVal[snpair_NP],
         "ClosePairs: Stat. AD on the N values (NP)");

   res->PVal[snpair_mNP] = InitTable (fam,
      "ClosePairs: A2 test on the values of A2 (m-NP)", Nr, j1, j2, jstep);

   res->PVal[snpair_mNP1] = InitTable (fam,
      "ClosePairs: Test on the Nm values of W_{n,i}(mNP1)",
      Nr, j1, j2, jstep);
   res->PVal[snpair_mNP2] = InitTable (fam,
      "ClosePairs: Stat. AD (mNP2)", Nr, j1, j2, jstep);
}


/*-------------------------------------------------------------------------*/

fnpair_Res1 *fnpair_CreateRes1 (void)
{
   fnpair_Res1 *res;

   res = util_Malloc (sizeof (fnpair_Res1));
   memset (res, 0, sizeof (fnpair_Res1));
   /*   for (j = 0; j < snpair_StatType_N; j++)
	res->PVal[j] = ftab_CreateTable (1, 0, 1, 1, "", ftab_pVal2, 0);*/
   return res;
}


/*-------------------------------------------------------------------------*/

void fnpair_DeleteRes1 (fnpair_Res1 * res)
{
   int j;
   if (res == NULL)
      return;
   for (j = 0; j < snpair_StatType_N; j++)
      ftab_DeleteTable (res->PVal[j]);
   util_Free (res);
}


/*=========================================================================*/

static void PrintHead (char *name, ffam_Fam * fam, TestType test, void *par1,
   int Nr, int j1, int j2, int jstep)
{
   long *Par = par1;

   printf
      ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", name);
   printf ("   N  = %ld,   r = %d,   t = %d", Par[0], (int) Par[1],
      (int) Par[2]);

   switch (test) {
   case A_BITMATCH:
      break;

   case A_CLOSEPAIRS:
      printf (",   p = %d,   m = %d", (int) Par[3], (int) Par[4]);
      break;

   case A_BICKEL:
      printf (",   p = %d,   Torus = ", (int) Par[3]);
      util_WriteBool (Par[4], -5);
      break;

   default:
      util_Error ("in fnpair, PrintHead:  no such case");
   }

   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);
}


/*=========================================================================*/

static void WriteM1 (void *vpar, long junk1, long junk2)
{
   int *Par = vpar;
   int maxm;
   maxm = Par[0];

   if (ftab_Style == ftab_Plain)
      printf ("Choose  m = Min{%d, sqrt(n/sqrt(N)) / 2}\n\n", maxm);
   else
      printf
         ("Choose  $m = \\min\\left\\{%d, \\sqrt{\\frac {n}{4\\sqrt{N}}}\\right\\}$\n\n",
         maxm);
}

/*-------------------------------------------------------------------------*/

static double ChooseM1 (void *vpar, long N, long n)
{
   int m;
   int *Par = vpar;
   int maxm;
   maxm = Par[0];

   WriteM1 (vpar, 0, 0);
   m = sqrt (n / sqrt ((double) N)) / 2.0;
   m = util_Min (m, maxm);
   if (m < 1.0)
      return -1.0;
   else
      return m;
}


/*-------------------------------------------------------------------------*/

fcho_Cho *fnpair_CreateM1 (int maxm)
{
   fcho_Cho *cho;
   int *Par;

   util_Assert (maxm <= snpair_MAXM,
      "fnpair_CreateM1:   maxm > snpair_MAXM");
   cho = util_Malloc (sizeof (fcho_Cho));
   Par = util_Calloc (1, sizeof (int));
   Par[0] = maxm;
   cho->param = Par;
   cho->Write = WriteM1;
   cho->Choose = ChooseM1;
   cho->name = util_Calloc (2, sizeof (char));
   strncpy (cho->name, "m", 1);
   return cho;
}

/*-------------------------------------------------------------------------*/

void fnpair_DeleteM1 (fcho_Cho * cho)
{
   if (NULL == cho)
      return;
   cho->name = util_Free (cho->name);
   cho->param = util_Free (cho->param);
   util_Free (cho);
}


/*=========================================================================*/

static void TabClosePairs (ffam_Fam * fam, void *vres, void *vcho,
   void *par1, int i, int j, int irow, int icol)
{
   int r, t, p, m;
   long N, n;
   long *Par = par1;
   fnpair_Res1 *fres = vres;
   snpair_Res *sres;
   fcho_Cho2 *cho = vcho;
   fcho_Cho *chon;
   fcho_Cho *chom;

   N = Par[0];
   r = Par[1];
   t = Par[2];
   p = Par[3];
   m = Par[4];

   util_Assert (cho, "fnpair:   cho is NULL");
   chon = cho->Chon;
   chom = cho->Chop2;
   n = fcho_ChooseParamL (chon, 2, fnpair_Maxn, i, j);
   if (n <= 0)
      return;
   if (m < 0) {
      util_Assert (chom, "fnpair:   chom is NULL");
      m = chom->Choose (chom->param, N, n);
      if (m < 1)
         return;
   }
   if (n < sqrt ((double) N) * 4 * m * m)
      return;

   sres = snpair_CreateRes ();
   snpair_ClosePairs (fam->Gen[irow], sres, N, n, r, t, p, m);
   fres->PVal[snpair_NP]->Mat[irow][icol] = sres->pVal[snpair_NP];
   if (m > 1)
      fres->PVal[snpair_mNP]->Mat[irow][icol] = sres->pVal[snpair_mNP];
   if ((m > 1) && (N > 1)) {
      fres->PVal[snpair_mNP1]->Mat[irow][icol] = sres->pVal[snpair_mNP1];
      fres->PVal[snpair_mNP2]->Mat[irow][icol] = sres->pVal[snpair_mNP2];
   }
   snpair_DeleteRes (sres);
}


/*------------------------------------------------------------------------*/

void fnpair_ClosePairs1 (ffam_Fam * fam, fnpair_Res1 * res, fcho_Cho2 * cho,
   long N, int r, int t, int p, int m, int Nr, int j1, int j2, int jstep)
{
   long Par[5];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = t;
   Par[3] = p;
   Par[4] = m;
   if (res == NULL) {
      localRes = TRUE;
      res = fnpair_CreateRes1 ();
   } else
      localRes = FALSE;

   PrintHead ("fnpair_ClosePairs1",
      fam, A_CLOSEPAIRS, Par, Nr, j1, j2, jstep);
   InitRes1 (fam, A_CLOSEPAIRS, res, N, Nr, j1, j2, jstep);
   ftab_MakeTables (fam, res, cho, Par, TabClosePairs, Nr, j1, j2, jstep);
   if (m < 0)
      PrintRes1 (res, A_CLOSEPAIRS, N, 2);
   else
      PrintRes1 (res, A_CLOSEPAIRS, N, m);
   if (localRes)
      fnpair_DeleteRes1 (res);
}


/*=========================================================================*/

static void TabBickel (ffam_Fam * fam, void *res1, void *cho,
   void *par1, int i, int j, int irow, int icol)
{
   int r, t, p;
   lebool Torus;
   long N, n;
   long *Par = par1;
   fnpair_Res1 *fres = res1;
   snpair_Res *sres;

   N = Par[0];
   r = Par[1];
   t = Par[2];
   p = Par[3];
   Torus = Par[4];

   n = fcho_ChooseParamL (cho, 2, fnpair_Maxn, i, j);
   if (n <= 0)
      return;

   sres = snpair_CreateRes ();
   snpair_BickelBreiman (fam->Gen[irow], sres, N, n, r, t, p, Torus);
   fres->PVal[snpair_BB]->Mat[irow][icol] = sres->pVal[snpair_BB];
   snpair_DeleteRes (sres);
}


/*------------------------------------------------------------------------*/

void fnpair_Bickel1 (ffam_Fam * fam, fnpair_Res1 * res, fcho_Cho * cho,
   long N, int r, int t, int p, lebool Torus,
   int Nr, int j1, int j2, int jstep)
{
   long Par[5];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = t;
   Par[3] = p;
   Par[4] = Torus;
   if (res == NULL) {
      localRes = TRUE;
      res = fnpair_CreateRes1 ();
   } else
      localRes = FALSE;

   PrintHead ("fnpair_Bickel1", fam, A_BICKEL, Par, Nr, j1, j2, jstep);
   InitRes1 (fam, A_BICKEL, res, N, Nr, j1, j2, jstep);
   ftab_MakeTables (fam, res, cho, Par, TabBickel, Nr, j1, j2, jstep);
   PrintRes1 (res, A_BICKEL, N, 0);
   if (localRes)
      fnpair_DeleteRes1 (res);
}


/*=========================================================================*/

static void TabBitMatch (ffam_Fam * fam, void *res1, void *cho,
   void *par1, int i, int j, int irow, int icol)
{
   int r, t;
   long N, n;
   long *Par = par1;
   fnpair_Res1 *fres = res1;
   snpair_Res *sres;

   N = Par[0];
   r = Par[1];
   t = Par[2];

   n = fcho_ChooseParamL (cho, 2, fnpair_Maxn, i, j);
   if (n <= 0)
      return;

   sres = snpair_CreateRes ();
   snpair_ClosePairsBitMatch (fam->Gen[irow], sres, N, n, r, t);
   fres->PVal[snpair_BM]->Mat[irow][icol] = sres->pVal[snpair_BM];
   snpair_DeleteRes (sres);
}


/*------------------------------------------------------------------------*/

void fnpair_BitMatch1 (ffam_Fam * fam, fnpair_Res1 * res, fcho_Cho * cho,
   long N, int r, int t, int Nr, int j1, int j2, int jstep)
{
   long Par[3];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = t;
   if (res == NULL) {
      localRes = TRUE;
      res = fnpair_CreateRes1 ();
   } else
      localRes = FALSE;

   PrintHead ("fnpair_BitMatch1", fam, A_BITMATCH, Par, Nr, j1, j2, jstep);
   InitRes1 (fam, A_BITMATCH, res, N, Nr, j1, j2, jstep);
   ftab_MakeTables (fam, res, cho, Par, TabBitMatch, Nr, j1, j2, jstep);
   PrintRes1 (res, A_BITMATCH, N, 0);
   if (localRes)
      fnpair_DeleteRes1 (res);
}


/*=========================================================================*/
