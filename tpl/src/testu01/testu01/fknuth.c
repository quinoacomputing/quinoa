/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fknuth.c
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
#include "gofs.h"
#include "gofw.h"

#include "fknuth.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"
#include "ftab.h"
#include "sknuth.h"

#include <string.h>
#include <stdio.h>


long fknuth_Maxn = 4194304;

enum {
   A_GAP,
   A_SIMPPOKER,
   A_COUPONCOLLECTOR,
   A_RUN,
   A_MAXOFT,
   A_N
};


/*------------------------------- Functions ------------------------------*/

static void InitRes1 (
   ffam_Fam *fam,
   fknuth_Res1 *res,         /* Results holder */
   int N,                     /* Number of replications */
   int Nr,
   int j1, int j2, int jstep,
   char *name1,
   char *name2
)
/* 
 * Initializes the fknuth_Res1 structure
 */
{
   fres_InitCont (fam, res->AD, N, Nr, j1, j2, jstep, name2);
   fres_InitCont (fam, res->Chi, N, Nr, j1, j2, jstep, name1);
}


/*-------------------------------------------------------------------------*/

fknuth_Res1 * fknuth_CreateRes1 (void)
{
   fknuth_Res1 *res;
   res = util_Malloc (sizeof (fknuth_Res1));
   res->Chi = fres_CreateCont ();
   res->AD = fres_CreateCont ();
   return res;
}


/*-------------------------------------------------------------------------*/

void fknuth_DeleteRes1 (fknuth_Res1 *res)
{
   if (res == NULL)
      return;
   fres_DeleteCont (res->AD);
   fres_DeleteCont (res->Chi);
   util_Free (res);
}


/*=========================================================================*/

static void PrintHead (char *name, ffam_Fam *fam, int test, void *par1,
   int Nr, int j1, int j2, int jstep)
{
   long *Par = par1;
   double *ParD = par1;

   printf
   ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", name);
   if (test == A_GAP)
      printf ("   N  = %ld,   r = %d", (long) ParD[0], (int) ParD[1]);
   else
      printf ("   N  = %ld,   r = %d", Par[0], (int) Par[1]);

   switch (test) {
   case A_GAP:
      printf (",   Alpha = %f,   Beta = %f", ParD[2], ParD[3]);
      break;
   case A_SIMPPOKER:
      printf (",   d = %d,   k = %d", (int) Par[2], (int) Par[3]);
      break;
   case A_COUPONCOLLECTOR:
      printf (",   d = %d", (int) Par[2]);
      break;
   case A_RUN:
      printf (",   Up = ");   util_WriteBool (Par[2], 5);
      printf (",   Indep = ");   util_WriteBool (Par[3], 5);
      break;
   case A_MAXOFT:
      printf (",   d = %d,   t = %d", (int) Par[2], (int) Par[3]);
      break;
   default:
      util_Error ("in fknuth, PrintHead:  no such case");
   }

   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);
}


/*=========================================================================*/

static void TabGap (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r;
   long N, n;
   double Alpha, Beta;
   double *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   r = Par[1];
   Alpha = Par[2];
   Beta = Par[3];

   n = fcho_ChooseParamL (cho, (long) (gofs_MinExpected /(Beta - Alpha)),
          fknuth_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateChi2 ();
   sknuth_Gap (fam->Gen[irow], sres, N, n, r, Alpha, Beta);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*========================================================================*/

void fknuth_Gap1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                 long N, int r, double Alpha, double Beta,
                 int Nr, int j1, int j2, int jstep)
{
   double Par[4];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = Alpha;
   Par[3] = Beta;

   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fknuth_Gap1", fam, A_GAP, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fknuth_Gap1");
   ftab_MakeTables (fam, res, cho, Par, TabGap, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabSimpPoker (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, d, k;
   long N, n;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   r = Par[1];
   d = Par[2];
   k = Par[3];

   n = fcho_ChooseParamL (cho, (long) (3.0 * gofs_MinExpected),
         fknuth_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateChi2 ();
   sknuth_SimpPoker (fam->Gen[irow], sres, N, n, r, d, k);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*========================================================================*/

void fknuth_SimpPoker1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                       long N, int r, int d, int k,
                       int Nr, int j1, int j2, int jstep)
{
   long Par[4];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = d;
   Par[3] = k;

   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fknuth_SimpPoker1", fam, A_SIMPPOKER, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fknuth_SimpPoker1");
   ftab_MakeTables (fam, res, cho, Par, TabSimpPoker, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabCouponCollector (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, d;
   long N, n;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   r = Par[1];
   d = Par[2];

   n = fcho_ChooseParamL (cho, (long) (3.0 * gofs_MinExpected),
          fknuth_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateChi2 ();
   sknuth_CouponCollector (fam->Gen[irow], sres, N, n, r, d);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*========================================================================*/

void fknuth_CouponCollector1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                             long N, int r, int d,
                             int Nr, int j1, int j2, int jstep)
{
   long Par[3];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = d;
   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead ("fknuth_CouponCollector1",
      fam, A_COUPONCOLLECTOR, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, "fknuth_CouponCollector1");
   ftab_MakeTables (fam, res, cho, Par, TabCouponCollector, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabRun (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r;
   lebool Up, Indep;
   long N, n;
   long *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;
   long nmin;

   N = Par[0];
   r = Par[1];
   Up = Par[2];
   Indep = Par[3];
   if (!Indep)
      nmin = 600;
   else
      nmin = 3.0 * gofs_MinExpected;

   n = fcho_ChooseParamL (cho, nmin, fknuth_Maxn, i, j);
   if (n <= 0)
      return;

   sres = sres_CreateChi2 ();
   if (Indep)
      sknuth_RunIndep (fam->Gen[irow], sres, N, n, r, Up);
   else
      sknuth_Run (fam->Gen[irow], sres, N, n, r, Up);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*========================================================================*/

void fknuth_Run1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho *cho,
                 long N, int r, lebool Up, lebool Indep,
                 int Nr, int j1, int j2, int jstep)
{
   long Par[4];
   lebool localRes;
   char Name[30];

   Par[0] = N;
   Par[1] = r;
   Par[2] = Up;
   Par[3] = Indep;

   if (Indep)
      strcpy (Name, "fknuth_RunIndep1");
   else
      strcpy (Name, "fknuth_Run1");

   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   PrintHead (Name, fam, A_RUN, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, Name);
   ftab_MakeTables (fam, res, cho, Par, TabRun, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*========================================================================*/

static void TabMaxOft (ffam_Fam *fam, void *res1, void *cho,
  void *par1, int i, int j, int irow, int icol)
{
   int r, d, t;
   long N, n;
   long *Par = par1;
   fknuth_Res1 *fres = res1;
   sknuth_Res1 *sres;

   N = Par[0];
   r = Par[1];
   d = Par[2];
   t = Par[3];

   n = fcho_ChooseParamL (cho, (long) (d * gofs_MinExpected), fknuth_Maxn,
           i, j);
   if (n <= 0)
      return;

   sres = sknuth_CreateRes1 ();
   sknuth_MaxOft (fam->Gen[irow], sres, N, n, r, d, t);
   fres_FillTableEntryC (fres->Chi, sres->Chi->pVal2, N, irow, icol);
   fres_FillTableEntryC (fres->AD, sres->Bas->pVal2, N, irow, icol);
   sknuth_DeleteRes1 (sres);
}


/*========================================================================*/

void fknuth_MaxOft1 (ffam_Fam *fam, fknuth_Res1 *res, fcho_Cho *cho,
                    long N, int r, int d, int t,
                    int Nr, int j1, int j2, int jstep)
{
   long Par[4];
   lebool localRes;

   Par[0] = N;
   Par[1] = r;
   Par[2] = d;
   Par[3] = t;
   if (res == NULL) {
      localRes = TRUE;
      res = fknuth_CreateRes1 ();
   } else
      localRes = FALSE;

   PrintHead ("fknuth_MaxOft1", fam, A_MAXOFT, Par, Nr, j1, j2, jstep);
   InitRes1 (fam, res, N, Nr, j1, j2, jstep, "fknuth_MaxOft1, Chi",
      "fknuth_MaxOft1, AD");
   ftab_MakeTables (fam, res, cho, Par, TabMaxOft, Nr, j1, j2, jstep);
   fres_PrintCont (res->Chi);
   fres_PrintCont (res->AD);
   if (localRes)
      fknuth_DeleteRes1 (res);
}


/*========================================================================*/

void fknuth_Serial1 (void)
{
   util_Error ("fknuth_Serial1:   use fmultin_Serial1 instead");
}


/*========================================================================*/

void fknuth_SerialSparse1 (void)
{
   util_Error ("fknuth_SerialSparse1:   use fmultin_Serial1 instead");
}


/*========================================================================*/

void fknuth_Collision1 (void)
{
   util_Error ("fknuth_Collision1:   use fmultin_Serial1 instead");
}


/*========================================================================*/

void fknuth_Permutation1 (void)
{
   util_Error ("fknuth_Permutation1:   use fmultin_Permut1 instead");
}


/*========================================================================*/

void fknuth_CollisionPermut1 (void)
{
   util_Error ("fknuth_CollisionPermut1:   use fmultin_Permut1 instead");
}


/*========================================================================*/
