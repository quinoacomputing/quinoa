/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fwalk.c
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

#include "fwalk.h"
#include "ffam.h"
#include "fres.h"
#include "fcho.h"
#include "ftab.h"
#include "swalk.h"

#include <string.h>

long fwalk_Maxn = 4194304;
long fwalk_MaxL = 4194304;
double fwalk_MinMu = 1.0 / 4194304;

typedef enum {
   A_VARGEO,
   A_RANDOMWALK1,
   A_RANDOMWALK1A
} TestType;


enum {
   B_ALGOP,
   B_ALGON
};

#define LEN 50



/*------------------------------ Functions --------------------------------*/


static void InitRes1 (ffam_Fam * fam, fwalk_Res1 * res, int N,
   int Nr, int j1, int j2, int jstep, char *name)
/* 
 * Initializes the fwalk_Res1 structure
 */
{
   char str[LEN + 1];
   size_t len;

   strncpy (str, name, (size_t) LEN);
   strncat (str, ", Statistic H", (size_t) LEN - 20);
   len = strlen (str);
   fres_InitCont (fam, res->H, N, Nr, j1, j2, jstep, str);
   str[len - 1] = 'M';
   fres_InitCont (fam, res->M, N, Nr, j1, j2, jstep, str);
   str[len - 1] = 'J';
   fres_InitCont (fam, res->J, N, Nr, j1, j2, jstep, str);
   str[len - 1] = 'R';
   fres_InitCont (fam, res->R, N, Nr, j1, j2, jstep, str);
   str[len - 1] = 'C';
   fres_InitCont (fam, res->C, N, Nr, j1, j2, jstep, str);
}


/*-------------------------------------------------------------------------*/

fwalk_Res1 *fwalk_CreateRes1 (void)
{
   fwalk_Res1 *res;
   res = util_Malloc (sizeof (fwalk_Res1));
   res->H = fres_CreateCont ();
   res->M = fres_CreateCont ();
   res->J = fres_CreateCont ();
   res->R = fres_CreateCont ();
   res->C = fres_CreateCont ();
   return res;
}


/*-------------------------------------------------------------------------*/

void fwalk_DeleteRes1 (fwalk_Res1 * res)
{
   if (res == NULL)
      return;
   fres_DeleteCont (res->H);
   fres_DeleteCont (res->M);
   fres_DeleteCont (res->J);
   fres_DeleteCont (res->R);
   fres_DeleteCont (res->C);
   util_Free (res);
}


/*-------------------------------------------------------------------------*/

static void PrintRes1 (fwalk_Res1 * res)
{
   fres_PrintCont (res->H);
   fres_PrintCont (res->M);
   fres_PrintCont (res->J);
   fres_PrintCont (res->R);
   fres_PrintCont (res->C);
}


/*=========================================================================*/

static void PrintHead (char *name, ffam_Fam * fam, TestType test, void *par1,
   int Nr, int j1, int j2, int jstep)
{
   long *Par = par1;
   double *ParD = par1;

   printf
   ("\n\n================================================================\n");
   printf ("Family:  %s\n\n", fam->name);
   printf ("Test:    %s\n", name);

   if (test == A_VARGEO)
      printf ("   N  = %ld,   n  = %ld,   r = %d",
         (long) ParD[0], (long) ParD[1], (int) ParD[2]);
   else
      printf ("   N  = %ld,   n  = %ld,   r = %d",
         Par[0], Par[1], (int) Par[2]);

   switch (test) {
   case A_VARGEO:
      printf (",   Mu = %f", ParD[3]);
      break;
   case A_RANDOMWALK1:
      printf (",   s = %d,   L  = %ld", (int) Par[3], (long) Par[4]);
      break;
   case A_RANDOMWALK1A:
      printf (",   s = %d,   t = %d,   L  = %ld,   C  = %lu",
         (int) Par[3], (int) Par[4], (long) Par[5], (bitset_BitSet) Par[6]);
      break;
   default:
      util_Error ("in fwalk, PrintHead:  no such case");
   }

   printf ("\n   Nr = %d,   j1 = %d,   j2 = %d,   jstep = %d\n\n",
      Nr, j1, j2, jstep);
}


/*=========================================================================*/

static int ChooseParamRW (int prec, void *cho,
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

   util_Assert (cho, "fwalk:   cho is NULL");
   chon = cho2->Chon;
   choL = cho2->Chop2;
   util_Assert (*pn < 0 || *pL < 0, "fwalk:   Either n or L must be < 0");

   if (*pn < 0) {
      util_Assert (chon, "fwalk:   n < 0 and chon is NULL");
      *pn = fcho_ChooseParamL (chon, (long) (3.0 * gofs_MinExpected),
               fwalk_Maxn, i, j);
      if (*pn <= 0)
         return -1;
   }

   *ps = fcho_Chooses (*pr, *ps, prec);
   if (*ps <= 0)
      return -1;

   if (*pL < 0) {
      util_Assert (choL, "fwalk:   L < 0 and choL is NULL");
      *pL = fcho_ChooseParamL (choL, LMin, fwalk_MaxL, i, j);
      if (*pL < 0)
         return -1;

      /* L must be even */
      if (*pL & 1)
	 ++(*pL);
   }
   return 0;
}


/*=========================================================================*/

static void TabRWalk1 (ffam_Fam * fam, void *res1, void *cho,
   void *par1, int i, int j, int irow, int icol)
{
   int r, s;
   long N, n, L;
   long *Par = par1;
   fwalk_Res1 *fres = res1;
   swalk_Res *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   s = Par[3];
   L = Par[4];

   if (ChooseParamRW (fam->Resol[irow], cho, &n, &r, &s, &L, 8, i, j))
      return;

   sres = swalk_CreateRes ();
   swalk_RandomWalk1 (fam->Gen[irow], sres, N, n, r, s, L, L);
   fres_FillTableEntryC (fres->H, sres->H[0]->pVal2, N, irow, icol);
   fres_FillTableEntryC (fres->M, sres->M[0]->pVal2, N, irow, icol);
   fres_FillTableEntryC (fres->J, sres->J[0]->pVal2, N, irow, icol);
   fres_FillTableEntryC (fres->R, sres->R[0]->pVal2, N, irow, icol);
   fres_FillTableEntryC (fres->C, sres->C[0]->pVal2, N, irow, icol);
   swalk_DeleteRes (sres);
}


/*-------------------------------------------------------------------------*/

void fwalk_RWalk1 (ffam_Fam * fam, fwalk_Res1 * res, fcho_Cho2 * cho,
   long N, long n, int r, int s, long L, int Nr, int j1, int j2, int jstep)
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
      res = fwalk_CreateRes1 ();
   } else
      localRes = FALSE;

   PrintHead ("fwalk_RWalk1", fam, A_RANDOMWALK1, Par, Nr, j1, j2, jstep);

   InitRes1 (fam, res, N, Nr, j1, j2, jstep, "fwalk_RWalk1");
   ftab_MakeTables (fam, res, cho, Par, TabRWalk1, Nr, j1, j2, jstep);
   PrintRes1 (res);
   if (localRes)
      fwalk_DeleteRes1 (res);
}


/*=========================================================================*/

static int ChooseParamVarGeo (void *cho, long *pn, double *pMu, int i, int j)
/*
 * Set the values of the parameters for the test. If a parameter is < 0,
 * will call a choose function to set it. Otherwise, will accept it as is.
 * Returns 0 if parameters are ok for the test, returns -1 if the test
 * should not be done for these parameters.
 */
{
   fcho_Cho2 *cho2 = cho;
   fcho_Cho *chon;
   fcho_Cho *choMu;

   util_Assert (cho, "fwalk:   cho is NULL");
   chon = cho2->Chon;
   choMu = cho2->Chop2;
   util_Assert (*pn < 0 || *pMu < 0, "fwalk:   Either n or Mu must be < 0");

   if (*pn < 0) {
      util_Assert (chon, "fwalk:   n < 0 and chon is NULL");
      *pn = fcho_ChooseParamL (chon, (long) (3.0 * gofs_MinExpected),
               fwalk_Maxn, i, j);
      if (*pn < 0)
         return -1;
   }

   if (*pMu < 0) {
      util_Assert (choMu, "fwalk:   Mu < 0 and choMu is NULL");
      *pMu = choMu->Choose (choMu->param, i, j);

      if (*pMu < fwalk_MinMu) {
         printf ("Mu < %.2g\n\n", fwalk_MinMu);
         return -1;
      }
   }

   return 0;
}


/*=========================================================================*/

static void TabVarGeo (ffam_Fam * fam, void *res1, void *cho,
   void *par1, int i, int j, int irow, int icol)
{
   int r, Algo;
   long N, n;
   double Mu;
   double *Par = par1;
   fres_Cont *fres = res1;
   sres_Chi2 *sres;

   N = Par[0];
   n = Par[1];
   r = Par[2];
   Mu = Par[3];
   Algo = Par[4];

   if (ChooseParamVarGeo (cho, &n, &Mu, i, j))
      return;

   sres = sres_CreateChi2 ();
   if (Algo == B_ALGOP)
      swalk_VarGeoP (fam->Gen[irow], sres, N, n, r, Mu);
   else
      swalk_VarGeoN (fam->Gen[irow], sres, N, n, r, Mu);
   fres_FillTableEntryC (fres, sres->pVal2, N, irow, icol);
   sres_DeleteChi2 (sres);
}


/*------------------------------------------------------------------------*/

static void InVarGeo (ffam_Fam * fam, fres_Cont * res, fcho_Cho2 * cho,
   long N, long n, int r, double Mu, int Algo,
   int Nr, int j1, int j2, int jstep)
{
   double Par[5];
   lebool localRes;
   char Name[30];

   Par[0] = N;
   Par[1] = n;
   Par[2] = r;
   Par[3] = Mu;
   Par[4] = Algo;

   if (res == NULL) {
      localRes = TRUE;
      res = fres_CreateCont ();
   } else
      localRes = FALSE;

   if (Algo == B_ALGOP)
      strcpy (Name, "fwalk_VarGeoP1");
   else
      strcpy (Name, "fwalk_VarGeoN1");

   PrintHead (Name, fam, A_VARGEO, Par, Nr, j1, j2, jstep);
   fres_InitCont (fam, res, N, Nr, j1, j2, jstep, Name);
   ftab_MakeTables (fam, res, cho, Par, TabVarGeo, Nr, j1, j2, jstep);
   fres_PrintCont (res);
   if (localRes)
      fres_DeleteCont (res);
}


/*------------------------------------------------------------------------*/

void fwalk_VarGeoP1 (ffam_Fam * fam, fres_Cont * res, fcho_Cho2 * cho,
   long N, long n, int r, double Mu, int Nr, int j1, int j2, int jstep)
{
   InVarGeo (fam, res, cho, N, n, r, Mu, B_ALGOP, Nr, j1, j2, jstep);
}


/*------------------------------------------------------------------------*/

void fwalk_VarGeoN1 (ffam_Fam * fam, fres_Cont * res, fcho_Cho2 * cho,
   long N, long n, int r, double Mu, int Nr, int j1, int j2, int jstep)
{
   InVarGeo (fam, res, cho, N, n, r, Mu, B_ALGON, Nr, j1, j2, jstep);
}


/*========================================================================*/
