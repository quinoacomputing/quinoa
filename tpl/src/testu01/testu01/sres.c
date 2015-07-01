/*************************************************************************\
 *
 * Package:        TestU01
 * File:           sres.c
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
#include "sres.h"
#include "fbar.h"

#include <string.h>
#include <math.h>




/*-------------------------------- Functions ------------------------------*/



void sres_InitChi2 (sres_Chi2 *res, long N, long jmax, char *nam)
{
   statcoll_Init (res->sVal1, N);
   statcoll_Init (res->pVal1, N);

   if (jmax < 0) {
      if (res->jmax > 0) {
         res->NbExp = util_Free (res->NbExp);
         res->Count = util_Free (res->Count);
         res->Loc = util_Free (res->Loc);
      }
   } else {
      if (res->jmax < 0) {
         res->NbExp = util_Calloc ((size_t) (jmax + 1), sizeof (double));
         res->Count = util_Calloc ((size_t) (jmax + 1), sizeof (long));
         res->Loc = util_Calloc ((size_t) (jmax + 1), sizeof (long));
      } else {
         int j;
         res->NbExp =
            util_Realloc (res->NbExp, (jmax + 1) * sizeof (double));
         res->Count = util_Realloc (res->Count, (jmax + 1) * sizeof (long));
         res->Loc = util_Realloc (res->Loc, (jmax + 1) * sizeof (long));
         for (j = 0; j <= jmax; j++) {
            res->NbExp[j] = 0.0;
            res->Count[j] = 0;
            res->Loc[j] = 0;
         }
      }
   }
   res->degFree = 0;
   res->jmin = 0;
   res->jmax = jmax;
   gofw_InitTestArray (res->sVal2, -1.0);
   gofw_InitTestArray (res->pVal2, -1.0);
   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
}


/*-------------------------------------------------------------------------*/

sres_Chi2 *sres_CreateChi2 (void)
{
   sres_Chi2 *res;
   res = util_Malloc (sizeof (sres_Chi2));
   memset (res, 0, sizeof (sres_Chi2));
   res->sVal1 = statcoll_Create (1, "");
   res->pVal1 = statcoll_Create (1, "");
   res->name = util_Calloc (1, sizeof (char));
   res->jmin = 0;
   res->jmax = -1;
   res->NbExp = NULL;
   res->Count = NULL;
   res->Loc = NULL;
   return res;
}


/*-------------------------------------------------------------------------*/

void sres_DeleteChi2 (sres_Chi2 * res)
{
   if (res == NULL)
      return;
   statcoll_Delete (res->sVal1);
   statcoll_Delete (res->pVal1);
   util_Free (res->NbExp);
   util_Free (res->Count);
   util_Free (res->Loc);
   util_Free (res->name);
   util_Free (res);
}


/*-------------------------------------------------------------------------*/

void sres_GetChi2SumStat (sres_Chi2 *res)
{
   const long N = res->sVal1->NObs;
   double sum = N * statcoll_Average (res->sVal1);
   res->sVal2[gofw_Sum] = sum;
   if (N <= 1) {
      res->pVal2[gofw_Sum] = res->sVal1->V[1];
      res->sVal2[gofw_Var] = 0;
      return;
   }
   res->pVal2[gofw_Sum] = fbar_ChiSquare2 (N*res->degFree, 12, sum);
}


/*=========================================================================*/

void sres_InitBasic (sres_Basic *res, long N, char *nam)
{
   statcoll_Init (res->sVal1, N);
   statcoll_Init (res->pVal1, N);
   gofw_InitTestArray (res->sVal2, -1.0);
   gofw_InitTestArray (res->pVal2, -1.0);
   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
}


/*-------------------------------------------------------------------------*/

sres_Basic *sres_CreateBasic (void)
{
   sres_Basic *res;
   res = util_Malloc (sizeof (sres_Basic));
   memset (res, 0, sizeof (sres_Basic));
   res->sVal1 = statcoll_Create (1, "");
   res->pVal1 = statcoll_Create (1, "");
   res->name = util_Calloc (1, sizeof (char));
   return res;
}


/*-------------------------------------------------------------------------*/

void sres_DeleteBasic (sres_Basic * res)
{
   if (res == NULL)
      return;
   statcoll_Delete (res->sVal1);
   statcoll_Delete (res->pVal1);
   util_Free (res->name);
   util_Free (res);
}


/*-------------------------------------------------------------------------*/

void sres_GetNormalSumStat (sres_Basic *res)
{
   const long N = res->sVal1->NObs;
   double sum = N * statcoll_Average (res->sVal1);
   res->sVal2[gofw_Sum] = sum;
   if (N <= 1) {
      res->pVal2[gofw_Sum] = res->sVal1->V[1];
      res->sVal2[gofw_Var] = 0;
      return;
   }
   res->pVal2[gofw_Sum] = fbar_Normal1 (sum/sqrt((double)N));
   sum = statcoll_Variance (res->sVal1);
   res->sVal2[gofw_Var] = sum;
   res->pVal2[gofw_Var] = fbar_ChiSquare2 (N - 1, 12, (N - 1)*sum);
}


/*=========================================================================*/

void sres_InitPoisson (sres_Poisson *res, long N, double Lambda, char *nam)
{
   statcoll_Init (res->sVal1, N);
   res->Lambda = Lambda;
   res->Mu = N * Lambda;
   res->sVal2 = -1.0;
   res->pLeft = -1.0;
   res->pRight = -1.0;
   res->pVal2 = -1.0;
   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
   
}


/*-------------------------------------------------------------------------*/

sres_Poisson * sres_CreatePoisson (void)
{
   sres_Poisson *res;
   res = util_Malloc (sizeof (sres_Poisson));
   memset (res, 0, sizeof (sres_Poisson));
   res->sVal1 = statcoll_Create (1, "");
   res->name = util_Calloc (1, sizeof (char));
   return res;
}


/*-------------------------------------------------------------------------*/

void sres_DeletePoisson (sres_Poisson *res)
{
   if (res == NULL)
      return;
   statcoll_Delete (res->sVal1);
   util_Free (res->name);
   util_Free (res);
}


/*=========================================================================*/

void sres_InitDisc (sres_Disc *res, long N, char *nam)
{
   statcoll_Init (res->sVal1, N);
   res->sVal2 = -1.0;
   res->pLeft = -1.0;
   res->pRight = -1.0;
   res->pVal2 = -1.0;
   res->name = util_Realloc (res->name, 1 + strlen (nam) * sizeof (char));
   strcpy (res->name, nam);
   
}


/*-------------------------------------------------------------------------*/

sres_Disc * sres_CreateDisc (void)
{
   sres_Disc *res;
   res = util_Malloc (sizeof (sres_Disc));
   memset (res, 0, sizeof (sres_Disc));
   res->sVal1 = statcoll_Create (1, "");
   res->name = util_Calloc (1, sizeof (char));
   return res;
}


/*-------------------------------------------------------------------------*/

void sres_DeleteDisc (sres_Disc *res)
{
   if (res == NULL)
      return;
   statcoll_Delete (res->sVal1);
   util_Free (res->name);
   util_Free (res);
}
