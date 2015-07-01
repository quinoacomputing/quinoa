/*************************************************************************\
 *
 * Package:        ProbDist
 * File:           statcoll.c
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

#include "statcoll.h"

#include "util.h"
#include "tables.h"

#include <stddef.h>
#include <string.h>


#define MAXLEN 127                /* Max number of chars in Desc field of a
                                     collector */



/*=========================================================================*/

void statcoll_SetDesc (statcoll_Collector * S, const char name[])
{
   size_t len;
   util_Assert (S != NULL,
      "statcoll_SetDesc: statcoll_Collector is a NULL pointer");
   if (S->Desc != NULL)
      S->Desc = (char *) util_Free (S->Desc);
   if (name == NULL)
      return;
   len = strlen (name);
   if (len > MAXLEN) {
      len = MAXLEN;
      util_Warning (1, "statcoll_Collector->Desc truncated to 127 chars");
   }
   S->Desc = (char *) util_Calloc (len + 1, sizeof (char));
   strncpy (S->Desc, name, (size_t) len);
   S->Desc[len] = '\0';
}


/*=========================================================================*/

statcoll_Collector *statcoll_Create (long N, const char name[])
{
   statcoll_Collector *S;

   util_Warning (N == 0,
      "statcoll_Create:   statcoll_Collector created with N = 0");
   S = (statcoll_Collector *) util_Malloc (sizeof (statcoll_Collector));
   /* We create an array V with N+1 elements, but will keep the N
      observations in elements [1..N]. */
   S->V = (double *) util_Calloc ((size_t) N + 1, sizeof (double));
   S->Dim = N;
   S->NObs = 0;
   S->Desc = NULL;
   statcoll_SetDesc (S, name);
   return S;
}


/*=========================================================================*/

statcoll_Collector *statcoll_Delete (statcoll_Collector * S)
{
   if (S == NULL) {
      util_Warning (S == NULL,
         "statcoll_Delete:   statcoll_Collector is a NULL pointer");
      return NULL;
   }
   S->V = (double *) util_Free (S->V);
   S->Desc = (char *) util_Free (S->Desc);
   util_Free (S);
   return NULL;
}


/*=========================================================================*/

void statcoll_Init (statcoll_Collector * S, long N)
{
   util_Assert (S != NULL,
      "statcoll_Init: statcoll_Collector is a NULL pointer");
   if (N > S->Dim) {
      S->V = (double *) util_Realloc (S->V, (N + 1) * sizeof (double));
      S->Dim = N;
   }
   S->NObs = 0;
}


/*=========================================================================*/

void statcoll_AddObs (statcoll_Collector * S, double x)
{
   util_Assert (S != NULL,
      "statcoll_AddObs:   statcoll_Collector is a NULL pointer");
   if (S->NObs >= S->Dim) {
      if (S->Dim > 0)
         S->Dim *= 2;
      else
         S->Dim = 8;
      S->V = (double *) util_Realloc (S->V, (S->Dim + 1) * sizeof (double));
   }
   ++S->NObs;
   S->V[S->NObs] = x;
}


/*=========================================================================*/

void statcoll_Write (statcoll_Collector * S, int k, int p1, int p2, int p3)
{
   tables_WriteTabD (S->V, 1, S->NObs, k, p1, p2, p3, S->Desc);
}


/*=========================================================================*/

double statcoll_Average (statcoll_Collector * S)
{
   long i;
   double Sum;
   util_Assert (S != NULL,
      "statcoll_Average:   statcoll_Collector is a NULL pointer");
   Sum = 0.0;
   if (S->NObs == 0) {
      util_Warning (1, "statcoll_Average:   NObs = 0");
      return 1.0;
   }
   for (i = 1; i <= S->NObs; i++)
      Sum += S->V[i];
   return Sum / S->NObs;
}


/*=========================================================================*/

double statcoll_Variance (statcoll_Collector * S)
{
   long i;
   double Av, Sum2, Diff;
   util_Assert (S != NULL,
      "statcoll_Variance:   statcoll_Collector is a NULL pointer");
   util_Assert (S->NObs > 1, "statcoll_Variance:   NObs <= 1");
   Av = statcoll_Average (S);
   Sum2 = 0.0;
   for (i = 1; i <= S->NObs; i++) {
      Diff = S->V[i] - Av;
      Sum2 += Diff * Diff;
   }
   return Sum2 / (S->NObs - 1);
}


/*=========================================================================*/

double statcoll_AutoCovar (statcoll_Collector * S, int k)
{
   long i;
   double Av2, Sum2;
   util_Assert (S != NULL,
      "statcoll_AutoCovar:   statcoll_Collector is a NULL pointer");
   util_Assert (k < S->NObs, "statcoll_AutoCovar:   k >= NObs");
   Av2 = statcoll_Average (S);
   Av2 = Av2 * Av2;
   Sum2 = 0.0;
   for (i = 1; i <= S->NObs - k; i++) {
      Sum2 += S->V[i] * S->V[i + k] - Av2;
   }
   return Sum2 / (S->NObs - k);
}


/*=========================================================================*/

double statcoll_Covar (statcoll_Collector * S1, statcoll_Collector * S2)
{
   long i;
   double Av1Av2, Sum;
   util_Assert (S1 != NULL,
      "statcoll_Covar:   statcoll_Collector S1 is a NULL pointer");
   util_Assert (S2 != NULL,
      "statcoll_Covar:   statcoll_Collector S2 is a NULL pointer");
   util_Assert (S1->NObs == S2->NObs,
      "statcoll_Covar:   S1->NObs != S2->NObs");
   util_Assert (S1->NObs > 1, "statcoll_Covar:   NObs <= 1");
   Av1Av2 = statcoll_Average (S1) * statcoll_Average (S2);
   Sum = 0.0;
   for (i = 1; i <= S1->NObs; i++) {
      Sum += S1->V[i] * S2->V[i] - Av1Av2;
   }
   return Sum / (S1->NObs - 1);
}
