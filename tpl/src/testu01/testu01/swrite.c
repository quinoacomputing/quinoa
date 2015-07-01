/*************************************************************************\
 *
 * Package:        TestU01
 * File:           swrite.c
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
#include "chrono.h"

#include "swrite.h"
#include "unif01.h"
#include "gofw.h"

#include <string.h>
#include <stdio.h>
#include <math.h>



#define LEN 100


lebool swrite_Basic = TRUE;
lebool swrite_Parameters = FALSE;
lebool swrite_Collectors = FALSE;
lebool swrite_Counters = FALSE;
lebool swrite_Classes = FALSE;
lebool swrite_Others = FALSE;

lebool swrite_Host = TRUE;

char swrite_ExperimentName[LEN + 1] = "";



/*=========================================================================*/

void swrite_SetExperimentName (char Name[])
{
   strncpy (swrite_ExperimentName, Name, (size_t) LEN);
}


/*=========================================================================*/

void swrite_Head (unif01_Gen *gen, char *TestName, long N, long n, int r)
{
   printf ("***********************************************************\n");
   printf ("HOST = ");
   if (swrite_Host) {
      gdef_WriteHostName ();
      printf ("\n");
   } else
      printf ("\n\n");
   util_Assert (gen != NULL, "No generator has been created");
   unif01_WriteNameGen (gen);
   printf ("\n");
   if (swrite_ExperimentName && strcmp (swrite_ExperimentName, "")) {
      printf ("%s", swrite_ExperimentName);
      printf (":\n\n");
   }
   printf ("%s", TestName);
   printf (":\n-----------------------------------------------\n");
   printf ("   N = %2ld,  n = %2ld,  r = %2d", N, n, r);
   util_Assert (N > 0, "   N <= 0");
   util_Assert (n > 0, "   n <= 0");
   util_Assert (r >= 0, "   r < 0");
}


/*=========================================================================*/

void swrite_Final (unif01_Gen *gen, chrono_Chrono *Timer)
{
   printf ("-----------------------------------------------\n");
   printf ("CPU time used                    :  ");
   chrono_Write (Timer, chrono_hms);
   printf ("\n");
   unif01_WriteState (gen);
   printf ("\n\n\n");
}


/*=========================================================================*/

void swrite_AddStrChi (char S[], int len, long d)
{
   char str[31];
   int j;
   strncpy (S, "Number of degrees of freedom          : ", len);
   j = strlen (S);
   util_Assert (len > j, "swrite_AddStrChi:   len <= j");
   sprintf (str, "%4ld", d);
   strncat (S, str, len - j);
   j = strlen (S);
   util_Assert (len > j, "swrite_AddStrChi *:   len <= j");
   strncat (S, "\nChi-square statistic                  :", len - j);
   S[len - 1] = '\0';
}


/*=========================================================================*/

void swrite_NormalSumTest (long N, sres_Basic *res)
{
   if (N <= 1)
      return;
   printf ("Tests on the sum of all N observations\n");
   printf ("Standardized normal statistic         :");
   gofw_Writep2 (res->sVal2[gofw_Sum]/sqrt((double)N), res->pVal2[gofw_Sum]);
   printf ("Sample variance                       :");
   gofw_Writep2 (res->sVal2[gofw_Var], res->pVal2[gofw_Var]);
}


/*=========================================================================*/
#define LENGTH 200

void swrite_Chi2SumTest (long N, sres_Chi2 *res)
{
   char str[LENGTH + 1];
   if (N <= 1)
      return;
   printf ("Test on the sum of all N observations\n");
   swrite_AddStrChi (str, LENGTH, N*res->degFree);
   printf (str);
   gofw_Writep2 (res->sVal2[gofw_Sum], res->pVal2[gofw_Sum]);
}


/*=========================================================================*/

void swrite_Chi2SumTestb (long N, double sval, double pval, long degFree)
{
   char str[LENGTH + 1];
   if (N <= 1)
      return;
   printf ("Test on the sum of all N observations\n");
   swrite_AddStrChi (str, LENGTH, N*degFree);
   printf (str);
   gofw_Writep2 (sval, pval);
}

