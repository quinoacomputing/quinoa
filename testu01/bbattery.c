/*************************************************************************\
 *
 * Package:        TestU01
 * File:           bbattery.c
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
#include "config.h"
#include "bbattery.h"
#include "smultin.h"
#include "sknuth.h"
#include "smarsa.h"
#include "snpair.h"
#include "svaria.h"
#include "sstring.h"
#include "swalk.h"
#include "scomp.h"
#include "sspectral.h"
#include "swrite.h"
#include "sres.h"
#include "unif01.h"
#include "ufile.h"

#include "gofs.h"
#include "gofw.h"
#include "fdist.h"
#include "fbar.h"
#include "num.h"
#include "chrono.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>



#define LEN 120
#define NAMELEN 30
#define NDIM 200                  /* Dimension of extern arrays */
#define THOUSAND 1000
#define MILLION (THOUSAND * THOUSAND)
#define BILLION (THOUSAND * MILLION)

/* The number of tests in each battery */
#define SMALLCRUSH_NUM 10
#define CRUSH_NUM 96
#define BIGCRUSH_NUM 106
#define RABBIT_NUM 26
#define ALPHABIT_NUM 9


double bbattery_pVal[1 + NDIM] = { 0 };
char *bbattery_TestNames[1 + NDIM] = { 0 };
int bbattery_NTests;

static char CharTemp[LEN + 1];

/* Gives the test number as enumerated in bbattery.tex. Some test applies
   more than one test, so the array of p-values does not correspond with 
   the test number in the doc. */
static int TestNumber[1 + NDIM] = { 0 };




/*-------------------------------- Functions ------------------------------*/


static void GetName (unif01_Gen * gen, char *genName)
{
   char *p;
   int len1, len2;

   if (NULL == gen) {
      genName[0] = '\0';
      return;
   }

   /* Print only the generator name, without the parameters or seeds. */
   /* The parameters start after the first blank; name ends with ':' */
   genName[LEN] = '\0';
   len1 = strcspn (gen->name, ":");
   len1 = util_Min (LEN, len1);
   strncpy (genName, gen->name, (size_t) len1);
   genName[len1] = '\0';
   /* For Filters or Combined generators */
   p = strstr (&gen->name[1 + len1], "unif01");
   while (p != NULL) {
      len1 += 2;
      if (len1 >= LEN)
         return;
      strcat (genName, ", ");
      len2 = strcspn (p, " \0");
      len2 = util_Min (LEN - len1, len2);
      if (len2 <= 0)
         return;
      strncat (genName, p, (size_t) len2);
      len1 = strlen (genName);
      genName[len1] = '\0';
      p += len2;
      p = strstr (p, "unif01");
   }
}


/*=========================================================================*/

static void WritepVal (double p)
/*
 * Write a p-value with a nice format.
 */
{
   if (p < gofw_Suspectp) {
      gofw_Writep0 (p);

   } else if (p > 1.0 - gofw_Suspectp) {
      if (p >= 1.0 - gofw_Epsilonp1) {
         printf (" 1 - eps1");
      } else if (p >= 1.0 - 1.0e-4) {
         printf (" 1 - ");
         num_WriteD (1.0 - p, 7, 2, 2);
         /* printf (" 1 - %.2g ", 1.0 - p); */
      } else if (p >= 1.0 - 1.0e-2)
         printf ("  %.4f ", p);
      else
         printf ("   %.2f", p);
   }
}


/*=========================================================================*/

static void WriteReport (
   char *genName,                 /* Generator or file name */
   char *batName,                 /* Battery name */
   int N,                         /* Max. number of tests */
   double pVal[],                 /* p-values of the tests */
   chrono_Chrono * Timer,         /* Timer */
   lebool Flag,                  /* = TRUE for a file, FALSE for a gen */
   lebool VersionFlag,           /* = TRUE: write the version number */
   double nb                      /* Number of bits in the random file */
   )
{
   int j, co;

   printf ("\n========= Summary results of ");
   printf ("%s", batName);
   printf (" =========\n\n");
   if (VersionFlag)
      printf (" Version:          %s\n", PACKAGE_STRING);
   if (Flag)
      printf (" File:             ");
   else
      printf (" Generator:        ");
   printf ("%s", genName);
   if (nb > 0)
      printf ("\n Number of bits:   %.0f", nb);
   co = 0;
   /* Some of the tests have not been done: their pVal[j] < 0. */
   for (j = 0; j < N; j++) {
      if (pVal[j] >= 0.0)
         co++;
   }
   printf ("\n Number of statistics:  %1d\n", co);
   printf (" Total CPU time:   ");
   chrono_Write (Timer, chrono_hms);

   co = 0;
   for (j = 0; j < N; j++) {
      if (pVal[j] < 0.0)          /* That test was not done: pVal = -1 */
         continue;
      if ((pVal[j] < gofw_Suspectp) || (pVal[j] > 1.0 - gofw_Suspectp)) {
         co++;
         break;
      }
   }
   if (co == 0) {
      printf ("\n\n All tests were passed\n\n\n\n");
      return;
   }

   if (gofw_Suspectp >= 0.01)
      printf ("\n The following tests gave p-values outside [%.4g, %.2f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   else if (gofw_Suspectp >= 0.0001)
      printf ("\n The following tests gave p-values outside [%.4g, %.4f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   else if (gofw_Suspectp >= 0.000001)
      printf ("\n The following tests gave p-values outside [%.4g, %.6f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   else
      printf ("\n The following tests gave p-values outside [%.4g, %.14f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   printf (":\n (eps  means a value < %6.1e)", gofw_Epsilonp);
   printf (":\n (eps1 means a value < %6.1e)", gofw_Epsilonp1);
   printf (":\n\n       Test                          p-value\n");
   printf (" ----------------------------------------------\n");

   co = 0;
   for (j = 0; j < N; j++) {
      if (pVal[j] < 0.0)          /* That test was not done: pVal = -1 */
         continue;
      if ((pVal[j] >= gofw_Suspectp) && (pVal[j] <= 1.0 - gofw_Suspectp))
         continue;                /* That test was passed */
      printf (" %2d ", TestNumber[j]);
      printf (" %-30s", bbattery_TestNames[j]);
      WritepVal (pVal[j]);
      printf ("\n");
      co++;
   }

   printf (" ----------------------------------------------\n");
   if (co < N - 1) {
      printf (" All other tests were passed\n");
   }
   printf ("\n\n\n");
}


/*=========================================================================*/

static void GetPVal_Walk (long N, swalk_Res * res, int *pj, char *mess, int j2)
/*
 * Get the p-values in a swalk_RandomWalk1 test
 */
{
   int j = *pj;
   const unsigned int len = 20;

   if (N == 1) {
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 H");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 M");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 J");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 R");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 C");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

   } else {
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 H");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 M");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 J");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 R");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);

      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (CharTemp, "RandomWalk1 C");
      strncat (CharTemp, mess, (size_t) len);
      strncpy (bbattery_TestNames[j], CharTemp, (size_t) LEN);
   }

   *pj = j;
}


/*=========================================================================*/

static void GetPVal_CPairs (long N, snpair_Res * res, int *pj, char *mess,
   int j2)
/*
 * Get the p-values in a snpair_ClosePairs test
 */
{
   int j = *pj;
   const unsigned int len = 20;

   if (N == 1) {
      bbattery_pVal[++j] = res->pVal[snpair_NP];
      TestNumber[j] = j2;
      strcpy (CharTemp, "ClosePairs NP");
      strncat (CharTemp, mess, (size_t) len);
      strcpy (bbattery_TestNames[j], CharTemp);

      bbattery_pVal[++j] = res->pVal[snpair_mNP];
      TestNumber[j] = j2;
      strcpy (CharTemp, "ClosePairs mNP");
      strncat (CharTemp, mess, (size_t) len);
      strcpy (bbattery_TestNames[j], CharTemp);

   } else {
      bbattery_pVal[++j] = res->pVal[snpair_NP];
      TestNumber[j] = j2;
      strcpy (CharTemp, "ClosePairs NP");
      strncat (CharTemp, mess, (size_t) len);
      strcpy (bbattery_TestNames[j], CharTemp);

      bbattery_pVal[++j] = res->pVal[snpair_mNP];
      TestNumber[j] = j2;
      strcpy (CharTemp, "ClosePairs mNP");
      strncat (CharTemp, mess, (size_t) len);
      strcpy (bbattery_TestNames[j], CharTemp);

      bbattery_pVal[++j] = res->pVal[snpair_mNP1];
      TestNumber[j] = j2;
      strcpy (CharTemp, "ClosePairs mNP1");
      strncat (CharTemp, mess, (size_t) len);
      strcpy (bbattery_TestNames[j], CharTemp);

      bbattery_pVal[++j] = res->pVal[snpair_mNP2];
      TestNumber[j] = j2;
      strcpy (CharTemp, "ClosePairs mNP2");
      strncat (CharTemp, mess, (size_t) len);
      strcpy (bbattery_TestNames[j], CharTemp);

      bbattery_pVal[++j] = res->pVal[snpair_NJumps];
      TestNumber[j] = j2;
      strcpy (CharTemp, "ClosePairs NJumps");
      strncat (CharTemp, mess, (size_t) len);
      strcpy (bbattery_TestNames[j], CharTemp);

      if (snpair_mNP2S_Flag) {
         bbattery_pVal[++j] = res->pVal[snpair_mNP2S];
         TestNumber[j] = j2;
         strcpy (CharTemp, "ClosePairs mNP2S");
         strncat (CharTemp, mess, (size_t) len);
         strcpy (bbattery_TestNames[j], CharTemp);
      }
   }

   *pj = j;
}


/*=========================================================================*/

static void InitBat (void)
/*
 * Initializes the battery of tests: sets all p-values to -1.
 */
{
   int j;
   static int flag = 0;
   for (j = 0; j < NDIM; j++)
      bbattery_pVal[j] = -1.0;
   if (0 == flag) {
      flag++;
      for (j = 0; j < NDIM; j++)
         bbattery_TestNames[j] = util_Calloc (LEN + 1, sizeof (char));
   }
}


/*=========================================================================*/

static void SmallCrush (unif01_Gen * gen, char *filename, int Rep[])
/*
 * A small battery of statistical tests for Random Number Generators 
 * used in simulation.
 * Rep[i] gives the number of times that test i will be done. The default
 * values are Rep[i] = 1 for all i.
 */
{
   const int r = 0;
   int i;
   int j = -1;
   int j2 = 0;
   char genName[LEN + 1] = "";
   chrono_Chrono *Timer;
   sres_Poisson *res1;
   sres_Chi2 *res2;
   sknuth_Res2 *res3;
   swalk_Res *res4;
   sknuth_Res1 *res5;
   sstring_Res *res6;
   lebool fileFlag;

   Timer = chrono_Create ();
   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "                 Starting SmallCrush\n"
         "                 Version: %s\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n",
         PACKAGE_STRING);
   }

   if (NULL == gen) {
      gen = ufile_CreateReadText (filename, 10 * MILLION);
      fileFlag = TRUE;
   } else
      fileFlag = FALSE;

   ++j2;
   if (fileFlag)
      ufile_InitReadText ();
   res1 = sres_CreatePoisson ();
   for (i = 0; i < Rep[j2]; ++i) {
#ifdef USE_LONGLONG
      smarsa_BirthdaySpacings (gen, res1, 1, 5 * MILLION, r, 1073741824,
         2, 1);
#else
      smarsa_BirthdaySpacings (gen, res1, 10, MILLION / 2, r, 67108864, 2, 1);
#endif
      bbattery_pVal[++j] = res1->pVal2;
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "BirthdaySpacings");
   }
   sres_DeletePoisson (res1);

   if (fileFlag)
      ufile_InitReadText ();
   ++j2;
   res3 = sknuth_CreateRes2 ();
   for (i = 0; i < Rep[j2]; ++i) {
      sknuth_Collision (gen, res3, 1, 5 * MILLION, 0, 65536, 2);
      bbattery_pVal[++j] = res3->Pois->pVal2;
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Collision");
   }
   sknuth_DeleteRes2 (res3);

   if (fileFlag)
      ufile_InitReadText ();
   ++j2;
   res2 = sres_CreateChi2 ();
   for (i = 0; i < Rep[j2]; ++i) {
      sknuth_Gap (gen, res2, 1, MILLION / 5, 22, 0.0, .00390625);
      bbattery_pVal[++j] = res2->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Gap");
   }

   ++j2;
   if (fileFlag)
      ufile_InitReadText ();
   for (i = 0; i < Rep[j2]; ++i) {
      sknuth_SimpPoker (gen, res2, 1, 2 * MILLION / 5, 24, 64, 64);
      bbattery_pVal[++j] = res2->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SimpPoker");
   }

   ++j2;
   if (fileFlag)
      ufile_InitReadText ();
   for (i = 0; i < Rep[j2]; ++i) {
      sknuth_CouponCollector (gen, res2, 1, MILLION / 2, 26, 16);
      bbattery_pVal[++j] = res2->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "CouponCollector");
   }

   if (fileFlag)
      ufile_InitReadText ();
   ++j2;
   res5 = sknuth_CreateRes1 ();
   for (i = 0; i < Rep[j2]; ++i) {
      sknuth_MaxOft (gen, res5, 1, 2 * MILLION, 0, MILLION / 10, 6);
      bbattery_pVal[++j] = res5->Chi->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft");
      bbattery_pVal[++j] = res5->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft AD");
   }
   sknuth_DeleteRes1 (res5);

   ++j2;
   if (fileFlag)
      ufile_InitReadText ();
   for (i = 0; i < Rep[j2]; ++i) {
      svaria_WeightDistrib (gen, res2, 1, MILLION / 5, 27, 256, 0.0, 0.125);
      bbattery_pVal[++j] = res2->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "WeightDistrib");
   }

   ++j2;
   if (fileFlag)
      ufile_InitReadText ();
   for (i = 0; i < Rep[j2]; ++i) {
      smarsa_MatrixRank (gen, res2, 1, 20 * THOUSAND, 20, 10, 60, 60);
      bbattery_pVal[++j] = res2->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank");
   }
   sres_DeleteChi2 (res2);

   if (fileFlag)
      ufile_InitReadText ();
   ++j2;
   res6 = sstring_CreateRes ();
   for (i = 0; i < Rep[j2]; ++i) {
      sstring_HammingIndep (gen, res6, 1, MILLION/2, 20, 10, 300, 0);
      bbattery_pVal[++j] = res6->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingIndep");
   }
   sstring_DeleteRes (res6);

   if (fileFlag)
      ufile_InitReadText ();
   ++j2;
   util_Assert (j2 <= SMALLCRUSH_NUM, "SmallCrush:   j2 > SMALLCRUSH_NUM");
   res4 = swalk_CreateRes ();
   for (i = 0; i < Rep[j2]; ++i) {
      swalk_RandomWalk1 (gen, res4, 1, MILLION, r, 30, 150, 150);
      GetPVal_Walk (1, res4, &j, "", j2);
   }
   swalk_DeleteRes (res4);

   bbattery_NTests = ++j;
   if (fileFlag) {
      WriteReport (filename, "SmallCrush", bbattery_NTests, bbattery_pVal,
         Timer, TRUE, TRUE, 0.0);
      ufile_DeleteReadBin (gen);
   } else {
      GetName (gen, genName);
      WriteReport (genName, "SmallCrush", bbattery_NTests, bbattery_pVal,
         Timer, FALSE, TRUE, 0.0);
   }
   chrono_Delete (Timer);
}


/*=========================================================================*/

void bbattery_SmallCrush (unif01_Gen * gen)
{
   int i;
   int Rep[1 + NDIM] = {0};
   for (i = 1; i <= SMALLCRUSH_NUM; ++i)
      Rep[i] = 1;
   SmallCrush (gen, NULL, Rep);
}


/*=========================================================================*/

void bbattery_SmallCrushFile (char *filename)
{
   int i;
   int Rep[1 + NDIM] = {0};
   for (i = 1; i <= SMALLCRUSH_NUM; ++i)
      Rep[i] = 1;
   SmallCrush (NULL, filename, Rep);
}


/*=========================================================================*/

void bbattery_RepeatSmallCrush (unif01_Gen * gen, int Rep[])
{
   SmallCrush (gen, NULL, Rep);
}


/*=========================================================================*/

static void Crush (unif01_Gen * gen, int Rep[])
/*
 * A battery of stringent statistical tests for Random Number Generators
 * used in simulation.
 * Rep[i] gives the number of times that test i will be done. The default
 * values are Rep[i] = 1 for all i.
 */
{
   const int s = 30;
   const int r = 0;
   int i;
   chrono_Chrono *Timer;
   char genName[LEN + 1] = "";
   int j = -1;
   int j2 = 0;

   Timer = chrono_Create ();
   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "                 Starting Crush\n"
         "                 Version: %s\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n",
         PACKAGE_STRING);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_SerialOver (gen, res, 1, 500 * MILLION, 0, 4096, 2);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SerialOver, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_SerialOver (gen, res, 1, 300 * MILLION, 0, 64, 4);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SerialOver, t = 4");
      }
      sres_DeleteBasic (res);
   }
   {
      smarsa_Res *res;
      res = smarsa_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 0, 1024 * 1024, 2);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 10, 1024 * 1024, 2);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 0, 1024, 4);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 20, 1024, 4);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 0, 32, 8);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 25, 32, 8);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 0, 4, 20);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 20");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, res, 10, 10 * MILLION, 28, 4, 20);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 20");
      }
      smarsa_DeleteRes (res);
   }
   {
      sres_Poisson *res;
      res = sres_CreatePoisson ();

#ifdef USE_LONGLONG
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         long d;
#if LONG_MAX <= 2147483647L
         d = 1073741824L;
         smarsa_BirthdaySpacings (gen, res, 10, 10 * MILLION, 0, d, 2, 1);
#else
         d = 2*1073741824L;
         smarsa_BirthdaySpacings (gen, res, 5, 20 * MILLION, 0, d, 2, 1);
#endif
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 5, 20 * MILLION, 0, 2097152, 3,
            1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 3");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 5, 20 * MILLION, 0, 65536, 4, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 3, 20 * MILLION, 0, 512, 7, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 7");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 3, 20 * MILLION, 7, 512, 7, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 7");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 3, 20 * MILLION, 14, 256, 8, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 3, 20 * MILLION, 22, 256, 8, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 8");
      }

#else
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 200, 4 * MILLION / 10, 0,
            67108864, 2, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 100, 4 * MILLION / 10, 0, 131072,
            3, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 3");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 200, 4 * MILLION / 10, 0,
            1024 * 8, 4, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 100, 4 * MILLION / 10, 0, 16, 13,
            1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 100, 4 * MILLION / 10, 10, 16,
            13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 100, 4 * MILLION / 10, 20, 16,
            13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 100, 4 * MILLION / 10, 26, 16,
            13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }
#endif

      sres_DeletePoisson (res);
   }
   {
      lebool flag = snpair_mNP2S_Flag;
      snpair_Res *res;
      res = snpair_CreateRes ();

      snpair_mNP2S_Flag = FALSE;
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairs (gen, res, 10, 2 * MILLION, 0, 2, 0, 30);
         GetPVal_CPairs (10, res, &j, ", t = 2", j2);
      }

      snpair_mNP2S_Flag = TRUE;
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairs (gen, res, 10, 2 * MILLION, 0, 3, 0, 30);
         GetPVal_CPairs (10, res, &j, ", t = 3", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairs (gen, res, 5, 2 * MILLION, 0, 7, 0, 30);
         GetPVal_CPairs (10, res, &j, ", t = 7", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairsBitMatch (gen, res, 4, 4 * MILLION, 0, 2);
         bbattery_pVal[++j] = res->pVal[snpair_BM];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "ClosePairsBitMatch, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairsBitMatch (gen, res, 2, 4 * MILLION, 0, 4);
         bbattery_pVal[++j] = res->pVal[snpair_BM];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "ClosePairsBitMatch, t = 4");
      }
      snpair_DeleteRes (res);
      snpair_mNP2S_Flag = flag;
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 40 * MILLION, 0, 16, 16);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, d = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 40 * MILLION, 26, 16, 16);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, d = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 10 * MILLION, 0, 64, 64);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, d = 64");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 10 * MILLION, 24, 64, 64);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, d = 64");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 40 * MILLION, 0, 4);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, d = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 40 * MILLION, 28, 4);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, d = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 10 * MILLION, 0, 16);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, d = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 10 * MILLION, 26, 16);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, d = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, 100 * MILLION, 0, 0.0, 0.125);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, 100 * MILLION, 27, 0.0, 0.125);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 27");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, 5 * MILLION, 0, 0.0, 1.0/256.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, 5 * MILLION, 22, 0.0, 1.0/256.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 22");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Run (gen, res, 1, 500 * MILLION, 0, TRUE);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of U01, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Run (gen, res, 1, 500 * MILLION, 15, FALSE);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of U01, r = 15");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Permutation (gen, res, 1, 50 * MILLION, 0, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Permutation, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Permutation (gen, res, 1, 50 * MILLION, 15, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Permutation, r = 15");
      }
      sres_DeleteChi2 (res);
   }
   {
      sknuth_Res2 *res;
      res = sknuth_CreateRes2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CollisionPermut (gen, res, 5, 10 * MILLION, 0, 13);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionPermut, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CollisionPermut (gen, res, 5, 10 * MILLION, 15, 13);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionPermut, r = 15");
      }
      sknuth_DeleteRes2 (res);
   }
   {
      sknuth_Res1 *res;
      res = sknuth_CreateRes1 ();

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 10, 10 * MILLION, 0, MILLION / 10, 5);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 5");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 5");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 5, 10 * MILLION, 0, MILLION / 10, 10);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 10");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 10");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 1, 10 * MILLION, 0, MILLION / 10, 20);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 20");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 20");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 1, 10 * MILLION, 0, MILLION / 10, 30);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 30");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 30");
      }
      sknuth_DeleteRes1 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleProd (gen, res, 1, 10 * MILLION, 0, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleProd, t = 10");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleProd (gen, res, 1, 10 * MILLION, 0, 30);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleProd, t = 30");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleMean (gen, res, 10*MILLION, 20, 0);
         bbattery_pVal[++j] = res->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleMean");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleCorr (gen, res, 1, 500 * MILLION, 0, 1);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleCorr");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_AppearanceSpacings (gen, res, 1, 10 * MILLION, 400 * MILLION,
            r, 30, 15);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AppearanceSpacings, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_AppearanceSpacings (gen, res, 1, 10 * MILLION, 100 * MILLION,
            20, 10, 15);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AppearanceSpacings, r = 20");
      }
      sres_DeleteBasic (res);
   }
   {
      smarsa_Res2 *res2;
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 2 * MILLION, 0, 256, 0.0, 0.125);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 2 * MILLION, 8, 256, 0.0, 0.125);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 2 * MILLION, 16, 256, 0.0, 0.125);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 2 * MILLION, 24, 256, 0.0, 0.125);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 24");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SumCollector (gen, res, 1, 20 * MILLION, 0, 10.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SumCollector");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, MILLION, r, s, 2 * s, 2 * s);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, 60 x 60");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, MILLION, 20, 10, 2 * s, 2 * s);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, 60 x 60");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 50 * THOUSAND, r, s, 10 * s, 10 * s);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, 300 x 300");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 50 * THOUSAND, 20, 10, 10 * s,
            10 * s);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, 300 x 300");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 2 * THOUSAND, r, s, 40 * s, 40 * s);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, 1200 x 1200");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 2 * THOUSAND, 20, 10, 40 * s, 40 * s);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, 1200 x 1200");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_Savir2 (gen, res, 1, 20 * MILLION, 0, 1024*1024, 30);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Savir2");
      }
      sres_DeleteChi2 (res);

      res2 = smarsa_CreateRes2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_GCD (gen, res2, 1, 100 * MILLION, 0, 30);
         bbattery_pVal[++j] = res2->GCD->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "GCD, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_GCD (gen, res2, 1, 40 * MILLION, 10, 20);
         bbattery_pVal[++j] = res2->GCD->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "GCD, r = 10");
      }
      smarsa_DeleteRes2 (res2);
   }
   {
      swalk_Res *res;
      res = swalk_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 50 * MILLION, r, s, 90, 90);
         GetPVal_Walk (1, res, &j, " (L = 90)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 10 * MILLION, 20, 10, 90, 90);
         GetPVal_Walk (1, res, &j, " (L = 90)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 5 * MILLION, r, s, 1000, 1000);
         GetPVal_Walk (1, res, &j, " (L = 1000)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, MILLION, 20, 10, 1000, 1000);
         GetPVal_Walk (1, res, &j, " (L = 1000)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, MILLION / 2, r, s, 10000, 10000);
         GetPVal_Walk (1, res, &j, " (L = 10000)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, MILLION / 10, 20, 10, 10000, 10000);
         GetPVal_Walk (1, res, &j, " (L = 10000)", j2);
      }
      swalk_DeleteRes (res);
   }
   {
      scomp_Res *res;
      res = scomp_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LinearComp (gen, res, 1, 120 * THOUSAND, r, 1);
         bbattery_pVal[++j] = res->JumpNum->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 0");
         bbattery_pVal[++j] = res->JumpSize->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LinearComp (gen, res, 1, 120 * THOUSAND, 29, 1);
         bbattery_pVal[++j] = res->JumpNum->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 29");
         bbattery_pVal[++j] = res->JumpSize->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 29");
      }
      scomp_DeleteRes (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LempelZiv (gen, res, 10, 25, r, s);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LempelZiv");
      }
      sres_DeleteBasic (res);
   }
   {
      sspectral_Res *res;
      res = sspectral_CreateRes ();

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sspectral_Fourier3 (gen, res, 50 * THOUSAND, 14, r, s);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Fourier3, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sspectral_Fourier3 (gen, res, 50 * THOUSAND, 14, 20, 10);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Fourier3, r = 20");
      }
      sspectral_DeleteRes (res);
   }
   {
      sstring_Res2 *res;
      res = sstring_CreateRes2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_LongestHeadRun (gen, res, 1, 1000, r, s, 20 + 10 * MILLION);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 0");
         bbattery_pVal[++j] = res->Disc->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_LongestHeadRun (gen, res, 1, 300, 20, 10, 20 + 10 * MILLION);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 20");
         bbattery_pVal[++j] = res->Disc->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 20");
      }
      sstring_DeleteRes2 (res);
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_PeriodsInStrings (gen, res, 1, 300 * MILLION, r, s);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "PeriodsInStrings, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_PeriodsInStrings (gen, res, 1, 300 * MILLION, 15, 15);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "PeriodsInStrings, r = 15");
      }
      sres_DeleteChi2 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingWeight2 (gen, res, 100, 100 * MILLION, r, s, MILLION);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingWeight2, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingWeight2 (gen, res, 30, 100 * MILLION, 20, 10, MILLION);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingWeight2, r = 20");
      }
      sres_DeleteBasic (res);
   }
   {
      sstring_Res *res;
      res = sstring_CreateRes ();
      /* sstring_HammingCorr will probably be removed: less sensitive than
         svaria_HammingIndep */
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingCorr (gen, res, 1, 500 * MILLION, r, s, s);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingCorr, L = 30");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingCorr (gen, res, 1, 50 * MILLION, r, s, 10 * s);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingCorr, L = 300");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingCorr (gen, res, 1, 10 * MILLION, r, s, 40 * s);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingCorr, L = 1200");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 300 * MILLION, r, s, s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L = 30");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 100 * MILLION, 20, 10, s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L = 30");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 30 * MILLION, r, s, 10 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L = 300");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 10 * MILLION, 20, 10, 10 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L = 300");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 10 * MILLION, r, s, 40 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L = 1200");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, MILLION, 20, 10, 40 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L = 1200");
      }
      sstring_DeleteRes (res);
   }
   {
      sstring_Res3 *res;
      res = sstring_CreateRes3 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_Run (gen, res, 1, 1 * BILLION, r, s);
         bbattery_pVal[++j] = res->NRuns->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 0");
         bbattery_pVal[++j] = res->NBits->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_Run (gen, res, 1, 1 * BILLION, 20, 10);
         bbattery_pVal[++j] = res->NRuns->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 20");
         bbattery_pVal[++j] = res->NBits->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 20");
      }
      sstring_DeleteRes3 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 10, 30 + BILLION, r, s, 1);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d = 1");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 5, 1 + BILLION, 20, 10, 1);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d = 1");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 10, 31 + BILLION, r, s, s);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d = 30");
      }

      ++j2;
 /*     util_Assert (j2 <= CRUSH_NUM, "Crush:   j2 > CRUSH_NUM");  */
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 5, 11 + BILLION, 20, 10, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d = 10");
      }
      sres_DeleteBasic (res);
   }

   bbattery_NTests = ++j;
   GetName (gen, genName);
   WriteReport (genName, "Crush", bbattery_NTests,
      bbattery_pVal, Timer, FALSE, TRUE, 0.0);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void bbattery_Crush (unif01_Gen * gen)
{
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= CRUSH_NUM; ++i)
      Rep[i] = 1;
   Crush (gen, Rep);
}


/*=========================================================================*/

void bbattery_RepeatCrush (unif01_Gen * gen, int Rep[])
{
   Crush (gen, Rep);
}


/*=========================================================================*/

static void BigCrush (unif01_Gen * gen, int Rep[])
/*
 * A battery of very stringent statistical tests for Random Number Generators
 * used in simulation.
 * Rep[i] gives the number of times that test i will be done. The default
 * values are Rep[i] = 1 for all i.
 */
{
   const int s = 30;
   const int r = 0;
   int i;
   chrono_Chrono *Timer;
   char genName[LEN + 1] = "";
   int j = -1;
   int j2 = 0;

   Timer = chrono_Create ();
   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "                 Starting BigCrush\n"
         "                 Version: %s\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n",
         PACKAGE_STRING);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_SerialOver (gen, res, 1, BILLION, 0, 256, 3);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SerialOver, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_SerialOver (gen, res, 1, BILLION, 22, 256, 3);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SerialOver, r = 22");
      }
      sres_DeleteBasic (res);
   }
   {
      smarsa_Res *resm;
      resm = smarsa_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 0, 1024*1024*2, 2);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 9, 1024*1024*2, 2);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 0, 1024*16, 3);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 3");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 16, 1024*16, 3);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 3");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 0, 64, 7);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 7");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 24, 64, 7);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 7");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 0, 8, 14);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 14");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 27, 8, 14);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 14");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 0, 4, 21);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 21");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_CollisionOver (gen, resm, 30, 20 * MILLION, 28, 4, 21);
         bbattery_pVal[++j] = resm->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionOver, t = 21");
      }
      smarsa_DeleteRes (resm);
   }
   {
      sres_Poisson *res;
      res = sres_CreatePoisson ();
#ifdef USE_LONGLONG
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         long d;
#if LONG_MAX <= 2147483647L
         d = 1073741824L;
         smarsa_BirthdaySpacings (gen, res, 250, 4 * MILLION, 0, d, 2, 1);
#else
         d = 2147483648L;
         smarsa_BirthdaySpacings (gen, res, 100, 10 * MILLION, 0, d, 2, 1);
#endif
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 20 * MILLION, 0, 2097152, 3,
            1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 3");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 30 * MILLION, 14, 65536, 4, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 20 * MILLION, 0, 512, 7, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 7");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 20 * MILLION, 7, 512, 7, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 7");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 30 * MILLION, 14, 256, 8, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 30 * MILLION, 22, 256, 8, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 30 * MILLION, 0, 16, 16, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 20, 30 * MILLION, 26, 16, 16, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 16");
      }

#else
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 0,
            67108864, 2, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 0,
            1024 * 8, 4, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 16,
            1024 * 8, 4, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 4");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 0, 16,
            13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 5, 16,
            13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 10,
            16, 13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 15,
            16, 13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 20,
            16, 13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_BirthdaySpacings (gen, res, 10 * THOUSAND, MILLION / 10, 26,
            16, 13, 1);
         bbattery_pVal[++j] = res->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings, t = 13");
      }
#endif
      sres_DeletePoisson (res);
   }
   {
      lebool flag = snpair_mNP2S_Flag;
      snpair_Res *res;
      res = snpair_CreateRes ();

      snpair_mNP2S_Flag = TRUE;
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairs (gen, res, 30, 6 * MILLION, 0, 3, 0, 30);
         GetPVal_CPairs (40, res, &j, ", t = 3", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairs (gen, res, 20, 4 * MILLION, 0, 5, 0, 30);
         GetPVal_CPairs (40, res, &j, ", t = 5", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairs (gen, res, 10, 3 * MILLION, 0, 9, 0, 30);
         GetPVal_CPairs (20, res, &j, ", t = 9", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairs (gen, res, 5, 2*MILLION, 0, 16, 0, 30);
         GetPVal_CPairs (10, res, &j, ", t = 16", j2);
      }
      snpair_DeleteRes (res);
      snpair_mNP2S_Flag =flag;
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 400 * MILLION, 0, 8, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 400 * MILLION, 27, 8, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, r = 27");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 100 * MILLION, 0, 32, 32);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_SimpPoker (gen, res, 1, 100 * MILLION, 25, 32, 32);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SimpPoker, r = 25");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 200 * MILLION, 0, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 200 * MILLION, 10, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, r = 10");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 200 * MILLION, 20, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, r = 20");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CouponCollector (gen, res, 1, 200 * MILLION, 27, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CouponCollector, r = 27");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, BILLION/2, 0, 0.0, 1.0/16.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, 300*MILLION, 25, 0.0, 1.0/32.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 25");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, BILLION/10, 0, 0.0, 1.0/128.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Gap (gen, res, 1, 10*MILLION, 20, 0.0, 1.0/1024.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Gap, r = 20");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Run (gen, res, 5, BILLION, 0, FALSE);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Run (gen, res, 10, BILLION, 15, TRUE);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run, r = 15");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Permutation (gen, res, 1, BILLION, 5, 3);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Permutation, t = 3" );
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Permutation (gen, res, 1, BILLION, 5, 5);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Permutation, t = 5");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Permutation (gen, res, 1, BILLION/2, 5, 7);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Permutation, t = 7");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_Permutation (gen, res, 1, BILLION/2, 10, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Permutation, t = 10");
      }
      sres_DeleteChi2 (res);
   }
   {
      sknuth_Res2 *res;
      res = sknuth_CreateRes2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CollisionPermut (gen, res, 20, 20 * MILLION, 0, 14);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionPermut, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_CollisionPermut (gen, res, 20, 20 * MILLION, 10, 14);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "CollisionPermut, r = 10");
      }
      sknuth_DeleteRes2 (res);
   }
   {
      sknuth_Res1 *res;
      res = sknuth_CreateRes1 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 40, 10 * MILLION, 0, MILLION / 10, 8);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 8");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 30, 10 * MILLION, 0, MILLION / 10, 16);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 16");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 20, 10 * MILLION, 0, MILLION / 10, 24);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 24");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 24");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sknuth_MaxOft (gen, res, 20, 10 * MILLION, 0, MILLION / 10, 32);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft, t = 32");
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MaxOft AD, t = 32");
      }
      sknuth_DeleteRes1 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleProd (gen, res, 40, 10 * MILLION, 0, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleProd, t = 8");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleProd (gen, res, 20, 10*MILLION, 0, 16);
         bbattery_pVal[++j] = res->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleProd, t = 16");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleProd (gen, res, 20, 10*MILLION, 0, 24);
         bbattery_pVal[++j] = res->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleProd, t = 24");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleMean (gen, res, 20*MILLION, 30, 0);
         bbattery_pVal[++j] = res->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleMean, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleMean (gen, res, 20*MILLION, 30, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleMean, r = 10");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleCorr (gen, res, 1, 2*BILLION, 0, 1);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleCorr, k = 1");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SampleCorr (gen, res, 1, 2*BILLION, 0, 2);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SampleCorr, k = 2");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_AppearanceSpacings (gen, res, 1, 10 * MILLION, BILLION,
            r, 3, 15);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AppearanceSpacings, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_AppearanceSpacings (gen, res, 1, 10 * MILLION, BILLION,
            27, 3, 15);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AppearanceSpacings, r = 27");
      }
      sres_DeleteBasic (res);
   }
   {
      smarsa_Res2 *res2;
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 20 * MILLION, 0, 256, 0.0, 0.25);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 20 * MILLION, 20, 256, 0.0, 0.25);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 20");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 20 * MILLION, 28, 256, 0.0, 0.25);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 28");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 20 * MILLION, 0, 256, 0.0, 0.0625);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 20 * MILLION, 10, 256, 0.0, 0.0625);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 10");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_WeightDistrib (gen, res, 1, 20 * MILLION, 26, 256, 0.0, 0.0625);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "WeightDistrib, r = 26");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         svaria_SumCollector (gen, res, 1, 500 * MILLION, 0, 10.0);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "SumCollector");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 10, MILLION, r, 5, 30, 30);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, L=30, r=0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 10, MILLION, 25, 5, 30, 30);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, L=30, r=26");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 5 * THOUSAND, r, 4, 1000, 1000);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, L=1000, r=0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 5 * THOUSAND, 26, 4, 1000, 1000);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, L=1000, r=26");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 80, 15, 15, 5000, 5000);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, L=5000");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_MatrixRank (gen, res, 1, 80, 0, 30, 5000, 5000);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank, L=5000");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_Savir2 (gen, res, 10, 10 * MILLION, 10, 1024*1024, 30);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Savir2");
      }
      sres_DeleteChi2 (res);

      res2 = smarsa_CreateRes2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smarsa_GCD (gen, res2, 10, 50 * MILLION, 0, 30);
         bbattery_pVal[++j] = res2->GCD->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "GCD");
      }
      smarsa_DeleteRes2 (res2);
   }
   {
      swalk_Res *res;
      res = swalk_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 100 * MILLION, r, 5, 50, 50);
         GetPVal_Walk (1, res, &j, " (L=50, r=0)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 100 * MILLION, 25, 5, 50, 50);
         GetPVal_Walk (1, res, &j, " (L=50, r=25)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 10 * MILLION, r, 10, 1000, 1000);
         GetPVal_Walk (1, res, &j, " (L=1000, r=0)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 10 * MILLION, 20, 10, 1000, 1000);
         GetPVal_Walk (1, res, &j, " (L=1000, r=20)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 1 * MILLION, r, 15, 10000, 10000);
         GetPVal_Walk (1, res, &j, " (L=10000, r=0)", j2);
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         swalk_RandomWalk1 (gen, res, 1, 1 * MILLION, 15, 15, 10000, 10000);
         GetPVal_Walk (1, res, &j, " (L=10000, r=15)", j2);
      }
      swalk_DeleteRes (res);
   }
   {
      scomp_Res *res;
      res = scomp_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LinearComp (gen, res, 1, 400 * THOUSAND + 20, r, 1);
         bbattery_pVal[++j] = res->JumpNum->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 0");
         bbattery_pVal[++j] = res->JumpSize->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LinearComp (gen, res, 1, 400 * THOUSAND + 20, 29, 1);
         bbattery_pVal[++j] = res->JumpNum->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 29");
         bbattery_pVal[++j] = res->JumpSize->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp, r = 0");
      }
      scomp_DeleteRes (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LempelZiv (gen, res, 10, 27, r, s);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LempelZiv, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LempelZiv (gen, res, 10, 27, 15, 15);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LempelZiv, r = 15");
      }
      sres_DeleteBasic (res);
   }
   {
      sspectral_Res *res;
      res = sspectral_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sspectral_Fourier3 (gen, res, 100 * THOUSAND, 14, r, 3);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Fourier3, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sspectral_Fourier3 (gen, res, 100 * THOUSAND, 14, 27, 3);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Fourier3, r = 27");
      }
      sspectral_DeleteRes (res);
   }
   {
      sstring_Res2 *res;
      res = sstring_CreateRes2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_LongestHeadRun (gen, res, 1, 1000, r, 3, 20 + 10 * MILLION);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 0");
         bbattery_pVal[++j] = res->Disc->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_LongestHeadRun (gen, res, 1, 1000, 27, 3, 20 + 10 * MILLION);
         bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 27");
         bbattery_pVal[++j] = res->Disc->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 27");
      }
      sstring_DeleteRes2 (res);
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_PeriodsInStrings (gen, res, 10, BILLION/2, r, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "PeriodsInStrings, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_PeriodsInStrings (gen, res, 10, BILLION/2, 20, 10);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "PeriodsInStrings, r = 20");
      }
      sres_DeleteChi2 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingWeight2 (gen, res, 10, BILLION, r, 3, MILLION);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingWeight2, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingWeight2 (gen, res, 10, BILLION, 27, 3, MILLION);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingWeight2, r = 27");
      }
      sres_DeleteBasic (res);
   }
   {
      sstring_Res *res;
      res = sstring_CreateRes ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingCorr (gen, res, 1, BILLION, 10, 10, s);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingCorr, L = 30");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingCorr (gen, res, 1, 100 * MILLION, 10, 10, 10 * s);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingCorr, L = 300");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingCorr (gen, res, 1, 100 * MILLION, 10, 10, 40 * s);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingCorr, L = 1200");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 10, 30 * MILLION, r, 3, s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L=30, r=0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 10, 30 * MILLION, 27, 3, s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L=30, r=27");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 30 * MILLION, r, 4, 10 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L=300, r=0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 30 * MILLION, 26, 4, 10 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L=300, r=26");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 10 * MILLION, r, 5, 40 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L=1200, r=0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingIndep (gen, res, 1, 10 * MILLION, 25, 5, 40 * s, 0);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingIndep, L=1200, r=25");
      }
      sstring_DeleteRes (res);
   }
   {
      sstring_Res3 *res;
      res = sstring_CreateRes3 ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_Run (gen, res, 1, 2*BILLION, r, 3);
         bbattery_pVal[++j] = res->NRuns->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 0");
         bbattery_pVal[++j] = res->NBits->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_Run (gen, res, 1, 2*BILLION, 27, 3);
         bbattery_pVal[++j] = res->NRuns->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 27");
         bbattery_pVal[++j] = res->NBits->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits, r = 27");
      }
      sstring_DeleteRes3 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 10, 30 + BILLION, r, 3, 1);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d=1, r=0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 10, 30 + BILLION, r, 3, 3);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d=3, r=0");
      }

      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 10, 30 + BILLION, 27, 3, 1);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d=1, r=27");
      }

      ++j2;
      util_Assert (j2 <= BIGCRUSH_NUM, "BigCrush:   j2 > BIGCRUSH_NUM");
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, 10, 30 + BILLION, 27, 3, 3);
         bbattery_pVal[++j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor, d=3, r=27");
      }
      sres_DeleteBasic (res);
   }

   bbattery_NTests = ++j;
   GetName (gen, genName);
   WriteReport (genName, "BigCrush", bbattery_NTests, bbattery_pVal,
      Timer, FALSE, TRUE, 0.0);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void bbattery_BigCrush (unif01_Gen * gen)
{
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= BIGCRUSH_NUM; ++i)
      Rep[i] = 1;
   BigCrush (gen, Rep);
}


/*=========================================================================*/

void bbattery_RepeatBigCrush (unif01_Gen * gen, int Rep[])
{
   BigCrush (gen, Rep);
}


/*=========================================================================*/
#if 0
static void WriteTime (time_t t0, time_t t1)
{
   int y1;
   double y = 0;

   y = difftime (t1, t0);
   /* printf (" Total time: %.2f sec\n\n", y); */
   printf (" Total time: ");
   y1 = y / 3600;
   printf ("%02d:", y1);
   y -= y1 * 3600.0;
   y1 = y / 60;
   printf ("%02d:", y1);
   y -= y1 * 60.0;
   printf ("%.2f\n\n", y);
}
#endif

/*-------------------------------------------------------------------------*/

static void Alphabit (unif01_Gen * gen, char *fname, double nb, int r, int s,
   lebool blocFlag, int w, int Rep[])
{
   chrono_Chrono *Timer;
 /*  time_t t0, t1; */
   int NbDelta = 1;
   double ValDelta[] = { 1 };
   long N = 1;
   long n, L;
   int j = 0;
   int j2 = 0;
   int i;
   lebool fileFlag;
   long bufsiz;
   char genName[LEN + 1] = "";
   double z;
   unif01_Gen *gen0;

   Timer = chrono_Create ();
 /*  t0 = time (NULL); */
   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "          Starting Alphabit:   nb = %.0f\n"
         "          Version: %s\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n",
         nb, PACKAGE_STRING);
   }
   util_Assert (nb > 0, "Alphabit:   nb <= 0");
   /* Bits will be read as 32-bit unsigned integers */
   nb -= fmod (nb, 32.0);
   bufsiz = nb / 32.0;

   if (blocFlag) {
      gen0 = ufile_CreateReadBin (fname, bufsiz);
      gen = unif01_CreateBitBlockGen (gen0, r, s, w);
      nb -= fmod (nb, 1024.0 / w);
      fileFlag = TRUE;
   } else if (NULL == gen) {
      gen = ufile_CreateReadBin (fname, bufsiz);
      fileFlag = TRUE;
   } else {
      fileFlag = FALSE;
   }

   {
      smultin_Param *par = NULL;
      smultin_Res *res;
      par = smultin_CreateParam (NbDelta, ValDelta, smultin_GenerCellSerial,
         3);
      res = smultin_CreateRes (par);
      if (fileFlag)
         ufile_InitReadBin ();

      if (nb > BILLION)
         N = 1 + nb / BILLION;
      else
         N = 1;
      n = nb / N;
      /* Set n as a multiple of s = 32 */
      n -= n % 32;
      j = -1;
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smultin_MultinomialBitsOver (gen, par, res, N, n, r, s, 2, FALSE);
         strcpy (bbattery_TestNames[++j], "MultinomialBitsOver, L = 2");
         if (N == 1)
            bbattery_pVal[j] = res->pVal2[0][gofw_Mean];
         else
            bbattery_pVal[j] = res->pVal2[0][gofw_AD];
         TestNumber[j] = j2;
      }

      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         smultin_MultinomialBitsOver (gen, par, res, N, n, r, s, 4, FALSE);
         strcpy (bbattery_TestNames[++j], "MultinomialBitsOver, L = 4");
         if (N == 1)
            bbattery_pVal[j] = res->pVal2[0][gofw_Mean];
         else
            bbattery_pVal[j] = res->pVal2[0][gofw_AD];
         TestNumber[j] = j2;
      }

      ++j2;
      if (n > 250) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            smultin_MultinomialBitsOver (gen, par, res, N, n, r, s, 8, FALSE);
            strcpy (bbattery_TestNames[++j], "MultinomialBitsOver, L = 8");
            if (N == 1)
               bbattery_pVal[j] = res->pVal2[0][gofw_Mean];
            else
               bbattery_pVal[j] = res->pVal2[0][gofw_AD];
           TestNumber[j] = j2;
         }
      }

      ++j2;
      if (n > 65000) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            smultin_MultinomialBitsOver (gen, par, res, N, n, r, s, 16, FALSE);
            strcpy (bbattery_TestNames[++j], "MultinomialBitsOver, L = 16");
            if (N == 1)
               bbattery_pVal[j] = res->pVal2[0][gofw_Mean];
            else
               bbattery_pVal[j] = res->pVal2[0][gofw_AD];
           TestNumber[j] = j2;
         }
      }

      smultin_DeleteRes (res);
      smultin_DeleteParam (par);
   }

   {
      sstring_Res *res;
      res = sstring_CreateRes ();

      if (fileFlag)
         ufile_InitReadBin ();
      z = nb / s;
      N = 1 + z / BILLION;
      n = z / N;
      ++j2;
      if (n >= 20) {
         for (i = 0; i < Rep[j2]; ++i) {
            sstring_HammingIndep (gen, res, N, n, r, s, 16, 0);
            j++;
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            strcpy (bbattery_TestNames[j], "HammingIndep, L = 16");
            TestNumber[j] = j2;
         }
      }

      if (fileFlag)
         ufile_InitReadBin ();
      n /= 2;
      ++j2;
      if (n >= 20) {
         for (i = 0; i < Rep[j2]; ++i) {
            sstring_HammingIndep (gen, res, N, n, r, s, 32, 0);
            j++;
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            strcpy (bbattery_TestNames[j], "HammingIndep, L = 32");
            TestNumber[j] = j2;
         }
      }

      if (fileFlag)
         ufile_InitReadBin ();
      n *= 2;
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_HammingCorr (gen, res, N, n, r, s, 32);
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "HammingCorr, L = 32");
      }
   }

   {
      swalk_Res *res;
      res = swalk_CreateRes ();

      if (fileFlag)
         ufile_InitReadBin ();
      L = 64;
      z = nb / L;
      N = 1 + z / BILLION;
      n = z / N;
      j2++;
      if (n >= 30) {
         for (i = 0; i < Rep[j2]; ++i) {
            swalk_RandomWalk1 (gen, res, N, n, r, s, L, L);
            GetPVal_Walk (N, res, &j, " (L = 64)", j2);
         }
      }

      if (fileFlag)
         ufile_InitReadBin ();
      L = 320;
      z = nb / L;
      N = 1 + z / BILLION;
      n = z / N;
      j2++;
      util_Assert (j2 <= ALPHABIT_NUM, "Alphabit:   j2 > ALPHABIT_NUM");
      if (n >= 30) {
         for (i = 0; i < Rep[j2]; ++i) {
            swalk_RandomWalk1 (gen, res, N, n, r, s, L, L);
            GetPVal_Walk (N, res, &j, " (L = 320)", j2);
         }
      }
      swalk_DeleteRes (res);
   }

   bbattery_NTests = ++j;
   if (blocFlag) {
      unif01_DeleteBitBlockGen (gen);
      gen = gen0;
   }
   if (fileFlag) {
      WriteReport (fname, "Alphabit", bbattery_NTests,
         bbattery_pVal, Timer, TRUE, TRUE, nb);
      ufile_DeleteReadBin (gen);
   } else {
      GetName (gen, genName);
      WriteReport (genName, "Alphabit", bbattery_NTests, bbattery_pVal,
         Timer, FALSE, TRUE, nb);
   }

   chrono_Delete (Timer);
  /*  t1 = time (NULL);
    WriteTime (t0, t1); */
}


/*=========================================================================*/

void bbattery_Alphabit (unif01_Gen * gen, double nb, int r, int s)
{
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= ALPHABIT_NUM; ++i)
      Rep[i] = 1;
   Alphabit (gen, NULL, nb, r, s, FALSE, 0, Rep);
}


/*=========================================================================*/

void bbattery_AlphabitFile (char *filename, double nb)
{
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= ALPHABIT_NUM; ++i)
      Rep[i] = 1;
   Alphabit (NULL, filename, nb, 0, 32, FALSE, 0, Rep);
}


/*=========================================================================*/

void bbattery_RepeatAlphabit (unif01_Gen * gen, double nb, int r, int s,
   int Rep[])
{
   Alphabit (gen, NULL, nb, r, s, FALSE, 0, Rep);
}


/*=========================================================================*/

void bbattery_BlockAlphabit (unif01_Gen * gen, double n, int r, int s)
{
   unif01_Gen *gen2;
   int L = 1;
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= ALPHABIT_NUM; ++i)
      Rep[i] = 1;
   while ((L <= 32) && (L <= s)) {
      gen2 = unif01_CreateBitBlockGen (gen, r, s, L);
      Alphabit (gen2, NULL, n, r, s, FALSE, 0, Rep);
      unif01_DeleteBitBlockGen (gen2);
      L *= 2;
   }
}


/*=========================================================================*/

void bbattery_RepeatBlockAlphabit (unif01_Gen * gen, double nb, int r, int s,
   int Rep[], int L)
{
   if ((L <= 32) && (L <= s)) {
      unif01_Gen *gen2;
      gen2 = unif01_CreateBitBlockGen (gen, r, s, L);
      Alphabit (gen2, NULL, nb, r, s, FALSE, 0, Rep);
      unif01_DeleteBitBlockGen (gen2);
   }
}


/*=========================================================================*/

void bbattery_BlockAlphabitFile (char *filename, double nb)
{
   int w = 1;
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= ALPHABIT_NUM; ++i)
      Rep[i] = 1;
   while (w <= 32) {
      Alphabit (NULL, filename, nb, 0, 32, TRUE, w, Rep);
      w *= 2;
   }
}


/*=========================================================================*/

static void DoMultinom (lebool fileFlag, /* */
   unif01_Gen * gen,              /* */
   double nb,                     /* Number of bits */
   int *pj,                       /* j */
   int j2,                        /* Test number in the battery */
   int Rep[]                      /* Number of replications */
   )
/*
 * Do the smultin_MultinomialBits in Rabbit
 */
{
   const long NLIM = 10000000;
   long n, N;
   int L, t;
   double x;
   int i;
   int j = *pj;
   smultin_Res *res;
   smultin_Param *par = NULL;
   double ValDelta[] = { -1 };

   util_Assert (nb > 0.0, "MultinomialBits:   nb <= 0");
   par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, -3);
   res = smultin_CreateRes (par);
   if (fileFlag)
      ufile_InitReadBin ();

#ifdef USE_LONGLONG
   /* Limit sample size n to NLIM because of memory limitations. */
   /* Determine number of replications N from this. */
   N = 1 + nb / NLIM;
   n = nb / N;
   /* Time limit on test: N = 30 */
   N = util_Min (30, N);
   /* Set n as a multiple of s = 32 */
   n -= n % 32;
   L = num_Log2 (n / 200.0 * n);
   L = util_Max (4, L);
   for (i = 0; i < Rep[j2]; ++i) {
      smultin_MultinomialBitsOver (gen, par, res, N, n, 0, 32, L, TRUE);
      strcpy (bbattery_TestNames[++j], "MultinomialBitsOver");
      bbattery_pVal[j] = res->pColl;
      TestNumber[j] = j2;
   }

#else
   x = nb / 32.0;
   N = 1 + x / NLIM;
   n = x / N;
   N = util_Min (30, N);
   L = 16;
   t = 32 / L;
   /* We want a number of collisions >= 2 */
   while ((L > 1) && (n / num_TwoExp[L] * n * t * t < 2.0)) {
      L /= 2;
      t = 32 / L;
   }
   n = n * (32 / L);
   /* We want a density n / k < 2 to use case Sparse = TRUE */
   if (n > 2 * num_TwoExp[L]) {
      N = n / num_TwoExp[L] * N;
      n /= N;
      while ((double) N * n * L > nb)
         n--;
   }
   while (n * L % 32 > 0)
      n--;
   if (n > 3) {
      for (i = 0; i < Rep[j2]; ++i) {
         smultin_MultinomialBits (gen, par, res, N, n, 0, 32, L, TRUE);
         strcpy (bbattery_TestNames[++j], "MultinomialBits");
         bbattery_pVal[j] = res->pColl;
         TestNumber[j] = j2;
      }
   }
#endif
   *pj = j;
   smultin_DeleteRes (res);
   smultin_DeleteParam (par);
}


/*-------------------------------------------------------------------------*/

static void DoAppear (lebool fileFlag, /* */
   unif01_Gen * gen, double nb,   /* Number of bits to test */
   int *pj,                       /* j */
   int j2,                        /* Test number in the battery */
   int Rep[]
   )
/*
 * Do the svaria_AppearanceSpacings test in Rabbit
 */
{
   sres_Basic *res;
   const long NLIM = 2000000000;
   int L;
   long N, Q;
   int i;
   int j = *pj;
   double temp = nb * (30.0 / 32.0) / 20.0;

   res = sres_CreateBasic ();
   if (num_TwoExp[30] < temp / 30.0)
      L = 30;
   else if (num_TwoExp[15] < temp / 15.0)
      L = 15;
   else if (num_TwoExp[10] < temp / 10.0)
      L = 10;
   else if (num_TwoExp[6] < temp / 6.0)
      L = 6;
   else if (num_TwoExp[5] < temp / 5.0)
      L = 5;
   else if (num_TwoExp[3] < temp / 3.0)
      L = 3;
   else
      L = 2;
   temp = nb / 2;
   temp *= 30.0 / 32.0;
   temp /= L;
   N = 1 + temp / NLIM;
   Q = temp / N;
   N = 1;

   if (Q < 50)
      return;
   if (fileFlag)
      ufile_InitReadBin ();
   for (i = 0; i < Rep[j2]; ++i) {
      svaria_AppearanceSpacings (gen, res, N, Q, Q, 0, 30, L);
      j++;
      if (N == 1)
         bbattery_pVal[j] = res->pVal2[gofw_Mean];
      else
         bbattery_pVal[j] = res->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "AppearanceSpacings");
   }
   sres_DeleteBasic (res);
   *pj = j;
}


/*-------------------------------------------------------------------------*/

static void DoWalk (lebool fileFlag, /* */
   unif01_Gen * gen,              /* */
   double nb,                     /* Number of bits to test */
   int *pj,                       /* j */
   int j2,                        /* Test number in the battery */
   int Rep[]
   )
/*
 * Do 3 swalk_RandomWalk1 tests in Rabbit
 */
{
   swalk_Res *res;
   long n, N, L;
   double z;
   int i;

   L = 128;
   z = nb / L;
   N = 1 + z / BILLION;
   n = z / N;
   N = 1;
   while (n < 100) {
      L /= 2;
      n *= 2;
   }
   if (L < 4)
      return;
   n = nb / (L * N);
   n = util_Min (n, 500 * MILLION);
   if (L < 32) {
      while (32 * n > nb)
         n--;
   }
   if (n < 30)
      return;

   res = swalk_CreateRes ();
   ++j2;
   if (fileFlag)
      ufile_InitReadBin ();
   for (i = 0; i < Rep[j2]; ++i) {
      swalk_RandomWalk1 (gen, res, N, n, 0, 32, L, L);
      GetPVal_Walk (N, res, pj, "", j2);
   }
   if (L < 96)
      return;

   L = 1024;
   z = nb / L;
   N = 1 + z / BILLION;
   n = z / N;
   n = util_Min (n, 50 * MILLION);
   N = 1;
   while ((double) n * L > nb)
      n--;
   if (n < 30)
      return;

   ++j2;
   if (fileFlag)
      ufile_InitReadBin ();
   for (i = 0; i < Rep[j2]; ++i) {
      swalk_RandomWalk1 (gen, res, N, n, 0, 32, L, L);
      GetPVal_Walk (N, res, pj, " (L = 1024)", j2);
   }

   L = 10016;
   z = nb / L;
   N = 1 + z / BILLION;
   n = z / N;
   n = util_Min (n, 5 * MILLION);
   N = 1;
   while ((double) n * L > nb)
      n--;
   if (n < 30)
      return;
   ++j2;
   if (fileFlag)
      ufile_InitReadBin ();
   for (i = 0; i < Rep[j2]; ++i) {
      swalk_RandomWalk1 (gen, res, N, n, 0, 32, L, L);
      GetPVal_Walk (N, res, pj, " (L = 10016)", j2);
   }

   swalk_DeleteRes (res);
}


/*-------------------------------------------------------------------------*/

static void Rabbit (unif01_Gen * gen, char *fname, double nb, int Rep[])
/*
 * A battery of statistical tests for a file of n random bits.
 */
{
   const int s = 32;
   int k, j = 0, j2 = 0;
   int i;
   long n, N, L;
   double nw, x;
   chrono_Chrono *Timer;
   long bufsiz;
   lebool fileFlag;
   char genName[LEN + 1] = "";

   Timer = chrono_Create ();
   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "          Starting Rabbit:   nb = %.0f\n"
         "          Version: %s\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n",
         nb, PACKAGE_STRING);
   }
   util_Assert (nb >= 500.0, "bbattery_Rabbit:   nb < 500");

   /* Bits will be read as 32-bit unsigned integers */
   nb -= fmod (nb, 32.0);
   nw = nb / 32.0;
   bufsiz = nw;

   if (NULL == gen) {
      gen = ufile_CreateReadBin (fname, bufsiz);
      fileFlag = TRUE;
   } else
      fileFlag = FALSE;

   j = -1;
   ++j2;
   DoMultinom (fileFlag, gen, nb, &j, j2, Rep);

   {
      const long NLIM = 4000000;
      snpair_Res *res;
      res = snpair_CreateRes ();
      N = 1 + nw / NLIM;
      n = nw / N;
      N = util_Min (N, 25);
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairsBitMatch (gen, res, N, n / 2, 0, 2);
         bbattery_pVal[++j] = res->pVal[snpair_BM];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "ClosePairsBitMatch, t = 2");
      }

      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         snpair_ClosePairsBitMatch (gen, res, N, n / 4, 0, 4);
         bbattery_pVal[++j] = res->pVal[snpair_BM];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "ClosePairsBitMatch, t = 4");
      }
      snpair_DeleteRes (res);
   }

   ++j2;
   DoAppear (fileFlag, gen, nb, &j, j2, Rep);

   {
      const long NLIM1 = 300000;
      const long NLIM2 = 10000;
      scomp_Res *res;
      res = scomp_CreateRes ();
      n = NLIM2 + 2.0 * sqrt (nb);
      n = util_Min (n, nb);
      n = util_Min (n, NLIM1);
      N = 1;
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LinearComp (gen, res, N, n, 0, s);
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->JumpSize->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->JumpSize->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp");
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->JumpNum->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->JumpNum->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LinearComp");
      }
      scomp_DeleteRes (res);
   }

   k = num_Log2 (nb + 0.5);
   if (k > 28)
      k = 28;
   N = 1;
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         scomp_LempelZiv (gen, res, N, k, 0, s);
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "LempelZiv");
      }
      sres_DeleteBasic (res);
   }
   {
      sspectral_Res *res;
      k = num_Log2 (nb + 0.5);
      k = util_Min (20, k);
      res = sspectral_CreateRes ();
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sspectral_Fourier1 (gen, res, 1, k, 0, s);
         j++;
         bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Fourier1");
      }

      x = sqrt (2.0 * nb);
      N = x / 2.0;
      if (N < 32) {
         k = 5;
         N = nb / 32.0;
      } else if (N >= 16384) {
         k = 14;
         N = nb / 16384.0;
      } else {
         k = num_Log2 (x / 2.0 + 0.5);
         N = nb / (num_TwoExp[k]);
      }
      N = util_Min (N, 300000);
      while ((num_TwoExp[k] + 32) * N > nb)
         N--;
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sspectral_Fourier3 (gen, res, N, k, 0, s);
         j++;
         bbattery_pVal[j] = res->Bas->pVal2[gofw_AD];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Fourier3");
      }
      sspectral_DeleteRes (res);
   }
   {
      sstring_Res2 *res;
      res = sstring_CreateRes2 ();
      x = util_Min (BILLION * 100.0, nb);
      n = 600;
      L = x / n;
      if (L <= 100000) {
         n /= 10;
         L *= 10;
      }
      if (L <= 10000) {
         n /= 2;
         L *= 2;
      }
      ++j2;
      if ((L >= 1032) && (n >= 30)) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            sstring_LongestHeadRun (gen, res, 1, n, 0, s, L);
            j++;
            bbattery_pVal[j] = res->Chi->pVal2[gofw_Mean];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "LongestHeadRun");
         }
      }
      sstring_DeleteRes2 (res);
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      nw = nb / 32.0;
      nw = util_Min (nw, 4.0 * BILLION);
      N = 1 + nw / BILLION;
      n = nw / N;
      ++j2;
      if (n >= 30) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            sstring_PeriodsInStrings (gen, res, N, n, 0, 31);
            ++j;
            if (N == 1)
               bbattery_pVal[j] = res->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->pVal2[gofw_Sum];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "PeriodsInStrings");
         }
      }

      nw = nb / s;
      N = 1 + nw / BILLION;
      n = nw / N;
      N = util_Min (10, N);
      ++j2;
      if (n > 29) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            sstring_HammingWeight (gen, res, N, n, 0, s, s);
            ++j;
            if (N == 1)
               bbattery_pVal[j] = res->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->pVal2[gofw_Sum];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "HammingWeight");
         }
      }
      sres_DeleteChi2 (res);
   }
   {
      sstring_Res *res;
      res = sstring_CreateRes ();
      nw = nb / s;
      N = 1 + nw / BILLION;
      n = nw / N;
      N = util_Min (10, N);
      ++j2;
      if (n > 2) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            sstring_HammingCorr (gen, res, N, n, 0, s, 32);
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "HammingCorr, L = 32");
         }
      }

      nw = nb / 64;
      N = 1 + nw / BILLION;
      n = nw / N;
      N = 1;
      ++j2;
      if (n > 2) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            sstring_HammingCorr (gen, res, N, n, 0, s, 2 * s);
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "HammingCorr, L = 64");
         }
      }

      nw = nb / (4 * s);
      N = 1 + nw / BILLION * 4;
      n = nw / N;
      N = 1;
      ++j2;
      if (n > 2) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            sstring_HammingCorr (gen, res, N, n, 0, s, 4 * s);
            j++;
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "HammingCorr, L = 128");
         }
      }

      nw = nb / s;
      N = 1 + nw / BILLION;
      n = nw / N;
      N = util_Min (5, N);
      ++j2;
      if (n > 29) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            sstring_HammingIndep (gen, res, N, n, 0, s, 16, 0);
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "HammingIndep, L = 16");
         }
      }

      nw = nb / (2 * s);
      N = 1 + nw / BILLION * 2;
      n = nw / N;
      N = 1;
      ++j2;
      if (n > 29) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            sstring_HammingIndep (gen, res, N, n, 0, s, s, 0);
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "HammingIndep, L = 32");
         }
      }

      nw = nb / (4 * s);
      N = 1 + nw / BILLION * 10;
      n = nw / N;
      N = 1;
      ++j2;
      if (n > 29) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            sstring_HammingIndep (gen, res, N, n, 0, s, 2 * s, 0);
            if (N == 1)
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Mean];
            else
               bbattery_pVal[j] = res->Bas->pVal2[gofw_Sum];
            strcpy (bbattery_TestNames[j], "HammingIndep, L = 64");
            TestNumber[j] = j2;
         }
      }
      sstring_DeleteRes (res);
   }
   {
      sres_Basic *res;
      int d;
      res = sres_CreateBasic ();

      d = 1;
      N = 1 + nb / BILLION;
      n = nb / N - d;
      n -= n % 32;
      N = util_Min (100, N);
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, N, n, 0, s, d);
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor");
      }

      d = 2;
      N = 1 + nb / BILLION;
      n = nb / N - d;
      n -= n % 32;
      N = util_Min (100, N);
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_AutoCor (gen, res, N, n, 0, s, d);
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "AutoCor");
      }

      sres_DeleteBasic (res);
   }
   {
      sstring_Res3 *res;
      res = sstring_CreateRes3 ();
      nw = nb / 5;
      N = 1 + nw / BILLION;
      n = nw / N;
      N = util_Min (20, N);
      if (fileFlag)
         ufile_InitReadBin ();
      ++j2;
      for (i = 0; i < Rep[j2]; ++i) {
         sstring_Run (gen, res, N, n, 0, s);
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->NRuns->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->NRuns->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits");
         j++;
         if (N == 1)
            bbattery_pVal[j] = res->NBits->pVal2[gofw_Mean];
         else
            bbattery_pVal[j] = res->NBits->pVal2[gofw_Sum];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "Run of bits");
       }
      sstring_DeleteRes3 (res);
   }

   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      n = nb / (s * s);
      n = util_Min (n, 50 * MILLION);
      ++j2;
      if (n >= 50) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            smarsa_MatrixRank (gen, res, 1, n, 0, s, s, s);
            bbattery_pVal[j] = res->pVal2[gofw_Mean];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "MatrixRank, 32 x 32");
         }
      }

      n = nb / (100.0 * s * s);
      n = util_Min (n, 300000);
      ++j2;
      if (n >= 50) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            smarsa_MatrixRank (gen, res, 1, n, 0, s, 10 * s, 10 * s);
            bbattery_pVal[j] = res->pVal2[gofw_Mean];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "MatrixRank, 320 x 320");
         }
      }

      n = nb / (1024.0 * s * s);
      n = util_Min (n, 20000);
      ++j2;
      if (n >= 50) {
         if (fileFlag)
            ufile_InitReadBin ();
         for (i = 0; i < Rep[j2]; ++i) {
            j++;
            smarsa_MatrixRank (gen, res, 1, n, 0, s, 32 * s, 32 * s);
            bbattery_pVal[j] = res->pVal2[gofw_Mean];
            TestNumber[j] = j2;
            strcpy (bbattery_TestNames[j], "MatrixRank, 1024 x 1024");
         }
      }
      sres_DeleteChi2 (res);
   }

   DoWalk (fileFlag, gen, nb, &j, j2, Rep);
   util_Assert (j2 <= RABBIT_NUM, "Rabbit:   j2 > RABBIT_NUM");

   bbattery_NTests = ++j;
   if (fileFlag) {
      WriteReport (fname, "Rabbit", bbattery_NTests,
         bbattery_pVal, Timer, TRUE, TRUE, nb);
      ufile_DeleteReadBin (gen);
   } else {
      GetName (gen, genName);
      WriteReport (genName, "Rabbit", bbattery_NTests, bbattery_pVal,
         Timer, FALSE, TRUE, nb);
   }
   chrono_Delete (Timer);
}


/*=========================================================================*/

void bbattery_Rabbit (unif01_Gen * gen, double nb)
{
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= RABBIT_NUM; ++i)
      Rep[i] = 1;
   Rabbit (gen, NULL, nb, Rep);
}


/*=========================================================================*/

void bbattery_RabbitFile (char *filename, double nb)
{
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= RABBIT_NUM; ++i)
      Rep[i] = 1;
   Rabbit (NULL, filename, nb, Rep);
}


/*=========================================================================*/

void bbattery_RepeatRabbit (unif01_Gen * gen, double nb, int Rep[])
{
   Rabbit (gen, NULL, nb, Rep);
}


/*=========================================================================*/

void bbattery_pseudoDIEHARD (unif01_Gen * gen)
/*
 * As close as possible to the DIEHARD test suite.
 */
{
   chrono_Chrono *Timer;
   smultin_Param *par = NULL;
   double ValDelta[] = { 1 };
   char genName[LEN + 1] = "";
   int k, i, j = -1;
   int j2 = 0;
   double x;
   long Count[7];
   double NumExp[7] = {
      67.668, 135.335, 135.335, 90.224, 45.112, 18.045, 8.282
   };

   Timer = chrono_Create ();
   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "                 Starting pseudoDIEHARD\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n");
   }
   {
      sres_Poisson *res;
      sres_Chi2 *Chi;
      Chi = sres_CreateChi2 ();
      sres_InitChi2 (Chi, 1, 6, "");
      res = sres_CreatePoisson ();
      printf ("smarsa_BirthdaySpacings test with r = 0, 1, 2, 3, 4, 5,"
         " 6, 7, 8,\n .....\n\n");
      swrite_Basic = FALSE;
      ++j2;
      for (i = 0; i <= 8; i++) {
         printf (" r = %d\n", i);
         for (k = 0; k <= 6; k++)
            Count[k] = 0;
         for (k = 0; k < 500; k++) {
            smarsa_BirthdaySpacings (gen, res, 1, 512, i, 16777216, 1, 1);
            if (res->sVal2 >= 6)
               ++Count[6];
            else
               ++Count[(int) res->sVal2];
         }
         x = gofs_Chi2 (NumExp, Count, 0, 6);
         printf ("ChiSquare statistic                   :");
         bbattery_pVal[++j] = fbar_ChiSquare2 (6, 12, x);
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "BirthdaySpacings");
         gofw_Writep2 (x, bbattery_pVal[j]);
      }
      printf ("\n\n\n\n");
      sres_DeletePoisson (res);
      sres_DeleteChi2 (Chi);
      swrite_Basic = TRUE;
   }
   ++j2;
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      smarsa_MatrixRank (gen, res, 1, 40000, 0, 31, 31, 31);
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = ++j2;
      strcpy (bbattery_TestNames[j], "MatrixRank");

      smarsa_MatrixRank (gen, res, 1, 40000, 0, 32, 32, 32);
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank");

      for (i = 0; i <= 24; i++) {
         smarsa_MatrixRank (gen, res, 1, 100000, i, 8, 6, 8);
         bbattery_pVal[++j] = res->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "MatrixRank");
      }
      sres_DeleteChi2 (res);
   }
   {
      smultin_Res *res;
      par = smultin_CreateParam (1, ValDelta, smultin_GenerCellSerial, 0);
      res = smultin_CreateRes (par);
      smultin_MultinomialBitsOver (gen, par, res, 20, 2097152, 0, 32, 20,
         TRUE);
      bbattery_pVal[++j] = res->pVal2[0][gofw_AD];
      TestNumber[j] = ++j2;
      strcpy (bbattery_TestNames[j], "MultinomialBitsOver");
      smultin_DeleteRes (res);
      smultin_DeleteParam (par);
   }
   {
      smarsa_Res *res;
      res = smarsa_CreateRes ();
      ++j2;
      for (i = 22; i >= 0; i--) {
         smarsa_Opso (gen, res, 1, i, 1);
         bbattery_pVal[++j] = res->Pois->pVal2;
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "OPSO");
      }
      ValDelta[0] = -1.0;
      ++j2;
      for (i = 27; i >= 0; i--) {
         if (swrite_Basic)
            printf ("***********************************************************\n"
               "Test OQSO calling smarsa_CollisionOver\n\n");
         smarsa_CollisionOver (gen, res, 1, 2097152, i, 32, 4);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "OQSO");
      }
      ++j2;
      for (i = 30; i >= 0; i--) {
         if (swrite_Basic)
            printf ("***********************************************************\n"
               "Test DNA calling smarsa_CollisionOver\n\n");
         smarsa_CollisionOver (gen, res, 1, 2097152, i, 4, 10);
         bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
         TestNumber[j] = j2;
         strcpy (bbattery_TestNames[j], "DNA");
      }
      smarsa_DeleteRes (res);
   }
   j2 += 2;
   {
      snpair_Res *res;
      res = snpair_CreateRes ();
      snpair_ClosePairs (gen, res, 100, 8000, 0, 2, 2, 1);
      bbattery_pVal[++j] = res->pVal[snpair_NP];
      TestNumber[j] = ++j2;
      strcpy (bbattery_TestNames[j], "ClosePairs");

      snpair_ClosePairs (gen, res, 20, 4000, 0, 3, 2, 1);
      bbattery_pVal[++j] = res->pVal[snpair_NP];
      TestNumber[j] = ++j2;
      strcpy (bbattery_TestNames[j], "ClosePairs");
      snpair_DeleteRes (res);
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      smarsa_Savir2 (gen, res, 1, 100000, 0, 90000, 18);
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = ++j2;
      strcpy (bbattery_TestNames[j], "Savir2");

      ++j2;
      sknuth_Run (gen, res, 10, 10000, 0, TRUE);
      bbattery_pVal[++j] = res->pVal2[gofw_Sum];
      TestNumber[j] = ++j2;
      strcpy (bbattery_TestNames[j], "Run of U01");

      sknuth_Run (gen, res, 10, 10000, 0, FALSE);
      bbattery_pVal[++j] = res->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of U01");

      sknuth_Run (gen, res, 10, 10000, 0, TRUE);
      bbattery_pVal[++j] = res->pVal2[gofw_Sum];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of U01");

      sknuth_Run (gen, res, 10, 10000, 0, FALSE);
      bbattery_pVal[++j] = res->pVal2[gofw_Sum];
      strcpy (bbattery_TestNames[j], "Run of U01");
      TestNumber[j] = j2;
      sres_DeleteChi2 (res);
   }

   bbattery_NTests = ++j;
   GetName (gen, genName);
   WriteReport (genName, "pseudoDIEHARD", bbattery_NTests, bbattery_pVal,
      Timer, FALSE, FALSE, 0.0);
   chrono_Delete (Timer);
}


/*=========================================================================*/

static double ProbabiliteLHR (long j, double Lnl)
/*
 * Returns the probability that the longest series of successive 1 has
 * length = j.
 */
{
   double x, temp;
   temp = (j + 1) * num_Ln2 - Lnl;
   x = exp (-exp (-temp));
   temp += num_Ln2;
   x = exp (-exp (-temp)) - x;
   return x;
}

/*-------------------------------------------------------------------------*/

static double GetPLongest (int longest)
/*
 * Get the probabilities for the longest run of 1 or 0 over 20000 bits.
 */
{
   double pLeft, pRight;
   double LnLen;
   int j;

   LnLen = log (20000.0);
   pLeft = 0.0;
   for (j = 0; j < longest; j++)
      pLeft += ProbabiliteLHR (j, LnLen);
   pRight = 1.0 - pLeft;
   pLeft += ProbabiliteLHR (longest, LnLen);
   return gofw_pDisc (pLeft, pRight);
}


/*-------------------------------------------------------------------------*/

static void WriteReportFIPS_140_2 (
   char *genName,                 /* Generator or file name */
   lebool Flag,                  /* = TRUE for a file, FALSE for a gen */
   int nbit,                      /* Number of bits */
   int longest0,                  /* Longest string of 0 */
   int longest1,                  /* Longest string of 1 */
   int nrun0[],                   /* Number of 0 runs */
   int nrun1[],                   /* Number of 1 runs */
   int ncount[]                   /* Number of 4 bits values */
   )
{
   int i, j;
   double X;
   fmass_INFO Q;
   double p, pLeft, pRight;
   lebool failFlag = FALSE;

   printf
      ("\n============== Summary results of FIPS-140-2 ==============\n\n");
   if (Flag) {
      printf (" File:             ");
   } else {
      printf (" Generator:        ");
   }
   printf ("%s", genName);
   printf ("\n Number of bits:   20000\n");

   printf ("\n       Test          s-value        p-value    FIPS Decision\n");
   printf (" --------------------------------------------------------\n");

   /* Monobit results */
   j = 0;
   printf (" %-20s", bbattery_TestNames[j]);
   printf (" %5d       ", nbit);
   Q = fmass_CreateBinomial (20000, 0.5, 0.5);
   pLeft = fdist_Binomial2 (Q, nbit);
   pRight = fbar_Binomial2 (Q, nbit);
   fmass_DeleteBinomial (Q);
   p = gofw_pDisc (pLeft, pRight);
   gofw_Writep0 (p);
   if ((nbit <= 9725) || nbit >= 10275) {
      printf (" %10s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %10s", "Pass");

   printf ("\n");

   /* Poker results */
   X = 0;
   for (i = 0; i < 16; i++)
      X += (double) ncount[i] * ncount[i];
   X = 16 * X / 5000 - 5000;
   j = 1;
   printf (" %-16s", bbattery_TestNames[j]);
   printf ("%10.2f       ", X);
   p = fbar_ChiSquare2 (15, 12, X);
   gofw_Writep0 (p);
   if ((X <= 2.16) || X >= 46.17) {
      printf (" %10s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %10s", "Pass");
   printf ("\n\n");

   /* Run results */
   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun0[1]);
   if ((nrun0[1] <= 2315) || nrun0[1] >= 2685) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun0[2]);
   if ((nrun0[2] <= 1114) || nrun0[2] >= 1386) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun0[3]);
   if ((nrun0[3] <= 527) || nrun0[3] >= 723) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun0[4]);
   if ((nrun0[4] <= 240) || nrun0[4] >= 384) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun0[5]);
   if ((nrun0[5] <= 103) || nrun0[5] >= 209) {
      failFlag = TRUE;
      printf (" %25s", "Fail");
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun0[6]);
   if ((nrun0[6] <= 103) || nrun0[6] >= 209) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun1[1]);
   if ((nrun1[1] <= 2315) || nrun1[1] >= 2685) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun1[2]);
   if ((nrun1[2] <= 1114) || nrun1[2] >= 1386) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun1[3]);
   if ((nrun1[3] <= 527) || nrun1[3] >= 723) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun1[4]);
   if ((nrun1[4] <= 240) || nrun1[4] >= 384) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun1[5]);
   if ((nrun1[5] <= 103) || nrun1[5] >= 209) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d", nrun1[6]);
   if ((nrun1[6] <= 103) || nrun1[6] >= 209) {
      printf (" %25s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %25s", "Pass");
   printf ("\n\n");

   /* Longest run results */
   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d       ", longest0);
   p = GetPLongest (longest0);
   gofw_Writep0 (p);
   if (longest0 >= 26) {
      printf (" %10s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %10s", "Pass");
   printf ("\n");

   printf (" %-20s", bbattery_TestNames[++j]);
   printf (" %5d       ", longest1);
   p = GetPLongest (longest1);
   gofw_Writep0 (p);
   if (longest1 >= 26) {
      printf (" %10s", "Fail");
      failFlag = TRUE;
   } else
      printf (" %10s", "Pass");
   printf ("\n");

   if (!failFlag) {
      printf (" ----------------------------------------------------------\n");
      printf (" All values are within the required intervals of FIPS-140-2\n");
   }
   printf ("\n\n\n");
}


/*-------------------------------------------------------------------------*/

#define SAMPLE 625                /* 625 * 32 = 20000 */
#define MASK4  15                 /* Mask of 4 bits */

static void FIPS_140_2 (unif01_Gen * gen, char *filename)
{
   int i, j;
   int nbit = 0;                  /* Number of bits */
   int longest0 = 0;              /* Longest string of 0 */
   int longest1 = 0;              /* Longest string of 1 */
   int nrun0[7] = { 0 };          /* Number of 0 runs */
   int nrun1[7] = { 0 };          /* Number of 1 runs */
   int ncount[16] = { 0 };        /* Number of 4 bits values */
   int prevBit;                   /* Previous bit */
   int len = 0;                   /* Length of run */
   unsigned long jBit;            /* Current bit */
   unsigned long Z;               /* Block of 32 bits */
   unsigned long Bits[SAMPLE + 1];
   lebool fileFlag = FALSE;
   char genName[LEN + 1] = "";

   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "                 Starting FIPS_140_2\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n");
   }
   util_Assert (NULL == gen || NULL == filename,
      "bbattery_FIPS_140_2:   one of gen or filename must be NULL");
   util_Assert (!(NULL == gen && NULL == filename),
      "bbattery_FIPS_140_2:   no generator and no file");
   util_Assert (!(NULL == gen && !(strcmp (filename, ""))),
      "bbattery_FIPS_140_2:   no generator and no file");

   if ((NULL == gen) && filename && strcmp (filename, "")) {
      gen = ufile_CreateReadBin (filename, SAMPLE);
      fileFlag = TRUE;
   }

   for (j = 0; j < SAMPLE; j++)
      Bits[j] = unif01_StripB (gen, 0, 32);

   if (fileFlag) {
      ufile_DeleteReadBin (gen);
      strncpy (genName, filename, (size_t) LEN);
   } else {
      GetName (gen, genName);
   }

   /* Make sure to count the first run; set prevBit != {0, 1} */
   prevBit = 2;

   for (j = 0; j < SAMPLE; j++) {
      /* Count the number of 1 */
      Z = Bits[j];
      while (Z > 0) {
         Z &= Z - 1;              /* Clear lowest 1 bit */
         ++nbit;
      }

      /* Count the number of 4 bits values */
      Z = Bits[j];
      for (i = 0; i < 8; i++) {
         (ncount[Z & MASK4])++;
         Z >>= 4;
      }

      /* Count the number of runs and get the longest runs */
      Z = Bits[j];
      jBit = bitset_maskUL[31];

      while (jBit > 0) {
         if (Z & jBit) {          /* bit 1 */
            if (prevBit != 1) {
               if (len < 6)
                  (nrun0[len])++;
               else
                  (nrun0[6])++;
               if (len > longest0)
                  longest0 = len;
               len = 1;
            } else {
               len++;
            }
            prevBit = 1;

         } else {                 /* bit 0 */
            if (prevBit != 0) {
               if (len < 6)
                  (nrun1[len])++;
               else
                  (nrun1[6])++;
               if (len > longest1)
                  longest1 = len;
               len = 1;
            } else {
               len++;
            }
            prevBit = 0;
         }
         jBit >>= 1;
      }
   }

   strcpy (bbattery_TestNames[0], "Monobit");
   strcpy (bbattery_TestNames[1], "Poker");
   j = 1;
   strcpy (bbattery_TestNames[++j], "0 Runs, length 1: ");
   strcpy (bbattery_TestNames[++j], "0 Runs, length 2: ");
   strcpy (bbattery_TestNames[++j], "0 Runs, length 3: ");
   strcpy (bbattery_TestNames[++j], "0 Runs, length 4: ");
   strcpy (bbattery_TestNames[++j], "0 Runs, length 5: ");
   strcpy (bbattery_TestNames[++j], "0 Runs, length 6+: ");
   strcpy (bbattery_TestNames[++j], "1 Runs, length 1: ");
   strcpy (bbattery_TestNames[++j], "1 Runs, length 2: ");
   strcpy (bbattery_TestNames[++j], "1 Runs, length 3: ");
   strcpy (bbattery_TestNames[++j], "1 Runs, length 4: ");
   strcpy (bbattery_TestNames[++j], "1 Runs, length 5: ");
   strcpy (bbattery_TestNames[++j], "1 Runs, length 6+: ");

   strcpy (bbattery_TestNames[++j], "Longest run of 0: ");
   strcpy (bbattery_TestNames[++j], "Longest run of 1: ");

   WriteReportFIPS_140_2 (genName, fileFlag, nbit, longest0, longest1,
      nrun0, nrun1, ncount);
}


/*-------------------------------------------------------------------------*/

void bbattery_FIPS_140_2 (unif01_Gen * gen)
{
   FIPS_140_2 (gen, NULL);
}


/*-------------------------------------------------------------------------*/

void bbattery_FIPS_140_2File (char *filename)
{
   FIPS_140_2 (NULL, filename);
}


/*=========================================================================*/
