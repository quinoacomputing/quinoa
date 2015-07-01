/*************************************************************************\
 *
 * Package:        TestU01
 * File:           fcong.c
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
#include "num.h"

#include "fcong.h"
#include "ffam.h"
#include "unif01.h"
#include "ulcg.h"
#include "uinv.h"
#include "ucubic.h"
#include "umrg.h"

#include <string.h>


#define MAXCAR  256                /* Max length of a line */
#define SEED    123                /* Default seed for generators */
#define BIGSEED "12345"            /* Default seed for big generators */
#define BIGCAR  "1024"        /* Max number of decimals in a big integer */


/* The kinds of predefined families */
typedef enum {
   LCG_a,
   LCGPow2_a,
   MRG2_a,
   MRG3_a,
   CombL2_a,
   CombWH2_a,
   InvImpl_a,
   InvImpl2a_a,
   InvImpl2b_a,
   InvExpl_a,
   InvExpl2a_a,
   InvExpl2b_a,
   InvMRG2_a,
   Cubic1_a,
   CombCubic2_a,
   CombCubLCG_a,
   GenN_a
}
GenType;



/*=========================================================================*/

static int ReadOneLCG (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m, a;

   /* Pass over LSize; it has been read already */
   if (fam->LSize[i] <= 31) {
      status = sscanf (Line, "%*d %ld %ld", &m, &a);
      util_Assert (status == 2,
         "ReadOneLCG:  Error in reading LCG parameter file");
      /* We use the float version since the a satisfy a*m < 2^{53} */
      fam->Gen[i] = ulcg_CreateLCGFloat (m, a, 0, SEED % m);

   } else {
#ifndef USE_GMP
      return -1;
#else
      char bigm[MAXCAR + 1];
      char biga[MAXCAR + 1];
      char format[20];

      /* Build reading format for big integers read as a string */
      strncpy (format, "%*d %", (size_t) 6);
      strncat (format, BIGCAR, (size_t) 10);
      strncat (format, "s %", (size_t) 5);
      strncat (format, BIGCAR, (size_t) 10);
      strncat (format, "s", (size_t) 2);
      status = sscanf (Line, format, bigm, biga);
      util_Assert (status == 2,
         "ReadOneLCG:   Error in reading LCG (Big) parameter file");
      fam->Gen[i] = ulcg_CreateBigLCG (bigm, biga, "0", BIGSEED);
#endif
   }
   return 0;
}


/*=========================================================================*/

static int ReadOneLCGPow2 (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long a;
   long mask = num_TwoExp[fam->LSize[i]] - 1.0;

   /* Pass over LSize; it has been read already */
   if (fam->LSize[i] <= 31) {
      status = sscanf (Line, "%*d %ld", &a);
      util_Assert (status == 1,
         "ReadOneLCG:   Error in reading LCGPow2 parameter file");
      fam->Gen[i] = ulcg_CreatePow2LCG (fam->LSize[i], a, 1, SEED & mask);

   } else {
#ifndef USE_GMP
      return -1;
#else
      char biga[MAXCAR + 1];
      char format[20];

      /* Build reading format for big integers read as a string */
      strncpy (format, "%*d %", (size_t) 6);
      strncat (format, BIGCAR, (size_t) 10);
      strncat (format, "s", (size_t) 2);
      status = sscanf (Line, format, biga);
      util_Assert (status == 1,
         "ReadOneLCG:   Error in reading LCGPow2 (Big) parameter file");
      fam->Gen[i] =
         ulcg_CreateBigPow2LCG (fam->LSize[i], biga, "1", BIGSEED);
#endif
   }
   return 0;
}


/*=========================================================================*/

static int ReadOneMRG2 (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m, A[2], S[2] = { SEED, SEED };

   /* Pass over LSize; it has been read already */
   status = sscanf (Line, "%*d %ld %ld %ld", &m, &A[0], &A[1]);
   util_Assert (status == 3,
      "ReadOneMRG2:   Error in reading MRG2 parameter file");
   S[0] %= m;
   S[1] %= m;
   fam->Gen[i] = umrg_CreateMRGFloat (m, 2, A, S);
   return 0;
}


/*=========================================================================*/

static int ReadOneInvMRG2 (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m, A[2], S[2] = { SEED, SEED };

   /* Pass over LSize; it has been read already */
   status = sscanf (Line, "%*d %ld %ld %ld", &m, &A[0], &A[1]);
   util_Assert (status == 3,
      "ReadOneInvMRG2:   Error in reading MRG2 parameter file");
   S[0] %= m;
   S[1] %= m;
   fam->Gen[i] = uinv_CreateInvMRGFloat (m, 2, A, S);
   return 0;
}


/*=========================================================================*/

static int ReadOneMRG3 (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m, A[3], S[3] = { SEED, SEED, SEED };

   /* Pass over LSize; it has been read already */
   status = sscanf (Line, "%*d %ld %ld %ld %ld", &m, &A[0], &A[1], &A[2]);
   util_Assert (status == 4,
      "ReadOneMRG3:   Error in reading MRG3 parameter file");
   S[0] %= m;
   S[1] %= m;
   S[2] %= m;
   fam->Gen[i] = umrg_CreateMRGFloat (m, 3, A, S);
   return 0;
}


/*=========================================================================*/

static int ReadOneCombL2 (char *Line, ffam_Fam * fam, int i)
{
   int status, LSize;
   long m1, m2, a1, a2;

   status = sscanf (Line, "%d %ld %ld %ld %ld", &LSize, &m1, &m2, &a1, &a2);
   util_Assert (status == 5,
      "ReadOneCombL2:   Error in reading CombL2 parameter file");
   if (LSize < 54)
      fam->Gen[i] = ulcg_CreateCombLEC2Float (m1, m2, a1, a2, 0, 0,
                    SEED % m2, SEED % m2);
   else
      fam->Gen[i] = ulcg_CreateCombLEC2 (m1, m2, a1, a2, 0, 0,
                    SEED % m2, SEED % m2);

   return 0;
}


/*=========================================================================*/

static int ReadOneCombWH2 (char *Line, ffam_Fam * fam, int i)
{
   int status, LSize;
   long m1, m2, a1, a2;

   status = sscanf (Line, "%d %ld %ld %ld %ld", &LSize, &m1, &m2, &a1, &a2);
   util_Assert (status == 5,
      "ReadOneCombWH2:   Error in reading CombWH2 parameter file");
   if (LSize < 54)
      fam->Gen[i] = ulcg_CreateCombWH2Float (m1, m2, a1, a2, 0, 0,
                    SEED % m2, SEED % m2);
   else
      fam->Gen[i] = ulcg_CreateCombWH2 (m1, m2, a1, a2, 0, 0,
                    SEED % m2, SEED % m2);

   return 0;
}


/*=========================================================================*/

static int ReadOneInvImpl (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m, a1, a2;

   /* Pass over LSize; it has been read already */
   status = sscanf (Line, "%*d %ld %ld %ld", &m, &a1, &a2);
   util_Assert (status == 3,
      "ReadOneInvImpl:   Error in reading InvImpl parameter file");
   fam->Gen[i] = uinv_CreateInvImpl (m, a1, a2, SEED % m);
   return 0;
}


/*=========================================================================*/

static int ReadOneInvExpl (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m;

   /* Pass over LSize; it has been read already */
   status = sscanf (Line, "%*d %ld", &m);
   util_Assert (status == 1,
      "ReadOneInvExpl:   Error in reading InvExpl parameter file");
   fam->Gen[i] = uinv_CreateInvExpl (m, SEED % m, 0);
   return 0;
}


/*=========================================================================*/

static int ReadOneCubic1 (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m, a;

   status = sscanf (Line, "%*d %ld %ld", &m, &a);
   util_Assert (status == 2,
      "ReadOneCubic1:   Error in reading Cubic1 parameter file");
   fam->Gen[i] = ucubic_CreateCubic1Float (m, a, SEED % m);
   return 0;
}


/*=========================================================================*/

static int ReadOneCombCubic2 (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m1, m2, a1, a2;

   status = sscanf (Line, "%*d %ld %ld %ld %ld", &m1, &m2, &a1, &a2);
   util_Assert (status == 4,
      "ReadOneCombCubic2:   Error in reading CombCubic2 parameter file");
   fam->Gen[i] = ucubic_CreateCombCubic2 (m1, m2, a1, a2,
                    SEED % m2, SEED % m2);
   return 0;
}


/*=========================================================================*/

static int ReadOneCombCubLCG (char *Line, ffam_Fam * fam, int i)
{
   int status;
   long m1, m2, a1, a2;
   unif01_Gen *gen1, *gen2;

   status = sscanf (Line, "%*d %ld %ld %ld %ld", &m1, &a1, &m2, &a2);
   util_Assert (status == 4,
      "ReadOneCombCubLCG:   Error in reading CombCubLCG parameter file");
   gen1 = ulcg_CreateLCGFloat (m1, a1, 0, SEED % m1);
   gen2 = ucubic_CreateCubic1Float (m2, a2, SEED % m2);
   fam->Gen[i] = unif01_CreateCombAdd2 (gen1, gen2, "CombCubLCG");
   return 0;
}



/*=========================================================================*/

static ffam_Fam *ReadAllGen (char *filename, char *deffile, GenType g,
   int i1, int i2, int istep)
/*
 * Read parameter file. If filename exists, it will be read; otherwise, the
 * default file deffile will be read. 
 * Keeps only generators whose LSize are in [i1, i2], and keep only generators
 * spaced istep apart. If i1 < the smallest LSize in the file, it will be
 * reset to the first LSize in the file; similarly if i2 > the largest LSize
 * in the file, it will be reset to the last LSize in the file.
 */
{
   FILE *f;
   char Line[MAXCAR + 1];
   ffam_Fam *fam;
   int i, j;
   int status, LSize;

   f = ffam_OpenFile (filename, deffile);

   /* Read the family name */
   util_GetLine (f, Line, '#');
   fam = ffam_CreateFam ((i2 - i1 + istep) / istep, Line);

   /* Read the parameters */
   i = j = LSize = 0;
   while (LSize <= i2) {
      status = util_GetLine (f, Line, '#');
      if (status)                 /* if EOF or error */
         break;

      status = sscanf (Line, "%d", &LSize);
      util_Assert (status == 1, "Error in reading LSize of generator");

      /* Consider only generators with (LSize >= i1) and (LSize <= i2) */
      if (LSize < i1)
         continue;
      if (LSize > i2)
         break;
      if (j++ % istep)
         continue;
      if (i >= fam->Ng)
         ffam_ReallocFam (fam, 2 * i);
      fam->LSize[i] = LSize;

      switch (g) {
      case LCG_a:
         fam->Resol[i] = LSize;
         status = ReadOneLCG (Line, fam, i);
         break;
      case LCGPow2_a:
         fam->Resol[i] = LSize;
         status = ReadOneLCGPow2 (Line, fam, i);
         break;
      case MRG2_a:
         fam->Resol[i] = LSize / 2;
         status = ReadOneMRG2 (Line, fam, i);
         break;
      case MRG3_a:
         fam->Resol[i] = LSize / 3;
         status = ReadOneMRG3 (Line, fam, i);
         break;
      case CombL2_a:
         fam->Resol[i] = LSize / 2;
         status = ReadOneCombL2 (Line, fam, i);
         break;
      case CombWH2_a:
         fam->Resol[i] = util_Min (LSize, 53);
         status = ReadOneCombWH2 (Line, fam, i);
         break;
      case InvImpl_a:
         fam->Resol[i] = LSize;
         status = ReadOneInvImpl (Line, fam, i);
         break;
      case InvExpl_a:
         fam->Resol[i] = LSize;
         status = ReadOneInvExpl (Line, fam, i);
         break;
      case InvMRG2_a:
         fam->Resol[i] = LSize / 2;
         status = ReadOneInvMRG2 (Line, fam, i);
         break;
      case Cubic1_a:
         fam->Resol[i] = LSize;
         status = ReadOneCubic1 (Line, fam, i);
         break;
      case CombCubic2_a:
         fam->Resol[i] = LSize;
         status = ReadOneCombCubic2 (Line, fam, i);
         break;
      case CombCubLCG_a:
         fam->Resol[i] = LSize;
         status = ReadOneCombCubLCG (Line, fam, i);
         break;
      default:
         util_Error ("fcong_Create...:   impossible case");
      }

      if (status)                 /* EOF? */
         break;
      i++;
   }
   util_Assert (i > 0, "fcong_Create...:   no generator");
   ffam_ReallocFam (fam, i);
   return fam;
}


/*=========================================================================*/

static ffam_Fam *ReadInvGen (char *filename, char *deffile, GenType g,
   int i1, int i2, int istep)
/*
 * Read parameter file for the inversive generators. If filename exists,
 * it will be read; otherwise, the default file deffile will be read.
 */
{
   FILE *f;
   char Line[MAXCAR + 1];
   ffam_Fam *fam;
   int i, j;
   int status;
   unsigned long a1, a2;

   f = ffam_OpenFile (filename, deffile);

   /* Read the family name */
   util_GetLine (f, Line, '#');
   fam = ffam_CreateFam ((i2 - i1 + istep) / istep, Line);

   /* Read the parameters */
   status = util_GetLine (f, Line, '#');
   util_Assert (status == 0, "ReadInvGen:   EOF or error");

   j = 0;
   switch (g) {
   case InvImpl2a_a:
      status = sscanf (Line, "%lu %lu ", &a1, &a2);
      util_Assert (status == 2, "ReadInvGen:   Error in reading.");
      i1 = util_Max (i1, 7);
      i2 = util_Min (i2, 31);
      for (i = i1; i <= i2; i += istep) {
	 fam->LSize[j] = i;
	 fam->Resol[j] = i + 1;
	 fam->Gen[j] = uinv_CreateInvImpl2a (i + 1, a1, a2, 1);
	 j++;
      }
      util_Assert (j > 0, "fcong_CreateInvImpl2a:    no generator!!");
      break;

   case InvImpl2b_a:
      status = sscanf (Line, "%lu %lu ", &a1, &a2);
      util_Assert (status == 2, "ReadInvGen:   Error in reading.");
      i1 = util_Max (i1, 7);
      i2 = util_Min (i2, 32);
      for (i = i1; i <= i2; i += istep) {
	 fam->LSize[j] = i;
	 fam->Resol[j] = i;
	 fam->Gen[j] = uinv_CreateInvImpl2b (i, a1, a2, 1);
	 j++;
      }
      util_Assert (j > 0, "fcong_CreateInvImpl2b:    no generator!!");
      break;

   case InvExpl2a_a:
      status = sscanf (Line, "%lu", &a1);
      util_Assert (status == 1, "ReadInvGen:   Error in reading.");
      i1 = util_Max (i1, 7);
      i2 = util_Min (i2, 32);
      for (i = i1; i <= i2; i += istep) {
	 fam->LSize[j] = i;
	 fam->Resol[j] = i;
	 fam->Gen[j] = uinv_CreateInvExpl2a (i, (long) a1, 1);
	 j++;
      }
      util_Assert (j > 0, "fcong_CreateInvExpl2a:    no generator!!");
      break;

   case InvExpl2b_a:
      status = sscanf (Line, "%lu", &a1);
      util_Assert (status == 1, "ReadInvGen:   Error in reading.");
      i1 = util_Max (i1, 7);
      i2 = util_Min (i2, 32);
      for (i = i1; i <= i2; i += istep) {
	 fam->LSize[j] = i;
	 fam->Resol[j] = i;
	 fam->Gen[j] = uinv_CreateInvExpl2b (i, (long) a1, 1);
	 j++;
      }
      util_Assert (j > 0, "fcong_CreateInvExpl2b:    no generator!!");
      break;

   default:
      util_Error ("ReadInvGen:   impossible case");
   }

   ffam_ReallocFam (fam, j);
   return fam;
}


/*=========================================================================*/

ffam_Fam *fcong_CreateLCG (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "LCGGood.par", LCG_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteLCG (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++) {
      if (fam->LSize[i] <= 31)
         ulcg_DeleteGen (fam->Gen[i]);
#ifdef USE_GMP
      else
         ulcg_DeleteBigLCG (fam->Gen[i]);
#endif
   }
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateLCGPow2 (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "LCGPow2.par", LCGPow2_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteLCGPow2 (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++) {
      if (fam->LSize[i] <= 31)
         ulcg_DeleteGen (fam->Gen[i]);
#ifdef USE_GMP
      else
         ulcg_DeleteBigPow2LCG (fam->Gen[i]);
#endif
   }
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateMRG2 (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "MRG2.par", MRG2_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteMRG2 (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++)
      umrg_DeleteMRGFloat (fam->Gen[i]);
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateMRG3 (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "MRG3.par", MRG3_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteMRG3 (ffam_Fam * fam)
{
   fcong_DeleteMRG2 (fam);
}


/*=========================================================================*/


ffam_Fam *fcong_CreateCombL2 (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "CombL2.par", CombL2_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteCombL2 (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++)
      ulcg_DeleteGen (fam->Gen[i]);
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateCombWH2 (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "CombWH2.par", CombWH2_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteCombWH2 (ffam_Fam * fam)
{
   fcong_DeleteCombL2 (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateInvImpl (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "InvImpl.par", InvImpl_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteInvImpl (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++)
      uinv_DeleteGen (fam->Gen[i]);
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateInvImpl2a (char *filename, int i1, int i2, int istep)
{
   return ReadInvGen (filename, "InvImpl2a.par", InvImpl2a_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteInvImpl2a (ffam_Fam * fam)
{
   fcong_DeleteInvImpl (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateInvImpl2b (char *filename, int i1, int i2, int istep)
{
   return ReadInvGen (filename, "InvImpl2b.par", InvImpl2b_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteInvImpl2b (ffam_Fam * fam)
{
   fcong_DeleteInvImpl (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateInvExpl (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "InvExpl.par", InvExpl_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteInvExpl (ffam_Fam * fam)
{
   fcong_DeleteInvImpl (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateInvExpl2a (char *filename, int i1, int i2, int istep)
{
   return ReadInvGen (filename, "InvExpl2a.par", InvExpl2a_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteInvExpl2a (ffam_Fam * fam)
{
   fcong_DeleteInvImpl (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateInvExpl2b (char *filename, int i1, int i2, int istep)
{
   return ReadInvGen (filename, "InvExpl2b.par", InvExpl2b_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteInvExpl2b (ffam_Fam * fam)
{
   fcong_DeleteInvImpl (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateInvMRG2 (char *filename, int i1, int i2, int istep)
{
   ffam_Fam *fam;
   size_t len;
   fam = ReadAllGen (filename, "MRG2.par", InvMRG2_a, i1, i2, istep);
   len = strlen ("InvMRG2");
   fam->name = util_Realloc (fam->name, (1 + len) * sizeof (char));
   strncpy (fam->name, "InvMRG2", (size_t) len);
   return fam;
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteInvMRG2 (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++)
      uinv_DeleteInvMRGFloat (fam->Gen[i]);
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateCubic1 (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "Cubic1.par", Cubic1_a, i1, i2, istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteCubic1 (ffam_Fam * fam)
{
   int i;
   for (i = 0; i < fam->Ng; i++)
      ucubic_DeleteGen (fam->Gen[i]);
   ffam_DeleteFam (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateCombCubic2 (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "CombCubic2.par", CombCubic2_a, i1, i2,
      istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteCombCubic2 (ffam_Fam * fam)
{
   fcong_DeleteCubic1 (fam);
}


/*=========================================================================*/

ffam_Fam *fcong_CreateCombCubLCG (char *filename, int i1, int i2, int istep)
{
   return ReadAllGen (filename, "CombCubLCG.par", CombCubLCG_a, i1, i2,
      istep);
}

/*-------------------------------------------------------------------------*/

void fcong_DeleteCombCubLCG (ffam_Fam * fam)
{
   unif01_Comb2_Param *param;
   int i;

   for (i = 0; i < fam->Ng; i++) {
      param = fam->Gen[i]->param;
      ulcg_DeleteGen (param->gen1);
      ucubic_DeleteGen (param->gen2);
      unif01_DeleteCombGen (fam->Gen[i]);
   }
   ffam_DeleteFam (fam);
}


/*=========================================================================*/
