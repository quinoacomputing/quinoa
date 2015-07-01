/*************************************************************************\
 *
 * Package:        TestU01
 * File:           usoft.c
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
#include "num.h"
#include "addstr.h"

#include "usoft.h"
#include "ulcg.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



#ifdef HAVE_MATHEMATICA
#include <pthread.h>
#include "mathlink.h"
#endif



/*============================= Constants ================================*/

#define DeuxExp31m1  2147483647             /* 2^31 - 1 */
#define MASK52       4503599627370495ULL    /* 2^52 - 1 */
#define MASK48       281474976710655ULL     /* 2^48 - 1 */
#define MASK32       4294967295             /* 2^32 - 1 */
#define MASK5        31                     /* 2^5  - 1 */

#define LEN    200                          /* Max length of strings */



/*================================ Types ================================*/


typedef struct {
   double Norm;
} SPlus_param;

typedef struct {
   unsigned long S1, S2;
} SPlus_state;

/*------------------------------------*/
#ifdef USE_LONGLONG

typedef struct {
   double Norm;
   int Flag;
} Java48_param;

typedef struct {
   ulonglong S;
} Java48_state;


typedef struct {
   double Z[32];
   double b;
   unsigned int i, j;
} MATLAB5_state;

#endif
/*------------------------------------*/

typedef struct {
   double U;
} Excel97_state;

/*------------------------------------*/

typedef struct {
   unsigned long S;
} VisualBasic_state;






/*============================== Functions ==============================*/

static double SPlus_U01 (void *vpar, void *vsta)
{
   SPlus_state *state = vsta;
   SPlus_param *param = vpar;
   unsigned int z;

   do {
      state->S1 *= 69069;
      state->S2 ^= state->S2 >> 15;
      state->S2 ^= state->S2 << 17;
#ifndef IS_ULONG32
      state->S1 &= unif01_MASK32;
      state->S2 &= unif01_MASK32;
#endif
      z = (state->S1 ^ state->S2) >> 1;
   } while (0 == z); 
   return param->Norm * z;      
}

/*-----------------------------------------------------------------------*/

static unsigned long SPlus_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (SPlus_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrSPlus (void *vsta)
{
   SPlus_state *state = vsta;
   printf (" S1 = %1lu,    S2 = %1lu\n", state->S1, state->S2);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * usoft_CreateSPlus (long s1, long s2)
{
   unif01_Gen *gen;
   SPlus_state *state;
   SPlus_param *param;
   size_t leng;
   char name[LEN + 1];

   util_Assert (s1 > 0, "usoft_CreateSPlus:   must have s1 > 0");
   util_Assert (s1 < DeuxExp31m1,
      "usoft_CreateSPlus:   must have s1 < 2^31 - 1");
   util_Assert (s2 > 0, "usoft_CreateSPlus:   must have s2 > 0");
   util_Assert (s2 < DeuxExp31m1,
      "usoft_CreateSPlus:   must have s2 < 2^31 - 1");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (SPlus_param));
   state = util_Malloc (sizeof (SPlus_state));

   strcpy (name, "usoft_CreateSPlus:");
   addstr_Long (name, "   s1 = ", s1);
   addstr_Long (name, ",   s2 = ", s2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = s1;
   state->S2 = s2;
   param->Norm = 1.0 / num_TwoExp[31];

   gen->GetBits = &SPlus_Bits;
   gen->GetU01  = &SPlus_U01;
   gen->Write   = &WrSPlus;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/
/* The 5 generators random of Unix-Linux */

#ifdef HAVE_RANDOM

static int coUnix = 0;                 /* Counter for Unix random RNG */

static unsigned long state1[] = {
   3UL,
   0x9a319039UL, 0x32d9c024UL, 0x9b663182UL, 0x5da1f342UL,
   0x7449e56bUL, 0xbeb1dbb0UL, 0xab5c5918UL, 0x946554fdUL,
   0x8c2e680fUL, 0xeb3d799fUL, 0xb11ee0b7UL, 0x2d436b86UL,
   0xda672e2aUL, 0x1588ca88UL, 0xe369735dUL, 0x904f35f7UL,
   0xd7158fd6UL, 0x6fa6f051UL, 0x616e6b96UL, 0xac94efdcUL,
   0xde3b81e0UL, 0xdf0a6fb5UL, 0xf103bc02UL, 0x48f340fbUL,
   0x36413f93UL, 0xc622c298UL, 0xf5a42ab8UL, 0x8a88d77bUL,
   0xf5ad9d0eUL, 0x8999220bUL, 0x27fb47b9UL,
   12345UL,
   0x9a319039UL, 0x32d9c024UL, 0x9b663182UL, 0x5da1f342UL,
   0x7449e56bUL, 0xbeb1dbb0UL, 0xab5c5918UL, 0x946554fdUL,
   0x8c2e680fUL, 0xeb3d799fUL, 0xb11ee0b7UL, 0x2d436b86UL,
   0xda672e2aUL, 0x1588ca88UL, 0xe369735dUL, 0x904f35f7UL,
   0xd7158fd6UL, 0x6fa6f051UL, 0x616e6b96UL, 0xac94efdcUL,
   0xde3b81e0UL, 0xdf0a6fb5UL, 0xf103bc02UL, 0x48f340fbUL,
   0x36413f93UL, 0xc622c298UL, 0xf5a42ab8UL, 0x8a88d77bUL,
   0xf5ad9d0eUL, 0x8999220bUL, 0x27fb47b9UL
};

/*------------------------------------------------------------------------*/

static void WrUnixRandom (void *junk)
{
}

/*------------------------------------------------------------------------*/

static double UnixRandom_U01 (void *junk1, void *junk2)
{
   return random () / (1.0 + RAND_MAX);
}

/*------------------------------------------------------------------------*/

static unsigned long UnixRandom_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (UnixRandom_U01 (vpar, vsta) * unif01_NORM32);
}

/*------------------------------------------------------------------------*/

unif01_Gen * usoft_CreateUnixRandom (unsigned int s)
{
   unif01_Gen *gen;
   unsigned seed;
   size_t leng;
   char name[LEN + 1];

   util_Assert (coUnix == 0,
      "usoft_CreateUnixRandom:   only 1 generator at a time can be in use");
   coUnix++;

   switch (s) {
   case 8:
   case 32:
   case 64:
   case 128:
   case 256:
      break;
   default:
      util_Error ("\nusoft_CreateUnixRandom:   "
                  "s must be in {8, 32, 64, 128, 256}\n\n");
   }
   gen = util_Malloc (sizeof (unif01_Gen));

   seed = 12345;
   initstate (seed, (char *) state1, (size_t) s);   /* Defined in stdlib */
   setstate ((char *) state1);

   strcpy (name, "usoft_CreateUnixRandom:");
   addstr_Uint (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &UnixRandom_Bits;
   gen->GetU01 = &UnixRandom_U01;
   gen->Write = &WrUnixRandom;
   gen->param = NULL;
   gen->state = NULL;
   return gen;
}

/*-----------------------------------------------------------------------*/

void usoft_DeleteUnixRandom (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   coUnix--;
}

#endif                            /* HAVE_RANDOM */



/*=========================================================================*/
#ifdef USE_LONGLONG

static double Java48_U01 (void *vpar, void *vsta)
{
   Java48_param *param = vpar;
   Java48_state *state = vsta;
   ulonglong temp;

   state->S = (25214903917ULL * state->S + 11) & MASK48;
   /* Keep only the 26 most significant bits of S; shift 5.  They will */
   /* make bits [27..52] of the next generated number, counting from 0. */
   temp = (state->S & 281474972516352ULL) << 5;

   state->S = (25214903917ULL * state->S + 11) & MASK48;
   /* Keep only the 27 most significant bits of S; shift 21. They will */
   /* make bits [0..26] of the next generated number */
   return (temp + (state->S >> 21)) * param->Norm;
}

/*-----------------------------------------------------------------------*/

static unsigned long Java48_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (Java48_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrJava48 (void *vsta)
{
   Java48_state *state = vsta;
   printf (" S =  %" PRIuLEAST64 "\n", state->S);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * usoft_CreateJava48 (ulonglong s, int jflag)
{
   unif01_Gen *gen;
   Java48_param *param;
   Java48_state *state;
   size_t leng;
   char name[LEN + 1];

   if (s > MASK48)
      util_Error ("usoft_CreateJava48:   s >= 281474976710656");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Java48_param));
   state = util_Malloc (sizeof (Java48_state));

   strcpy (name, "usoft_CreateJava48:");
   addstr_ULONG (name, "   s = ", s);
   addstr_Long (name, ",   jflag = ", (long) jflag);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   /* This bit XOR is used to get the same numbers as those generated by the 
      SUN Random.nextDouble */
   if (jflag)
      state->S = s ^ 0x5DEECE66DULL;
   else
      state->S = s;

   param->Norm = 1.0 / num_TwoExp[53];
   gen->GetBits = &Java48_Bits;
   gen->GetU01  = &Java48_U01;
   gen->Write   = &WrJava48;
   gen->param   = param;
   gen->state   = state;
   return gen;
}

#endif                            /* USE_LONGLONG */


/*=========================================================================*/

static double Excel97_U01 (void *junk, void *vsta)
{
   Excel97_state *state = vsta;
   state->U = 9821.0 * state->U + 0.211327;
   state->U -= (int) state->U;
   return state->U;
}

/*-----------------------------------------------------------------------*/

static unsigned long Excel97_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (Excel97_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrExcel97 (void *vsta)
{
   Excel97_state *state = vsta;
   printf (" R = %20.16f\n", state->U);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * usoft_CreateExcel97 (double r)
{
   unif01_Gen *gen;
   Excel97_state *state;
   size_t leng;
   char name[LEN + 1];

   if (r < 0.0 || r >= 1.0)
      util_Error ("usoft_CreateExcel97:   r must be in [0, 1)");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (Excel97_state));

   strcpy (name, "usoft_CreateExcel97:");
   addstr_Double (name, "   r = ", r);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->U = r;
   gen->GetBits = &Excel97_Bits;
   gen->GetU01  = &Excel97_U01;
   gen->Write   = &WrExcel97;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}


/*=========================================================================*/

unif01_Gen * usoft_CreateExcel2003 (int x0, int y0, int z0)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   gen = ulcg_CreateCombWH3 (30323, 30307, 30269, 170, 172, 171,
                             0, 0, 0, x0, y0, z0);
   strcpy (name, "usoft_CreateExcel2003:");
   addstr_Uint (name, "   x0 = ", x0);
   addstr_Uint (name, ",   y0 = ", y0);
   addstr_Uint (name, ",   z0 = ", z0);
   leng = strlen (name);
   gen->name = util_Free (gen->name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);
   return gen;
}


/*=========================================================================*/

static unsigned long VisualBasic_Bits (void *junk, void *vsta)
{
   VisualBasic_state *state = vsta;
   state->S = (16598013 * state->S + 12820163) & 0xffffff;
   return (state->S << 8);
}

/*-----------------------------------------------------------------------*/

static double VisualBasic_U01 (void *vpar, void *vsta)
{
   return VisualBasic_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrVisualBasic (void *vsta)
{
   VisualBasic_state *state = vsta;
   printf (" S = %1lu\n", state->S);
}

/*-----------------------------------------------------------------------*/

unif01_Gen * usoft_CreateVisualBasic (unsigned long s)
{
   unif01_Gen *gen;
   VisualBasic_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (VisualBasic_state));

   strcpy (name, "usoft_CreateVisualBasic:");
   addstr_Ulong (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S = s;
   gen->GetBits = &VisualBasic_Bits;
   gen->GetU01  = &VisualBasic_U01;
   gen->Write   = &WrVisualBasic;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}



/*=========================================================================*/
#if defined(USE_GMP) && defined(USE_LONGLONG)

unif01_Gen * usoft_CreateMaple_9 (longlong s)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   char str[20];

   sprintf (str, "%1" PRIdLEAST64, s);
   gen = ulcg_CreateBigLCG ("999999999989", "427419669081", "0", str);
   gen->name = util_Free (gen->name);
   strcpy (name, "usoft_CreateMaple_9:");
   addstr_LONG (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);
   return gen;
}


void usoft_DeleteMaple_9 (unif01_Gen *gen)
{
   ulcg_DeleteBigLCG (gen);
}

#endif


/*=========================================================================*/
#ifdef USE_LONGLONG
#define IMAX 32

static double MATLAB5_U01 (void *junk, void *vsta)
{
   static const double ULP = 1.0 / 9007199254740992.0;    /* 1 / 2^53 */
   MATLAB5_state *state = vsta;
   ulonglong t, mask = state->j;
   int n;
   double x;

   x = state->Z[(state->i + 20) & MASK5] - state->Z[(state->i + 5) & MASK5]
          - state->b;
   if (x < 0.0) {
      x += 1.0;
      state->b = ULP;
   } else
      state->b = 0.0;
   state->Z[state->i] = x;
   state->i = (state->i + 1) & MASK5;

   state->j ^= (state->j<<13);
   state->j ^= (state->j>>17);
   state->j ^= (state->j<<5);
   mask ^= ((((ulonglong) state->j) << 32) & MASK52);
   x = frexp (x, &n);
   t = ldexp (x, 53);
   t ^= mask;
   x = ldexp ((double) t, n - 53);
   return x;
}


/*-----------------------------------------------------------------------*/

static unsigned long MATLAB5_Bits (void *vpar, void *vsta)
{
   return unif01_NORM32 * MATLAB5_U01 (vpar, vsta);
}

/*-----------------------------------------------------------------------*/

static void WrMATLAB5 (void *vsta)
{
   unsigned int j;
   MATLAB5_state *state = vsta;
   printf (" i = %1u,", state->i);
   printf ("   j = %1u,", state->j);
   printf ("   b = %d,\n Z = ", state->b > 0.0 ? 1 : 0);
   if (unif01_WrLongStateFlag) {
      printf (" {\n ");
      for (j = 0; j < IMAX; j++) {
         printf ("   %.16f", state->Z[j]);
         if (j < IMAX - 1)
            printf (",");
         if ((j % 3) == 2)
            printf ("\n ");
      };
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen * usoft_CreateMATLAB (int i, unsigned int j0, int b, double Z[])
{
   unif01_Gen *gen;
   MATLAB5_state *state;
   size_t leng;
   char name[LEN + 1];
   int r;

   strcpy (name, "usoft_CreateMATLAB:");
   addstr_Int (name, "   i = ", i);
   if (i >= 0) {
      addstr_Uint (name, ",   j = ", j0);
      addstr_Int (name, ",   b = ", b);
      util_Assert (Z != NULL, "usoft_CreateMATLAB:   Z is NULL");
      addstr_ArrayDouble (name, ",   Z = ", IMAX, Z);
   }

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MATLAB5_state));

   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (i < 0) {
      unsigned int s, j;
      double x;
      state->b = 0;
      state->i = 0;
      j = state->j = 1U << 31;
      /* RANDSETUP in randtx: Generate n reals in [0,1] bit by bit. */
      for (r = 0; r < IMAX; r++) {
         x = 0;
         for (s = 0; s < 53; s++) {
            j ^= (j<<13);  j ^= (j>>17);  j ^= (j<<5);
            x = 2*x + ((j >> 19) & 1);
         }
         state->Z[r] = ldexp (x, -53);
      }

   } else {
      double junk;
      for (r = 0; r < IMAX; r++) {
         util_Assert (state->Z[r] >= 0.0,
                      "usoft_CreateMATLAB:   negative Z[r]");
         state->Z[r] = modf (Z[r], &junk);
      }
      state->b = b > 0 ? 1.0/num_TwoExp[53] : 0.0;
      state->i = i & MASK5;
      state->j = j0 > 0 ? j0 : (1U << 31);
   }

   gen->param = NULL;
   gen->state = state;
   gen->GetBits = &MATLAB5_Bits;
   gen->GetU01  = &MATLAB5_U01;
   gen->Write   = &WrMATLAB5;
   return gen;
}


/*-----------------------------------------------------------------------*/

void usoft_DeleteMATLAB (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}
#undef IMAX
#endif


/*=========================================================================*/
#ifdef HAVE_MATHEMATICA

/* The Mathematica RNG for the Real */
static int coMath = 0;                /* Counter */
static double *Math_A;
static long Math_n = 262144;
static MLINK lp;
static MLEnvironment env;

#if 0
static double Mathematica00 (void)
/* 
 * Very slow function: generate 1 number at a time
 */
{
   double x;
   MLPutFunction (lp, "Random", 0);
   MLEndPacket (lp);
   while (MLNextPacket (lp) != RETURNPKT)
      MLNewPacket (lp);
   MLGetDouble (lp, &x);
   return x;
}
#endif

static void GetArray (void)
{
   /* Ask Mathematica to generate Math_n random numbers at a time and put
      them in array Math_A */
   static int flag = 0;

   /* free memory taken up by previous array Math_A */
   if (flag)
      MLDisownRealList (lp, Math_A, Math_n);
   flag = 1;

   /* send the command Table[Random[], {Math_n}] to the Mathematica kernel */
   MLPutFunction (lp, "Table", 2);
   MLPutFunction (lp, "Random", 0);
   MLPutFunction (lp, "List", 1);
   MLPutInteger (lp, Math_n);
   MLEndPacket (lp);
   while (MLNextPacket (lp) != RETURNPKT)
      MLNewPacket (lp);

   /* Get the Math_n numbers in Math_A */
   MLGetRealList (lp, &Math_A, &Math_n);
}

/*------------------------------------------------------------------------*/

static double Mathematica_U01 (void *junk1, void *junk2)
{
   /* return one random number; when at top of array, get another bunch */
   static long i = 0;
   if (i == Math_n) {
      GetArray ();
      i = 0;
   }
   return Math_A[i++];
}

/*------------------------------------------------------------------------*/

static unsigned long Mathematica_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Mathematica_U01 (vpar, vsta));
}

/*------------------------------------------------------------------------*/

static void WrMathematica (void *junk)
{
   /* The state is expressed as a very big integer (many decimal digits);
      do not write it */
}

/*------------------------------------------------------------------------*/

unif01_Gen *usoft_CreateMathematicaReal (int argc, char *argv[], long s)
{
   long errno;
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   if (argc == 1) {
      util_Error ("usoft_CreateMathematicaReal:\n\n   Usage:  <prog> -linkname 'math -mathlink' -linklaunch");
   }
   util_Assert (coMath == 0,
      "usoft_CreateMathematicaReal:   only 1 generator at a time can be used");
   coMath++;

   env = MLInitialize (NULL);
   if (env == NULL) {
      util_Assert (1, "usoft_CreateMathematicaReal:   MathLink MLInitialize returns NULL");
      return NULL;
   }
   lp = MLOpenArgv (env, argv, argv + argc, &errno);
   if (lp == NULL) {
      util_Assert (1, "usoft_CreateMathematicaReal:   MathLink MLOpenArgv returns NULL");
      return NULL;
   }

   MLPutFunction (lp, "SeedRandom", 1);
   MLPutInteger (lp, s);
   MLEndPacket (lp);
   while (MLNextPacket (lp) != RETURNPKT)
      MLNewPacket (lp);
   MLNewPacket (lp);

   gen = util_Malloc (sizeof (unif01_Gen));

   strcpy (name, "usoft_CreateMathematicaReal:");
   addstr_Long (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   GetArray ();
   gen->GetBits = &Mathematica_Bits;
   gen->GetU01  = &Mathematica_U01;
   gen->Write   = &WrMathematica;
   gen->param   = NULL;
   gen->state   = NULL;
   return gen;
}

/*------------------------------------------------------------------------*/

void usoft_DeleteMathematicaReal (unif01_Gen *gen)
{
   MLDisownRealList (lp, Math_A, Math_n);
   MLClose (lp);
   MLDeinitialize (env);
   if (NULL == gen)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   coMath--;
}


/*=========================================================================*/

/* The Mathematica RNG for the Integer */
static int coMathI = 0;                /* Counter */
static int *Math_AI;
static long Math_nI = 262144;
static MLINK lpI;
static MLEnvironment envI;


static void GetArrayI (void)
{
/* Ask Mathematica to generate Math_nI random numbers at a time and put
   them in array Math_AI */
   static int flag = 0;

   /* free memory taken up by previous array Math_AI */
   if (flag)
      MLDisownIntegerList (lpI, Math_AI, Math_nI);
   flag = 1;

   /* send the command Table[Random[Integer, {0, 2^30 - 1}], {Math_nI}] to
      the Mathematica kernel */
   MLPutFunction (lpI, "Table", 2);
    MLPutFunction (lpI, "Random", 2);
      MLPutSymbol (lpI, "Integer");
      MLPutFunction (lpI, "List", 2);
         MLPutInteger (lpI, 0);
         MLPutInteger (lpI, 1073741823);  /* 2^30 - 1; CA gives 30 bits max */
    MLPutFunction (lpI, "List", 1);
      MLPutInteger (lpI, Math_nI);
   MLEndPacket (lpI);
   while (MLNextPacket (lpI) != RETURNPKT)
      MLNewPacket (lpI);

   /* Get the Math_nI numbers in Math_AI */
   MLGetIntegerList (lpI, &Math_AI, &Math_nI);
}


/*------------------------------------------------------------------------*/

static unsigned long MathematicaI_Bits (void *junk1, void *junk2)
{
   /* return one random number; when at top of array, get another bunch */
   static long i = 0;
   if (i == Math_nI) {
      GetArrayI ();
      i = 0;
   }
   /* Math_AI: Random numbers in [0, 2^30 - 1] */
   return ((unsigned int) Math_AI[i++] << 2);
}


/*------------------------------------------------------------------------*/

static double MathematicaI_U01 (void *vpar, void *vsta)
{
   return unif01_INV32 * MathematicaI_Bits (vpar, vsta);
}


/*------------------------------------------------------------------------*/

unif01_Gen *usoft_CreateMathematicaInteger (int argc, char *argv[], long s)
{
   long errno;
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1] = {0};

   if (argc == 1) {
      util_Error ("usoft_CreateMathematicaInteger:\n\n   Usage:  <prog> -linkname 'math -mathlink' -linklaunch");
   }
   util_Assert (coMathI == 0,
   "usoft_CreateMathematicaInteger:   only 1 generator at a time can be used");
   coMathI++;

   envI = MLInitialize (NULL);
   if (envI == NULL) {
      util_Assert (1, "usoft_CreateMathematicaInteger:   MathLink MLInitialize returns NULL");
      return NULL;
   }
   lpI = MLOpenArgv (envI, argv, argv + argc, &errno);
   if (lpI == NULL) {
      util_Assert (1, "usoft_CreateMathematicaInteger:   MathLink MLOpenArgv returns NULL");
      return NULL;
   }

   MLPutFunction (lpI, "SeedRandom", 1);
   MLPutInteger (lpI, s);
   MLEndPacket (lpI);
   while (MLNextPacket (lpI) != RETURNPKT)
      MLNewPacket (lpI);
   MLNewPacket (lpI);

   gen = util_Malloc (sizeof (unif01_Gen));

   strcpy (name, "usoft_CreateMathematicaInteger:");
   addstr_Long (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   GetArrayI ();
   gen->GetBits = &MathematicaI_Bits;
   gen->GetU01  = &MathematicaI_U01;
   gen->Write   = &WrMathematica;
   gen->param   = NULL;
   gen->state   = NULL;
   return gen;
}

/*------------------------------------------------------------------------*/

void usoft_DeleteMathematicaInteger (unif01_Gen *gen)
{
   MLDisownIntegerList (lpI, Math_AI, Math_nI);
   MLClose (lpI);
   MLDeinitialize (envI);
   if (NULL == gen)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   coMathI--;
}

#endif


/*=========================================================================*/
/* a simple example of using the Kernel of Mathematica:
   request to Mathematica to compute the sum of 2 integers interactively
   on the terminal from a C program.

 Compile with
       gcc umath.c -I$MYMLPATH -L$MYMLPATH -lML -lm
 Run with
       a.out -linkname 'math -mathlink' -linkmode launch

*/

#if 0
int main (int argc, char *argv[])
{

   int i, j, sum;
   MLINK lp;
   int pkt;
   MLEnvironment env;

   printf ("Enter two integers: \t\n");
   scanf ("%d %d", &i, &j);
   env = MLInitialize (NULL);
   if (env == NULL)
      return 1;

   lp = MLOpen (argc, argv);
   if (lp == NULL)
      return 1;

   /* Send Plus[i, j] */
   MLPutFunction (lp, "Plus", 2);
   MLPutInteger (lp, i);
   MLPutInteger (lp, j);

   /* skip any packets before the first ReturnPacket */
   while (MLNextPacket (lp) != RETURNPKT)
      MLNewPacket (lp);

   /* inside the ReturnPacket we expect an integer */
   MLGetInteger (lp, &sum);
   printf ("sum = %d\n\n", sum);

   MLClose (lp);
   MLDeinitialize (env);
   return 0;
}
#endif

/*=========================================================================*/

void usoft_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
