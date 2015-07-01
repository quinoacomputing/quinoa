/*************************************************************************\
 *
 * Package:        TestU01
 * File:           unif01.c
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
#include "chrono.h"
#include "swrite.h"
#include "unif01.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define LEN0 500                   /* Length of strings */
#define LEN1 100                   /* Length of strings */

#define MASK32 0xffffffffUL        /* 2^32 - 1 */




/*------------------------- extern variables ------------------------------*/

lebool unif01_WrLongStateFlag = FALSE;





/* ========================== functions ================================== */

void unif01_WriteNameGen (unif01_Gen *gen)
{
   if (gen->name)
      printf ("%s\n\n", gen->name);
}

void unif01_WriteState (unif01_Gen *gen)
{
   printf ("\nGenerator state:\n");
   gen->Write (gen->state);
   printf ("\n");
}

void unif01_WrLongStateDef (void)
{
   printf ("  Not shown here ... takes too much space\n");
}


/**************************************************************************/

double unif01_StripD (unif01_Gen *gen, int r)
{
   if (r == 0) {
      return (gen->GetU01) (gen->param, gen->state);
   } else {
      double u = num_TwoExp[r] * (gen->GetU01) (gen->param, gen->state);
      return (u - (long) u);
   }
}

long unif01_StripL (unif01_Gen *gen, int r, long d)
{
   if (r == 0)
      return (long) (d * gen->GetU01 (gen->param, gen->state));
   else {
      double u = num_TwoExp[r] * (gen->GetU01) (gen->param, gen->state);
      return (long) (d * (u - (long) u));
   }  
}

unsigned long unif01_StripB (unif01_Gen *gen, int r, int s)
{
   if (r == 0) {
      return gen->GetBits (gen->param, gen->state) >> (32 - s);
   } else {
      unsigned long u = gen->GetBits (gen->param, gen->state);
      return ((u << r) & MASK32) >> (32 - s);
   }
}


/*************************************************************************/

/* Dummy generator, always return 0.  */

static double DummyGen_U01 (void *junk1, void *junk2)
{
   return 0.0;
}

static unsigned long DummyGen_Bits (void *junk1, void *junk2)
{
   return 0;
}

static void WrDummyGen (void *junk)
{
   printf ("   Empty Generator (no state)\n");
}

unif01_Gen * unif01_CreateDummyGen (void)
{
   unif01_Gen *gen;
   size_t len;

   gen = util_Malloc (sizeof (unif01_Gen));
   len = strlen ("Dummy generator that always returns 0");
   gen->name    = util_Calloc (len + 1, sizeof (char));
   strncpy (gen->name, "Dummy generator that always returns 0", len);
   gen->param   = NULL;
   gen->state   = NULL;
   gen->Write   = &WrDummyGen;
   gen->GetBits = &DummyGen_Bits;
   gen->GetU01  = &DummyGen_U01;
   return gen;
}

void unif01_DeleteDummyGen (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
}

void unif01_DeleteGen (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/**************************************************************************/
/*
 * The original generator is gen0. The position of the bit at which the 
 * increased precision is applied is s, counting from the most significant 
 * bit;   v = 1 / 2^s.
 */
typedef struct {
   unif01_Gen *gen0;
   double v;
   int s;
} DoubleGen_param;


static double DoubleGen_U01 (void *vpar, void *junk)
{
   double U;
   DoubleGen_param *paramD = vpar;
   unif01_Gen *gen = paramD->gen0;

   U = gen->GetU01 (gen->param, gen->state);
   U += paramD->v * gen->GetU01 (gen->param, gen->state);
   if (U < 1.0)
      return U;
   else
      return U - 1.0;
}

static unsigned long DoubleGen_Bits (void *vpar, void *junk)
{
   return (unsigned long) (unif01_NORM32 * DoubleGen_U01 (vpar, junk));
}


unif01_Gen * unif01_CreateDoubleGen2 (unif01_Gen *gen, double v)
{
   unif01_Gen *genD;
   DoubleGen_param *paramD;
   char *name;
   char str[20];
   size_t len, len2, len3;

   util_Assert (v > 0.0, "unif01_CreateDoubleGen2:   h <= 0");
   util_Assert (v < 1.0, "unif01_CreateDoubleGen2:   h >= 1");
   genD = util_Malloc (sizeof (unif01_Gen));
   paramD = util_Malloc (sizeof (DoubleGen_param));
   paramD->s = -num_Log2(v);
   paramD->v = v;
   paramD->gen0 = gen;

   len = strlen (gen->name);
   len2 = strlen ("\nunif01_CreateDoubleGen2 with h = ");
   len += len2;
   sprintf (str, "%-g", v);
   len3 = strlen (str);
   len += len3;
   name = util_Calloc (len + 1, sizeof (char));
   strncpy (name, gen->name, len);
   strncat (name, "\nunif01_CreateDoubleGen2 with h = ", len2);
   strncat (name, str, len3);

   /* The state of the double generator is simply the state of the original
      generator */
   genD->name    = name;
   genD->param   = paramD;
   genD->state   = gen->state;
   genD->Write   = gen->Write;
   genD->GetBits = &DoubleGen_Bits;
   genD->GetU01  = &DoubleGen_U01;
   return genD;
}

unif01_Gen * unif01_CreateDoubleGen (unif01_Gen *gen, int s)
{
   unif01_Gen *genD;
   DoubleGen_param *paramD;
   char *name;
   char str[8];
   size_t len, len2, len3;

   util_Assert (s > 0, "unif01_CreateDoubleGen:   s <= 0");
   genD = unif01_CreateDoubleGen2 (gen, 1.0 / num_TwoExp[s]);
   paramD = genD->param;
   paramD->s = s;

   len = strlen (gen->name);
   len2 = strlen ("\nunif01_CreateDoubleGen with s = ");
   len += len2;
   sprintf (str, "%-d", paramD->s);
   len3 = strlen (str);
   len += len3;
   name = util_Calloc (len + 1, sizeof (char));
   strncpy (name, gen->name, len);
   strncat (name, "\nunif01_CreateDoubleGen with s = ", len2);
   strncat (name, str, len3);
   genD->name    = name;
   return genD;
}

void unif01_DeleteDoubleGen (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/**************************************************************************/

typedef struct {
   unif01_Gen *gen0;                /* Original generator */
   long *ILac;                      /* Table of lacunary indices */
   int k;                           /* Size of ILac */
   int cur;                         /* Current index in the table */
   long n;
} LacGen_param;

static double LacGen_U01 (void *vpar, void *junk)
{
   LacGen_param *paramL = vpar;
   unif01_Gen *gen = paramL->gen0;
   int cur = paramL->cur;
   long *ILac = paramL->ILac;
   long j;
#if 1
   if (cur > 0) {
     for (j = 2; j <= ILac[cur] - ILac[cur - 1]; j++)
        gen->GetU01 (gen->param, gen->state);
   } else {
     for (j = 0; j < ILac[0]; j++)
        gen->GetU01 (gen->param, gen->state);
   }
   cur++;
   if (cur >= paramL->k)
      cur = 0;
   paramL->cur = cur;

#else
   /* For debugging: write the lacunary indices of the random numbers outputted */
   if (cur > 0) {
      for (j = 2; j <= ILac[cur] - ILac[cur - 1]; j++) {
         gen->GetU01 (gen->param, gen->state);
         paramL->n++;
      }
   } else {
      for (j = 0; j < ILac[0]; j++) {
         gen->GetU01 (gen->param, gen->state);
         paramL->n++;
      }
   }
   cur++;
   if (cur >= paramL->k)
      cur = 0;
   paramL->cur = cur;
   printf ("Lac = %ld\n", paramL->n);
   paramL->n++;
#endif

   return gen->GetU01 (gen->param, gen->state);
}

static unsigned long LacGen_Bits (void *vpar, void *junk)
{
   LacGen_param *paramL = vpar;
   unif01_Gen *gen = paramL->gen0;
   int cur = paramL->cur;
   long *ILac = paramL->ILac;
   long j;

   if (cur > 0) {
     for (j = 2; j <= ILac[cur] - ILac[cur - 1]; j++)
        gen->GetBits (gen->param, gen->state);
   } else {
     for (j = 0; j < ILac[0]; j++)
        gen->GetBits (gen->param, gen->state);
   }
   cur++;
   if (cur >= paramL->k)
      cur = 0;
   paramL->cur = cur;
   return gen->GetBits (gen->param, gen->state);
}

unif01_Gen * unif01_CreateLacGen (unif01_Gen *gen, int k, long I[])
{
   unif01_Gen *genL;
   LacGen_param *paramL;
   char name[LEN0 + 1] = "";
   char str[16];
   size_t len, len2;
   int j;

   genL = util_Malloc (sizeof (unif01_Gen));
   paramL = util_Malloc (sizeof (LacGen_param));
   paramL->gen0 = gen;
   paramL->k = k;
   paramL->cur = 0;
   paramL->n = 0;
   paramL->ILac = util_Calloc ((size_t) k, sizeof (long));
   for (j = 0; j < k; j++)
      paramL->ILac[j] = I[j];

   len = strlen (gen->name);
   strncpy (name, gen->name, len);
   len2 = strlen ("\nunif01_CreateLacGen with k = ");
   len += len2;
   strncat (name, "\nunif01_CreateLacGen with k = ", len2);
   sprintf (str, "%-d", k);
   strncat (name, str, 16);
   strncat (name, ", I = (", 8);

   for (j = 0; j < k; j++) {
      sprintf (str, "%-ld", I[j]);
      strncat (name, str, 16);
      if (j < k - 1)
         strncat (name, ", ", 2);
      else
         strncat (name, ")", 1);
   }

   len = strlen (name);
   genL->name = util_Calloc (1 + len, sizeof (char));
   strncpy (genL->name, name, len);

   /* The state of the lacunary generator is simply the state of the original
      generator */
   genL->param   = paramL;
   genL->state   = gen->state;
   genL->Write   = gen->Write;
   genL->GetBits = &LacGen_Bits;
   genL->GetU01  = &LacGen_U01;
   return genL;
}

void unif01_DeleteLacGen (unif01_Gen *gen)
{
   LacGen_param *param;
   if (NULL == gen) return;
   param = gen->param;
   param->ILac = util_Free (param->ILac);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/**************************************************************************/

typedef struct {
   unif01_Gen *gen0;       /* The original generator */
   double R;               /* Total probability over [0, a) */
   double S;               /* (R - a) / (1 - a) */
   double invp;            /* Inverse of probability density over [0, a) */
   double invq;            /* Inverse of probability density over [a, 1) */
} BiasGen_param;


static double BiasGen_U01 (void *vpar, void *junk)
{
   double U;
   BiasGen_param *paramB = vpar;
   unif01_Gen *gen = paramB->gen0;

   U = gen->GetU01 (gen->param, gen->state);
   if (U < paramB->R)
      return (U * paramB->invp);
   else
      return (U - paramB->S) * paramB->invq;
}


static unsigned long BiasGen_Bits (void *vpar, void *junk)
{
   return (unsigned long) (unif01_NORM32 * BiasGen_U01 (vpar, junk));
}


unif01_Gen * unif01_CreateBiasGen (unif01_Gen *gen, double a, double R)
{
   const double Epsilon = 2.0E-16;
   unif01_Gen *genB;
   BiasGen_param *paramB;
   double p;                 /* probability density over [0, a) */
   double q;                 /* probability density over [a, 1) */
   char name[LEN0 + 1] = "";
   char str[16];
   size_t len;

   util_Assert (R >= 0.0 && R <= 1.0,
                "unif01_CreateBiasGen:   P must be in [0, 1]");
   util_Assert (a > 0.0 && a < 1.0,
                "unif01_CreateBiasGen:   a must be in (0, 1)");

   genB = util_Malloc (sizeof (unif01_Gen));
   paramB = util_Malloc (sizeof (BiasGen_param));
   paramB->gen0 = gen;

   p = R / a;
   q = (1.0 - R) / (1.0 - a);
   if (p < Epsilon)
      paramB->invp = 0.0;
   else
      paramB->invp = 1.0 / p;
   if (q < Epsilon)
      paramB->invq = 0.0;
   else
      paramB->invq = 1.0 / q;
   paramB->R = R;
   paramB->S = (p - q) * a;

   strncpy (name, gen->name, LEN0);
   len = strlen ("\nunif01_CreateBiasGen with  P = ");
   strncat (name, "\nunif01_CreateBiasGen with  P = ", len);
   sprintf (str, "%.4f", R);
   len = strlen (str);
   strncat (name, str, len);
   strncat (name, ",  a = ", 8);
   sprintf (str, "%.4f", a);
   len = strlen (str);
   strncat (name, str, len);

   len = strlen (name);
   genB->name = util_Calloc (1 + len, sizeof (char));
   strncpy (genB->name, name, len);

   /* The state of the bias generator is simply the state of the original
      generator */
   genB->param   = paramB;
   genB->state   = gen->state;
   genB->Write   = gen->Write;
   genB->GetBits = &BiasGen_Bits;
   genB->GetU01  = &BiasGen_U01;
   return genB;
}

void unif01_DeleteBiasGen (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/*************************************************************************/
typedef struct {
   unif01_Gen *gen0;               /* The original generator */
   int k;                          /* keep k numbers */
   int s;                          /* skip s numbers */
   int n;                          /* state */
} LuxGen_param;


static unsigned long LuxGen_Bits (void *vpar, void *junk)
{
   LuxGen_param *paramL = vpar;
   unif01_Gen *gen = paramL->gen0;
   if (0 == paramL->n) {
      int i;
      for (i = paramL->s; i > 0; --i)
         gen->GetBits (gen->param, gen->state);
      paramL->n = paramL->k;
   }
   --paramL->n;
   return gen->GetBits (gen->param, gen->state);
}


static double LuxGen_U01 (void *vpar, void *junk)
{
   LuxGen_param *paramL = vpar;
   unif01_Gen *gen = paramL->gen0;
   if (0 == paramL->n) {
      int i;
      for (i = paramL->s; i > 0; --i)
         gen->GetU01 (gen->param, gen->state);
      paramL->n = paramL->k;
   }
   --paramL->n;
   return gen->GetU01 (gen->param, gen->state);
}


unif01_Gen * unif01_CreateLuxGen (unif01_Gen *gen, int k, int L)
{
   unif01_Gen *genL;
   LuxGen_param *paramL;
   char name[LEN0 + 1] = "";
   char str[26];
   size_t len;
   const int s = L - k;

   util_Assert (k > 0, "unif01_CreateLuxGen:   k <= 0");
   util_Assert (k <= L, "unif01_CreateLuxGen:   L < k");

   genL = util_Malloc (sizeof (unif01_Gen));
   paramL = util_Malloc (sizeof (LuxGen_param));
   paramL->gen0 = gen;
   paramL->s = s;
   paramL->k = k;
   paramL->n = k;

   strncpy (name, gen->name, LEN0);
   len = strlen ("\nunif01_CreateLuxGen:   k = ");
   strncat (name, "\nunif01_CreateLuxGen:   k = ", len);
   sprintf (str, "%-d,   L = %-d", k, L);
   len = strlen (str);
   strncat (name, str, len);
   len = strlen (name);
   genL->name = util_Calloc (1 + len, sizeof (char));
   strncpy (genL->name, name, len);

   genL->param   = paramL;
   genL->state   = gen->state;
   genL->Write   = gen->Write;
   genL->GetBits = &LuxGen_Bits;
   genL->GetU01  = &LuxGen_U01;
   return genL;
}


void unif01_DeleteLuxGen (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/*************************************************************************/
typedef struct {
   unif01_Gen *gen0;               /* The original generator */
   unsigned long mask;             /* s most significant bits */
} TruncGen_param;


static unsigned long TruncGen_Bits (void *vpar, void *junk)
{
   TruncGen_param *paramT = vpar;
   unif01_Gen *gen = paramT->gen0;

   return paramT->mask & gen->GetBits (gen->param, gen->state);
}

static double TruncGen_U01 (void *vpar, void *vsta)
{
   return TruncGen_Bits (vpar, vsta) * unif01_INV32;
}

unif01_Gen * unif01_CreateTruncGen (unif01_Gen *gen, int b)
{
   unif01_Gen *genT;
   TruncGen_param *paramT;
   char name[LEN0 + 1] = "";
   char str[16];
   size_t len;

   if (b < 0)
      util_Error ("unif01_CreateTruncGen:   s < 0");
   if (b > 32)
      util_Error ("unif01_CreateTruncGen:   s > 32");

   genT = util_Malloc (sizeof (unif01_Gen));
   paramT = util_Malloc (sizeof (TruncGen_param));
   paramT->gen0 = gen;
   if (b >= 32)
      paramT->mask = 0xffffffffU;
   else
      paramT->mask = (0xffffffffU >> (32 - b)) << (32 - b);

   strncpy (name, gen->name, LEN0);
   len = strlen ("\nunif01_CreateTruncGen with b = ");
   strncat (name, "\nunif01_CreateTruncGen with b = ", len);
   sprintf (str, "%-d", b);
   len = strlen (str);
   strncat (name, str, len);
   strncat (name, "  bits:", 8);

   len = strlen (name);
   genT->name = util_Calloc (1 + len, sizeof (char));
   strncpy (genT->name, name, len);

   /* The state of the trunc generator is simply the state of the original
      generator */
   genT->param   = paramT;
   genT->state   = gen->state;
   genT->Write   = gen->Write;
   genT->GetBits = &TruncGen_Bits;
   genT->GetU01  = &TruncGen_U01;
   return genT;
}

void unif01_DeleteTruncGen (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/************************************************************************/
/*
 * The original generator is gen0. The r most significant bits of each random
 * number are dropped, and the s following bits are kept. 
 */
typedef struct {
   unif01_Gen *gen0;
   int nrows;   /* Number of integers used in making a new random number */
   int B;       /* Number of blocks in 1 s-bits group */
   int w;       /* Number of bits in a block */
   unsigned long maskw;   /* Mask of w bits = 2^w - 1 */
   int r;
   int s;
} BitBlock_param;


typedef struct {
   unsigned long *Z;
   int n;        /* Build n random numbers at a time, n <= 32 */
   BitBlock_param *param;
} BitBlock_state;


static unsigned long BitBlock_Bits (void *vpar, void *vsta)
{
   BitBlock_state *state = vsta;

   if (state->n <= 0) {
      BitBlock_param *param = vpar;
      unsigned long X;
      int i, j;

      /* Generate B random integers Z from the bits of nrows random integers
         X from the original gen0 */
      for (j = 0; j < param->B; j++) 
	 state->Z[j] = 0;
      for (i = 0; i < param->nrows; i++) {
 	 /* Get a random integer X of s bits */
 	 X = unif01_StripB (param->gen0, param->r, param->s);
         /* Take w of the s bits of X to make bits of Z[j] */
         for (j = 0; j < param->B; j++) {
	    state->Z[j] <<= param->w;
	    state->Z[j] |= X & param->maskw;
	    X >>= param->w;
	 }
      }
      state->n = param->B;
   }
   return state->Z[--state->n];
}


static double BitBlock_U01 (void *vpar, void *vsta)
{
   return BitBlock_Bits (vpar, vsta) * unif01_INV32;
}


static void WrBitBlock (void *vsta)
{
   BitBlock_state *state = vsta;
   state->param->gen0->Write (state->param->gen0->state);
}


unif01_Gen * unif01_CreateBitBlockGen (unif01_Gen *gen, int r, int s, int w)
{
   unif01_Gen *genV;
   BitBlock_param *paramV;
   BitBlock_state *stateV;
   char *name;
   char str[64];
   size_t len, len2, len3;
 
   util_Assert (s > 0, "unif01_CreateBitBlockGen:   s <= 0");
   util_Assert (r >= 0, "unif01_CreateBitBlockGen:   r < 0");
   util_Assert (r + s <= 32, "unif01_CreateBitBlockGen:   r + s > 32");
   util_Assert (w > 0, "unif01_CreateBitBlockGen:   w < 1");
   util_Assert (32 % w == 0, "unif01_CreateBitBlockGen:   w must divide 32");

   genV = util_Malloc (sizeof (unif01_Gen));
   paramV = util_Malloc (sizeof (BitBlock_param));
   stateV = util_Malloc (sizeof (BitBlock_state));
   paramV->gen0 = gen;
   paramV->s = s;
   paramV->r = r;
   paramV->w = w;
   paramV->B = s / w;
   paramV->maskw = num_TwoExp[paramV->w] - 1.0;
   paramV->nrows = 32 / w;
   stateV->param = paramV;
   stateV->n = 0;
   stateV->Z = util_Calloc ((size_t) paramV->B, sizeof (unsigned long));

   len = strlen (gen->name);
   len2 = strlen ("\nunif01_CreateBitBlockGen:   ");
   len += len2;
   sprintf (str, "r = %1d,   s = %1d,   w = %1d", r, s, w);
   len3 = strlen (str);
   len += len3;
   name = util_Calloc (len + 1, sizeof (char));
   strncpy (name, gen->name, len);
   strncat (name, "\nunif01_CreateBitBlockGen:   ", len2);
   strncat (name, str, len3);
   genV->name    = name;
   genV->param   = paramV;
   genV->state   = stateV;
   genV->Write   = &WrBitBlock;
   genV->GetBits = &BitBlock_Bits;
   genV->GetU01  = &BitBlock_U01;
   return genV;
}


void unif01_DeleteBitBlockGen (unif01_Gen *gen)
{
   BitBlock_state *state;
   if (NULL == gen) return;
   state = gen->state;
   state->Z = util_Free (state->Z);
   gen->param = util_Free (gen->param);
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/************************************************************************/

static double CombGen2_U01_Add (void *vpar, void *junk)
{
   unif01_Comb2_Param *g = vpar;
   unif01_Gen *gen1 = g->gen1;
   unif01_Gen *gen2 = g->gen2;
   double U;

   U = gen1->GetU01 (gen1->param, gen1->state) +
       gen2->GetU01 (gen2->param, gen2->state);
   if (U >= 1.0)
      return (U - 1.0);
   else
      return U;
}


static unsigned long CombGen2_Bits_Add (void *vpar, void *junk)
{
   return (unsigned long) (unif01_NORM32 * CombGen2_U01_Add (vpar, junk));
}


static unsigned long CombGen2_Bits_Xor (void *vpar, void *junk)
{
   unif01_Comb2_Param *g = vpar;
   unif01_Gen *gen1 = g->gen1;
   unif01_Gen *gen2 = g->gen2;

   return gen1->GetBits (gen1->param, gen1->state) ^
          gen2->GetBits (gen2->param, gen2->state);
}


static double CombGen2_U01_Xor (void *vpar, void *junk)
{
   return CombGen2_Bits_Xor (vpar, junk) * unif01_INV32;
}


static void WrCombGen2 (void *vsta)
{
   unif01_Comb2_Param *g = vsta;
   printf ("2 Combined Generators:\n");
   g->gen1->Write (g->gen1->state);
   g->gen2->Write (g->gen2->state);
}


static unif01_Gen * CreateCombGen2 (unif01_Gen *g1, unif01_Gen *g2, 
   char *mess, char *name)
{
   unif01_Gen *gen;
   unif01_Comb2_Param *paramC;
   size_t len, L;

   gen = util_Malloc (sizeof (unif01_Gen));
   paramC = util_Malloc (sizeof (unif01_Comb2_Param));
   paramC->gen1 = g1;
   paramC->gen2 = g2;

   len = strlen (g1->name) + strlen (g2->name) + strlen (name) + strlen (mess);
   len += 5;
   gen->name = util_Calloc (len + 1, sizeof (char));
   L = strlen (mess);
   if (L > 0) {
      strncpy (gen->name, mess, len);
      if (mess[L - 1] != ':')
         strncat (gen->name, ":", 3);
      strncat (gen->name, "\n", 3);
   }
   strncat (gen->name, g1->name, len);
   strncat (gen->name, "\n", 3);
   strncat (gen->name, g2->name, len);
   strncat (gen->name, name, len);

   gen->param  = paramC;
   gen->state  = paramC;
   gen->Write  = &WrCombGen2;
   return gen;
}


unif01_Gen * unif01_CreateCombAdd2 (unif01_Gen *g1, unif01_Gen *g2, char *Mess)
{
   unif01_Gen *gen;
   gen = CreateCombGen2 (g1, g2, Mess, "\nunif01_CreateCombAdd2");
   gen->GetU01 = &CombGen2_U01_Add;
   gen->GetBits = &CombGen2_Bits_Add;
   return gen;
}


unif01_Gen * unif01_CreateCombXor2 (unif01_Gen *g1, unif01_Gen *g2,
   char *Mess)
{
   unif01_Gen *gen;
   gen = CreateCombGen2 (g1, g2, Mess, "\nunif01_CreateCombXor2");
   gen->GetU01 = &CombGen2_U01_Xor;
   gen->GetBits = &CombGen2_Bits_Xor;
   return gen;
}

void unif01_DeleteCombGen (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/************************************************************************/

typedef struct {
  unif01_Gen *gen1;
  unif01_Gen *gen2;
  unif01_Gen *gen3;
} Comb3_Param;


static double CombGen3_U01_Add (void *vpar, void *junk)
{
   Comb3_Param *g = vpar;
   unif01_Gen *gen1 = g->gen1;
   unif01_Gen *gen2 = g->gen2;
   unif01_Gen *gen3 = g->gen3;
   double U;

   /*
   When the combined generator is used to generate random integers, in rare
   cases, an integer may differ by 1 unit depending on the order of
   addition of the 3 terms (one from each component). This is due
   to the last bit (bit 53) of the value returned which may be affected by
   floating-point numerical errors. Furthermore, the result
   may be different if the addition is done without function calls
   (as in the pre-programmed version of Wichmann-Hill for example in
   {\tt ulcg\_CreateCombWH3}), in which case, the 2 extra guard bits
   required by the IEEE-754 standard in floating-point arithmetic
   operations may give a more precise result.
   */
   U = gen1->GetU01 (gen1->param, gen1->state) +
       gen2->GetU01 (gen2->param, gen2->state) +
       gen3->GetU01 (gen3->param, gen3->state);

   if (U < 1.0)
      return U;
   if (U < 2.0)
      return (U - 1.0);
   return U - 2.0;
}


static unsigned long CombGen3_Bits_Add (void *vpar, void *junk)
{
   return (unsigned long) (CombGen3_U01_Add (vpar, junk) * unif01_NORM32);
}


static unsigned long CombGen3_Bits_Xor (void *vpar, void *junk)
{
   Comb3_Param *g = vpar;
   unif01_Gen *gen1 = g->gen1;
   unif01_Gen *gen2 = g->gen2;
   unif01_Gen *gen3 = g->gen3;

   return  gen1->GetBits (gen1->param, gen1->state) ^ 
           gen2->GetBits (gen2->param, gen2->state) ^
           gen3->GetBits (gen3->param, gen3->state);
}


static double CombGen3_U01_Xor (void *vpar, void *junk)
{
   return CombGen3_Bits_Xor (vpar, junk) * unif01_INV32;
}


static void WrCombGen3 (void *vsta )
{
   Comb3_Param *g = vsta;
   printf ("3 Combined Generators:\n");
   g->gen1->Write (g->gen1->state);
   g->gen2->Write (g->gen2->state);
   g->gen3->Write (g->gen3->state);
}


static unif01_Gen * CreateCombGen3 (unif01_Gen *g1, unif01_Gen *g2,
   unif01_Gen *g3, const char *mess, const char *name)
{
   unif01_Gen *gen;
   Comb3_Param *paramC;
   size_t len, L;

   gen = util_Malloc (sizeof (unif01_Gen));
   paramC = util_Malloc (sizeof (Comb3_Param));
   paramC->gen1 = g1;
   paramC->gen2 = g2;
   paramC->gen3 = g3;

   len = strlen (g1->name) + strlen (g2->name) + strlen (g3->name) +
         strlen (name) + strlen (mess);
   len += 5;
   gen->name = util_Calloc (len + 1, sizeof (char));
   L = strlen (mess);
   if (L > 0) {
      strncpy (gen->name, mess, len);
      if (mess[L - 1] != ':')
         strncat (gen->name, ":", 3);
      strncat (gen->name, "\n", 3);
   }
   strncat (gen->name, g1->name, len);
   strncat (gen->name, "\n", 3);
   strncat (gen->name, g2->name, len);
   strncat (gen->name, "\n", 3);
   strncat (gen->name, g3->name, len);
   strncat (gen->name, name, len);

   gen->param  = paramC;
   gen->state  = paramC;
   gen->Write  = &WrCombGen3;
   return gen;
}


unif01_Gen * unif01_CreateCombAdd3 (unif01_Gen *g1, unif01_Gen *g2,
   unif01_Gen *g3, char *mess)
{
   unif01_Gen *gen;
   gen = CreateCombGen3 (g1, g2, g3, mess, "\nunif01_CreateCombAdd3");
   gen->GetU01 = &CombGen3_U01_Add;
   gen->GetBits = &CombGen3_Bits_Add;
   return gen;
}


unif01_Gen * unif01_CreateCombXor3 (unif01_Gen *g1, unif01_Gen *g2,
   unif01_Gen *g3, char *mess)
{
   unif01_Gen *gen;
   gen = CreateCombGen3 (g1, g2, g3, mess, "\nunif01_CreateCombXor3");
   gen->GetU01 = &CombGen3_U01_Xor;
   gen->GetBits = &CombGen3_Bits_Xor;
   return gen;
}


/*=========================================================================*/

typedef struct {
   int j;                           /* Which random number */
   int i;                           /* Which generator */
   int L;
   int k;                           /* Number of parallel generators */
   unif01_Gen **agen;               /* Parallel generators */
} ParallelGen_state;


static double ParallelGen_U01 (void *junk, void *vsta)
{
   ParallelGen_state *stateP = vsta;
   unif01_Gen *g;

   if (++stateP->j >= stateP->L) {
      stateP->j = 0;
      if (++stateP->i >= stateP->k)
         stateP->i = 0;
   }
   g = stateP->agen[stateP->i];
   return g->GetU01 (g->param, g->state);
}


static unsigned long ParallelGen_Bits (void *junk, void *vsta)
{
   ParallelGen_state *stateP = vsta;
   unif01_Gen *g;

   if (++stateP->j >= stateP->L) {
      stateP->j = 0;
      if (++stateP->i >= stateP->k)
         stateP->i = 0;
   }
   g = stateP->agen[stateP->i];
   return g->GetBits (g->param, g->state);
}


static void WrParallelGen (void *vsta)
{
   int i;
   ParallelGen_state *state = vsta;
   printf ("   i = %d,    j = %d\n\nParallel Generators:\n", state->i, state->j);
   for (i = 0; i < state->k; ++i)
      unif01_WriteNameGen(state->agen[i]);
}


unif01_Gen * unif01_CreateParallelGen (int k, unif01_Gen *gen[], int L)
{
#define NCAT 16
   unif01_Gen *genP;
   ParallelGen_state *stateP;
   char name[LEN0 + 1] = {0};
   char str[NCAT + 1];
   size_t len;
   int j;

   genP = util_Malloc (sizeof (unif01_Gen));
   stateP = util_Malloc (sizeof (ParallelGen_state));
   stateP->k = k;
   stateP->L = L;
   stateP->i = k;
   stateP->j = L;
   stateP->agen = util_Calloc ((size_t) k, sizeof (unif01_Gen *));
   for (j = 0; j < k; j++)
      stateP->agen[j] = gen[j];

   len = strlen ("unif01_CreateParallelGen:   k = ");
   strncpy (name, "unif01_CreateParallelGen:   k = ", len);
   sprintf (str, "%-d", k);
   strncat (name, str, NCAT);
   strncat (name, ",   L = ", NCAT);
   sprintf (str, "%-d", L);
   strncat (name, str, NCAT);
   len = strlen (name);
   genP->name = util_Calloc (1 + len, sizeof (char));
   strncpy (genP->name, name, len);

   genP->state   = stateP;
   genP->Write   = &WrParallelGen;
   genP->GetBits = &ParallelGen_Bits;
   genP->GetU01  = &ParallelGen_U01;
   return genP;
#undef NCAT 
}


void unif01_DeleteParallelGen (unif01_Gen *gen)
{
   ParallelGen_state *state;
   if (NULL == gen) return;
   state = gen->state;
   state->agen = util_Free (state->agen);
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/*=========================================================================*/

static double (*externGen_U01)(void);  /* The external generator U01 */
static int coGU = 0;                       /* Counter for GU_U01 */


static double GU_U01 (void *junk1, void *junk2)
{
   return externGen_U01 ();
}


static unsigned long GU_Bits (void *junk1, void *junk2)
{
   return (unsigned long) (externGen_U01 () * unif01_NORM32);
}


static void WrExternGen (void *junk2)
{
}


unif01_Gen *unif01_CreateExternGen01 (char *name, double (*f_U01)(void*,void*),
                                      unsigned long (*b_U01)(void*,void*))
{
   unif01_Gen *gen;
   size_t leng;

   /*util_Assert (coGU == 0,
      "unif01_CreateExternGen01:   only 1 such generator can be in use");*/
   coGU++;
   gen = util_Malloc (sizeof (unif01_Gen));
   gen->state = NULL;
   gen->param = NULL;
   gen->Write = WrExternGen;
   //externGen_U01 = f_U01;
   gen->GetU01 = f_U01; //GU_U01;
   gen->GetBits = b_U01; //GU_Bits;

   if (name) {
      leng = strlen (name);
      gen->name = util_Calloc (leng + 2, sizeof (char));
      strncpy (gen->name, name, leng);
   } else {
      gen->name = util_Calloc (1, sizeof (char));
      gen->name[0] = '\0';
   }
   return gen;
}


void unif01_DeleteExternGen01 (unif01_Gen * gen)
{
   if (NULL == gen)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   coGU--;
}


/*=========================================================================*/

static unsigned int (*externGen_Bits)(void);
static int coGB = 0;                        /* Counter for GB_U01 */

static double GB_U01 (void *junk1, void *junk2)
{
   return externGen_Bits () / unif01_NORM32;
}


static unsigned long GB_Bits (void *junk1, void *junk2)
{
   return externGen_Bits ();
}


unif01_Gen *unif01_CreateExternGenBits (char *name,
    unsigned int (*f_Bits)(void))
{
   unif01_Gen *gen;
   size_t leng;

   util_Assert (coGB == 0,
      "unif01_CreateExternGenBits:   only 1 such generator can be in use");
   coGB++;
   gen = util_Malloc (sizeof (unif01_Gen));
   gen->state = NULL;
   gen->param = NULL;
   gen->Write = WrExternGen;
   externGen_Bits = f_Bits;
   gen->GetU01 = GB_U01;
   gen->GetBits = GB_Bits;

   if (name) {
      leng = strlen (name);
      gen->name = util_Calloc (leng + 2, sizeof (char));
      strncpy (gen->name, name, leng);
   } else {
      gen->name = util_Calloc (1, sizeof (char));
      gen->name[0] = '\0';
   }
   return gen;
}


void unif01_DeleteExternGenBits (unif01_Gen * gen)
{
   if (NULL == gen)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   coGB--;
}


/*=========================================================================*/

static unsigned long (*externGenLong_Bits)(void);
static int coGBL = 0;                        /* Counter for GBLong_U01 */

static double GBLong_U01 (void *junk1, void *junk2)
{
   return externGenLong_Bits () / unif01_NORM32;
}


static unsigned long GBLong_Bits (void *junk1, void *junk2)
{
   return externGenLong_Bits ();
}


unif01_Gen *unif01_CreateExternGenBitsL (char *name,
    unsigned long (*f_Bits)(void))
{
   unif01_Gen *gen;
   size_t leng;

   util_Assert (coGBL == 0,
      "unif01_CreateExternGenBitsL:   only 1 such generator can be in use");
   coGBL++;
   gen = util_Malloc (sizeof (unif01_Gen));
   gen->state = NULL;
   gen->param = NULL;
   gen->Write = WrExternGen;
   externGenLong_Bits = f_Bits;
   gen->GetU01 = GBLong_U01;
   gen->GetBits = GBLong_Bits;

   if (name) {
      leng = strlen (name);
      gen->name = util_Calloc (leng + 2, sizeof (char));
      strncpy (gen->name, name, leng);
   } else {
      gen->name = util_Calloc (1, sizeof (char));
      gen->name[0] = '\0';
   }
   return gen;
}


void unif01_DeleteExternGenBitsL (unif01_Gen * gen)
{
   if (NULL == gen)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   coGBL--;
}


/**************************************************************************/

void unif01_TimerGen (unif01_Gen *gen, unif01_TimerRec * pt, long n,
    lebool fU01)
{
   chrono_Chrono *C1;
   double U;
   unsigned long V;
   long i;

   C1 = chrono_Create ();
   if (fU01)
      for (i = 0; i < n; i++)
         U = gen->GetU01 (gen->param, gen->state);
   else
      for (i = 0; i < n; i++)
         V = gen->GetBits (gen->param, gen->state);
   pt->time = chrono_Val (C1, chrono_sec);
   pt->mean = 0.0;
   pt->n = n;
   pt->fU01 = fU01;
   pt->gen = gen;
   chrono_Delete (C1);
}

void unif01_TimerSumGen (unif01_Gen *gen, unif01_TimerRec * pt, long n,
    lebool fU01)
{
   chrono_Chrono *C1;
   double Sum = 0.0;
   unsigned long Y = 0;
   long i;

   C1 = chrono_Create ();
   if (fU01)
      for (i = 0; i < n; i++)
         Sum += gen->GetU01 (gen->param, gen->state);
   else
      for (i = 0; i < n; i++)
         Y += gen->GetBits (gen->param, gen->state);
   pt->time = chrono_Val (C1, chrono_sec);
   if (fU01)
      pt->mean = Sum / n;
   else
      pt->mean = (double) Y / n;  
   pt->n = n;
   pt->gen = gen;
   pt->fU01 = fU01;
   chrono_Delete (C1);
}

void unif01_WriteTimerRec (unif01_TimerRec *R)
{
   unif01_Gen *gen = R->gen;
   char stri [LEN1 + 1] = "";
   char *p;
   size_t len;

   printf ("\n-------------  Results of speed test  ---------------");
   printf ("\n\n Host:        ");
   if (swrite_Host)
      gdef_WriteHostName ();
   else
      printf ("\n");

   /* Print only the generator name, without the parameters or seeds. */
   /* The parameters start after the first blank; name ends with ':' */
   printf (" Generator:   ");
   len = strcspn (gen->name, ":");
   strncpy (stri, gen->name, len);
   stri [len] = '\0';
   printf ("%s", stri);
   p = strstr (gen->name, "unif01");
   while (p != NULL) {
      /* For Filters or Combined generators */
      len = strcspn (p, " \0");
      strncpy (stri, p, len);
      stri [len] = '\0';
      printf (",   %s", stri);
      p += len;
      p = strstr (p, "unif01");
   }
   if (R->fU01) {
      printf ("\n Method:      GetU01");
      if (R->mean > 0.0)
         printf ("\n Mean =       %.15f", R->mean);
   } else {
      printf ("\n Method:      GetBits");
      if (R->mean > 0.0)
         printf ("\n Mean =       %.16g", R->mean);
   }
   printf ("\n Number of calls:  %ld", R->n);
   printf ("\n Total CPU time: ");
   printf ("%6.2f sec\n\n", R->time);
}

void unif01_TimerGenWr (unif01_Gen *gen, long n, lebool fU01)
{
   unif01_TimerRec timer;
   unif01_TimerGen (gen, &timer, n, fU01);
   unif01_WriteTimerRec (&timer);
}

void unif01_TimerSumGenWr (unif01_Gen *gen, long n, lebool fU01)
{
   unif01_TimerRec timer;
   unif01_TimerSumGen (gen, &timer, n, fU01);
   unif01_WriteTimerRec (&timer);
}
