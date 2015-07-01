/*************************************************************************\
 *
 * Package:        TestU01
 * File:           uxorshift.c
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
#include "addstr.h"
#include "uxorshift.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>

#define  LEN 200

#define  MASK32 0xffffffffUL       /* Mask of 32 bits */


typedef struct {
   unsigned long *Z;
   int N;
} XorshiftC_state;

typedef struct {
   int a, b, c;
} XorshiftC_param;

typedef XorshiftC_state XorshiftD_state;


typedef struct {
   unsigned int a, b, c;
} Shift3_param;


typedef struct {
   int *a;
} XorshiftD_param;

typedef struct {
   unsigned int X[8];
   unsigned int k;
} Xorshift7_state;

typedef Xorshift7_state  Xorshift13_state;



/***************************************************************************/

static unsigned long Shift32RLL_Bits (void *vpar, void *vsta)
{
   unsigned long *p = vsta;
   unsigned long t = *p;
   Shift3_param *param = vpar;

   t ^= t >> param->a;
   t ^= t << param->b;
   t ^= t << param->c;
   *p = t & MASK32;
   return *p;
}

/*-----------------------------------------------------------------------*/

static double Shift32RLL_U01 (void *par, void *sta)
{
   return Shift32RLL_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* LRL */
static unsigned long Shift32LRL_Bits (void *vpar, void *vsta)
{
   unsigned long *p = vsta;
   unsigned long t = *p;
   Shift3_param *param = vpar;

   t ^= t << param->a;
#ifndef IS_ULONG32
   t &= MASK32;
#endif
   t ^= t >> param->b;
   t ^= t << param->c;
   *p = t & MASK32;
   return *p;
}

/*-----------------------------------------------------------------------*/

static double Shift32LRL_U01 (void *par, void *sta)
{
   return Shift32LRL_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* LLR */
static unsigned long Shift32LLR_Bits (void *vpar, void *vsta)
{
   unsigned long *p = vsta;
   unsigned long t = *p;
   Shift3_param *param = vpar;

   t ^= t << param->a;
   t ^= t << param->b;
#ifndef IS_ULONG32
   t &= MASK32;
#endif
   t ^= t >> param->c;
   *p = t;
   return *p;
}

/*-----------------------------------------------------------------------*/

static double Shift32LLR_U01 (void *par, void *sta)
{
   return Shift32LLR_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* RRL */
static unsigned long Shift32RRL_Bits (void *vpar, void *vsta)
{
   unsigned long *p = vsta;
   unsigned long t = *p;
   Shift3_param *param = vpar;

   t ^= t >> param->a;
   t ^= t >> param->b;
   t ^= t << param->c;
   *p = t & MASK32;
   return *p;
}

/*-----------------------------------------------------------------------*/

static double Shift32RRL_U01 (void *par, void *sta)
{
   return Shift32RRL_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* RLR */
static unsigned long Shift32RLR_Bits (void *vpar, void *vsta)
{
   unsigned long *p = vsta;
   unsigned long t = *p;
   Shift3_param *param = vpar;

   t ^= t >> param->a;
   t ^= t << param->b;
#ifndef IS_ULONG32
   t &= MASK32;
#endif
   t ^= t >> param->c;
   *p = t;
   return *p;
}

/*-----------------------------------------------------------------------*/

static double Shift32RLR_U01 (void *par, void *sta)
{
   return Shift32RLR_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* LRR */
static unsigned long Shift32LRR_Bits (void *vpar, void *vsta)
{
   unsigned long *p = vsta;
   unsigned long t = *p;
   Shift3_param *param = vpar;

   t ^= t << param->a;
#ifndef IS_ULONG32
   t &= MASK32;
#endif
   t ^= t >> param->b;
   t ^= t >> param->c;
   *p = t;
   return *p;
}

/*-----------------------------------------------------------------------*/

static double Shift32LRR_U01 (void *par, void *sta)
{
   return Shift32LRR_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrShift32 (void *sta)
{
   unsigned long *p = sta;
   printf (" y = %lu\n", *p);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *uxorshift_CreateXorshift32 (int a, int b, int c, unsigned int y)
{
   unif01_Gen *gen;
   unsigned long *p1;
   size_t leng;
   char name[LEN + 1];
   Shift3_param *param;

   util_Assert ((a < 32) && (a > -32),
       "uxorshift_CreateXorshift32:   a must be in [-32..32]");
   util_Assert ((b < 32) && (b > -32),
       "uxorshift_CreateXorshift32:   b must be in [-32..32]");
   util_Assert ((c < 32) && (c > -32),
       "uxorshift_CreateXorshift32:   c must be in [-32..32]");

   gen = util_Malloc (sizeof (unif01_Gen));
   p1 = gen->state = util_Malloc (sizeof (unsigned long));
   param = util_Malloc (sizeof (Shift3_param));

   *p1 = y;
   param->a = a > 0 ? a : -a;
   param->b = b > 0 ? b : -b;
   param->c = c > 0 ? c : -c;
   gen->param = param;
   gen->Write = &WrShift32;

   strcpy (name, "uxorshift_CreateXorshift32:");
   addstr_Int (name, "   a = ", a);
   addstr_Int (name, ",   b = ", b);
   addstr_Int (name, ",   c = ", c);
   addstr_Uint (name, ",   y = ", y);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (a > 0) {
      if (b > 0) {
         if (c > 0) {
	   util_Error (
           "uxorshift_CreateXorshift32:   case { <<, <<, << } not programmed");
         } else {
            gen->GetBits = &Shift32LLR_Bits;
            gen->GetU01 = &Shift32LLR_U01;
         }
      } else {
         if (c > 0) {
            gen->GetBits = &Shift32LRL_Bits;
            gen->GetU01 = &Shift32LRL_U01;
         } else {
            gen->GetBits = &Shift32LRR_Bits;
            gen->GetU01 = &Shift32LRR_U01;
         }
      }
   } else {
      if (b > 0) {
         if (c > 0) {
            gen->GetBits = &Shift32RLL_Bits;
            gen->GetU01 = &Shift32RLL_U01;
         } else {
            gen->GetBits = &Shift32RLR_Bits;
            gen->GetU01 = &Shift32RLR_U01;
         }
      } else {
         if (c > 0) {
            gen->GetBits = &Shift32RRL_Bits;
            gen->GetU01 = &Shift32RRL_U01;
         } else {
	   util_Error (
            "uxorshift_CreateXorshift32:   case { >>, >>, >> } not programmed");
         }
      }
   }
   return gen;
}


/***************************************************************************/
#ifdef USE_LONGLONG

static unsigned long Shift64RLL_Bits (void *vpar, void *vsta)
{
   ulonglong *state = vsta;
   ulonglong t = *state;
   Shift3_param *param = vpar;

   t ^= t >> param->a;
   t ^= t << param->b;
   *state = t ^ (t << param->c);
   return *state >> 32;
}

/*-----------------------------------------------------------------------*/

static double Shift64RLL_U01 (void *par, void *sta)
{
   return Shift64RLL_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* LRL */
static unsigned long Shift64LRL_Bits (void *vpar, void *vsta)
{
   ulonglong *state = vsta;
   ulonglong t = *state;
   Shift3_param *param = vpar;

   t ^= t << param->a;
   t ^= t >> param->b;
   *state = t ^ (t << param->c);
   return *state >> 32;
}

/*-----------------------------------------------------------------------*/

static double Shift64LRL_U01 (void *par, void *sta)
{
   return Shift64LRL_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* LLR */
static unsigned long Shift64LLR_Bits (void *vpar, void *vsta)
{
   ulonglong *state = vsta;
   ulonglong t = *state;
   Shift3_param *param = vpar;

   t ^= t << param->a;
   t ^= t << param->b;
   *state = t ^ (t >> param->c);
   return *state >> 32;
}

/*-----------------------------------------------------------------------*/

static double Shift64LLR_U01 (void *par, void *sta)
{
   return Shift64LLR_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* RRL */
static unsigned long Shift64RRL_Bits (void *vpar, void *vsta)
{
   ulonglong *state = vsta;
   ulonglong t = *state;
   Shift3_param *param = vpar;

   t ^= t >> param->a;
   t ^= t >> param->b;
   *state = t ^ (t << param->c);
   return *state >> 32;
}

/*-----------------------------------------------------------------------*/

static double Shift64RRL_U01 (void *par, void *sta)
{
   return Shift64RRL_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* RLR */
static unsigned long Shift64RLR_Bits (void *vpar, void *vsta)
{
   ulonglong *state = vsta;
   ulonglong t = *state;
   Shift3_param *param = vpar;

   t ^= t >> param->a;
   t ^= t << param->b;
   *state = t ^ (t >> param->c);
   return *state >> 32;
}

/*-----------------------------------------------------------------------*/

static double Shift64RLR_U01 (void *par, void *sta)
{
   return Shift64RLR_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/
/* LRR */
static unsigned long Shift64LRR_Bits (void *vpar, void *vsta)
{
   ulonglong *state = vsta;
   ulonglong t = *state;
   Shift3_param *param = vpar;

   t ^= t << param->a;
   t ^= t >> param->b;
   *state = t ^ (t >> param->c);
   return *state >> 32;
}

/*-----------------------------------------------------------------------*/

static double Shift64LRR_U01 (void *par, void *sta)
{
   return Shift64LRR_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrShift64 (void *sta)
{
   ulonglong *state = sta;
   printf (" S = %" PRIuLEAST64 "\n", *state);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *uxorshift_CreateXorshift64 (int a, int b, int c, ulonglong y)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   Shift3_param *param;
   ulonglong *state;

   util_Assert ((a < 64) && (a > -64),
       "uxorshift_CreateXorshift64:   a must be in [-64..64]");
   util_Assert ((b < 64) && (b > -64),
       "uxorshift_CreateXorshift64:   b must be in [-64..64]");
   util_Assert ((c < 64) && (c > -64),
       "uxorshift_CreateXorshift64:   c must be in [-64..64]");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (ulonglong));
   param = util_Malloc (sizeof (Shift3_param));

   *state = y;
   param->a = a > 0 ? a : -a;
   param->b = b > 0 ? b : -b;
   param->c = c > 0 ? c : -c;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrShift64;

   strcpy (name, "uxorshift_CreateXorshift64:");
   addstr_Int (name, "   a = ", a);
   addstr_Int (name, ",   b = ", b);
   addstr_Int (name, ",   c = ", c);
   addstr_ULONG (name, ",   y = ", y);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (a > 0) {
      if (b > 0) {
         if (c > 0) {
	   util_Error (
           "uxorshift_CreateXorshift64:   case { <<, <<, << } not programmed");
         } else {
            gen->GetBits = &Shift64LLR_Bits;
            gen->GetU01 = &Shift64LLR_U01;
         }
      } else {
         if (c > 0) {
            gen->GetBits = &Shift64LRL_Bits;
            gen->GetU01 = &Shift64LRL_U01;
         } else {
            gen->GetBits = &Shift64LRR_Bits;
            gen->GetU01 = &Shift64LRR_U01;
         }
      }
   } else {
      if (b > 0) {
         if (c > 0) {
            gen->GetBits = &Shift64RLL_Bits;
            gen->GetU01 = &Shift64RLL_U01;
         } else {
            gen->GetBits = &Shift64RLR_Bits;
            gen->GetU01 = &Shift64RLR_U01;
         }
      } else {
         if (c > 0) {
            gen->GetBits = &Shift64RRL_Bits;
            gen->GetU01 = &Shift64RRL_U01;
         } else {
	   util_Error (
            "uxorshift_CreateXorshift64:   case { >>, >>, >> } not programmed");
         }
      }
   }
   return gen;
}


#endif
/***************************************************************************/


static unsigned long XorshiftC_Bits (void *vpar, void *vsta)
/*
 * Only elements Z[j], j = 1, 2, ..., N are used for the state.
 * Elements Z[0] is unused.
 */
{
   XorshiftC_state *state = vsta;
   XorshiftC_param *param = vpar;
   unsigned long t;
   int j;

   if (param->a > 0)
      t = state->Z[1] ^ (state->Z[1] << param->a);
   else
      t = state->Z[1] ^ (state->Z[1] >> -param->a);

   for (j = 1; j < state->N; j++)
      state->Z[j] = state->Z[j + 1];

   if (param->b > 0) 
      t = t ^ (t << param->b);
   else {
#ifndef IS_ULONG32
      t &= MASK32;
#endif
      t = t ^ (t >> -param->b);
   }

   if (param->c > 0)
      state->Z[state->N] = state->Z[state->N] ^
                           (state->Z[state->N] << param->c) ^ t;
   else 
      state->Z[state->N] = state->Z[state->N] ^
                           (state->Z[state->N] >> -param->c) ^ t;

#ifndef IS_ULONG32
   state->Z[state->N] &= MASK32;
#endif

   return state->Z[state->N];
}


/*-----------------------------------------------------------------------*/

static double XorshiftC_U01 (void *par, void *sta)
{
   return XorshiftC_Bits (par, sta) * unif01_INV32;
}


/*-----------------------------------------------------------------------*/

static void WrXorshiftC (void *vsta)
{
   XorshiftC_state *state = vsta;
   int i;
   if (unif01_WrLongStateFlag || (state->N < 8)) {
      printf (" S = {\n ");
      for (i = 1; i <= state->N; i++) {
         printf ("   %12lu", state->Z[i]);
         if (i < state->N)
            printf (",");
         if ((i % 4) == 0)
            printf ("\n ");
      }
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}


/*-----------------------------------------------------------------------*/

unif01_Gen *uxorshift_CreateXorshiftC (int a, int b, int c, int N,
              unsigned int S[])
{
   unif01_Gen *gen;
   int i;
   size_t leng;
   char name[LEN + 1];
   XorshiftC_param *param;
   XorshiftC_state *state;

   util_Assert ((a < 32) && (a > -32),
       "uxorshift_CreateXorshiftC:   a must be in [-31..31]");
   util_Assert ((b < 32) && (b > -32),
       "uxorshift_CreateXorshiftC:   b must be in [-31..31]");
   util_Assert ((c < 32) && (c > -32),
       "uxorshift_CreateXorshiftC:   c must be in [-31..31]");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (XorshiftC_state));
   param = util_Malloc (sizeof (XorshiftC_param));

   param->a = a;
   param->b = b;
   param->c = c;

   strcpy (name, "uxorshift_CreateXorshiftC:");
   addstr_Int (name, "   a = ", a);
   addstr_Int (name, ",   b = ", b);
   addstr_Int (name, ",   c = ", c);
   addstr_Int (name, ",   N = ", N);
   addstr_ArrayUint (name, ",   S = ", N, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->N = N;
   state->Z = util_Calloc ((size_t) N + 1, sizeof (unsigned long));
   for (i = 0; i < N; i++) 
      state->Z[i + 1] = S[i];

   gen->GetBits = &XorshiftC_Bits;
   gen->GetU01 = &XorshiftC_U01;
   gen->state = state;
   gen->param = param;
   gen->Write = &WrXorshiftC;

   return gen;
}


/*-------------------------------------------------------------------------*/

void uxorshift_DeleteXorshiftC (unif01_Gen * gen)
{
   XorshiftC_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->Z);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/***************************************************************************/

static unsigned long XorshiftD_Bits (void *vpar, void *vsta)
/*
 * Only elements Z[j], j = 1, 2, ..., N are used for the state.
 * Elements Z[0] is unused.
 */
{
   XorshiftD_state *state = vsta;
   XorshiftD_param *param = vpar;
   unsigned long t = 0;
   int j;

   for (j = 1; j <= state->N; j++) {
      if (param->a[j] > 0)
	 t ^= state->Z[j] ^ (state->Z[j] << param->a[j]);
      else
	 t ^= state->Z[j] ^ (state->Z[j] >> -param->a[j]);
   }

   for (j = 1; j < state->N; j++)
      state->Z[j] = state->Z[j + 1];

#ifndef IS_ULONG32
      t &= MASK32;
#endif
   state->Z[state->N] = t;
   return t;
}


/*-----------------------------------------------------------------------*/

static double XorshiftD_U01 (void *par, void *sta)
{
   return XorshiftD_Bits (par, sta) * unif01_INV32;
}


/*-----------------------------------------------------------------------*/

static void WrXorshiftD (void *vsta)
{
   WrXorshiftC (vsta);
}


/*-----------------------------------------------------------------------*/

unif01_Gen *uxorshift_CreateXorshiftD (int N, int b[], unsigned int S[])
{
   unif01_Gen *gen;
   int i;
   size_t leng;
   char name[LEN + 1];
   XorshiftD_param *param;
   XorshiftD_state *state;


   for (i = 0; i < N; i++) {
      util_Assert ((b[i] < 32) && (b[i] > -32),
          "uxorshift_CreateXorshiftD:   all b[i] must be in [-31..31]");
   }

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (XorshiftD_state));
   param = util_Malloc (sizeof (XorshiftD_param));

   strcpy (name, "uxorshift_CreateXorshiftD:");
   addstr_Int (name, "   r = ", N);
   addstr_ArrayInt(name, ",   b = ", N, b);
   addstr_ArrayUint (name, ",   S = ", N, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->N = N;
   state->Z = util_Calloc ((size_t) N + 1, sizeof (unsigned long));
   param->a = util_Calloc ((size_t) N + 1, sizeof (int));
   for (i = 0; i < N; i++) {
      state->Z[i + 1] = S[i];
      param->a[i + 1] = b[i];
   }

   gen->GetBits = &XorshiftD_Bits;
   gen->GetU01 = &XorshiftD_U01;
   gen->state = state;
   gen->param = param;
   gen->Write = &WrXorshiftD;

   return gen;
}


/*-------------------------------------------------------------------------*/

void uxorshift_DeleteXorshiftD (unif01_Gen * gen)
{
   uxorshift_DeleteXorshiftC (gen);
}


/***************************************************************************/

static unsigned long Xorshift7_Bits (void *junk, void *vsta)
{
   Xorshift7_state *state = vsta;
   unsigned int y, t;
   t = state->X[(state->k + 7) & 0x7U];
   t = t ^ (t << 13);
   y = t ^ (t << 9);
   t = state->X[(state->k + 4) & 0x7U];
   y ^= t ^ (t << 7);
   t = state->X[(state->k + 3) & 0x7U];
   y ^= t ^ (t >> 3);
   t = state->X[(state->k + 1) & 0x7U];
   y ^= t ^ (t >> 10);
   t = state->X[state->k];
   t = t ^ (t >> 7);
   y ^= t ^ (t << 24);
   state->X[state->k] = y;
   state->k = (state->k + 1) & 0x7U;
   return y;
}


/*-------------------------------------------------------------------------*/

static double Xorshift7_U01 (void *vpar, void *vsta)
{
   return (unif01_INV32 * Xorshift7_Bits (vpar, vsta));
}


/*-------------------------------------------------------------------------*/

static void WrXorshift7 (void *vsta)
{
   Xorshift7_state *state = vsta;
   int j;
   printf (" k = %1u\n", state->k);
   printf (" X = {");
   for (j = 0; j < 8; j++) {
      printf ("  %10u", state->X[j]);
      if (j < 7)
         printf (",");
      if (j == 3)
         printf ("\n      ");
   }
   printf (" }\n\n");
}


/*-------------------------------------------------------------------------*/

unif01_Gen *uxorshift_CreateXorshift7 (unsigned int S[])
{
   unif01_Gen *gen;
   Xorshift7_state *state;
   size_t leng;
   char name[LEN + 1];
   int j;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (Xorshift7_state));
   for (j = 0; j < 8; j++)
      state->X[j] = S[j];
   state->k = 0;

   strncpy (name, "uxorshift_CreateXorshift7:", LEN);
   addstr_ArrayUint (name, "   S = ", 8, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &Xorshift7_Bits;
   gen->GetU01 = &Xorshift7_U01;
   gen->Write = &WrXorshift7;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*=========================================================================*/

static unsigned long Xorshift13_Bits (void *junk, void *vsta)
{
   Xorshift13_state *state = vsta;
   unsigned int y, t;

   t = state->X[(state->k + 7) & 0x7U];
   y = t ^ (t << 17);
   t = state->X[(state->k + 6) & 0x7U];
   y ^= t ^ (t << 10);
   t = state->X[(state->k + 4) & 0x7U];
   t = t ^ (t << 17);
   y ^= t ^ (t >> 9);
   t = state->X[(state->k + 4) & 0x7U];
   y ^= t ^ (t >> 3);
   t = state->X[(state->k + 3) & 0x7U];
   y ^= t ^ (t >> 12);
   t = state->X[(state->k + 3) & 0x7U];
   y ^= t ^ (t >> 25);
   t = state->X[(state->k + 2) & 0x7U];
   t = t ^ (t >> 2);
   y ^= t ^ (t >> 3);
   t = state->X[(state->k + 1) & 0x7U];
   y ^= t ^ (t >> 27);
   t = state->X[(state->k + 1) & 0x7U];
   y ^= t ^ (t >> 22);
   t = state->X[state->k];
   t = t ^ (t >> 3);
   y ^= t ^ (t << 24);
   state->X[state->k] = y;
   state->k = (state->k + 1) & 0x7U;
   return y;
}


/*-------------------------------------------------------------------------*/

static double Xorshift13_U01 (void *vpar, void *vsta)
{
   return (unif01_INV32 * Xorshift13_Bits (vpar, vsta));
}


/*-------------------------------------------------------------------------*/

unif01_Gen *uxorshift_CreateXorshift13 (unsigned int S[])
{
   unif01_Gen *gen;
   Xorshift13_state *state;
   size_t leng;
   char name[LEN + 1];
   int j;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (Xorshift13_state));
   for (j = 0; j < 8; j++)
      state->X[j] = S[j];
   state->k = 0;

   strncpy (name, "uxorshift_CreateXorshift13:", LEN);
   addstr_ArrayUint (name, "   S = ", 8, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &Xorshift13_Bits;
   gen->GetU01 = &Xorshift13_U01;
   gen->Write = &WrXorshift7;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*=========================================================================*/

void uxorshift_DeleteGen (unif01_Gen * gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}
