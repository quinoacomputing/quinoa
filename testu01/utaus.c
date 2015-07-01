/*************************************************************************\
 *
 * Package:        TestU01
 * File:           utaus.c
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
#include "mystr.h"
#include "addstr.h"

#include "utaus.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>




/*============================= constants ================================*/

#define LEN0 300                   /* Max length of strings */
#define LEN1 100                   /* Max length of strings */




/*================================ Types ================================*/

typedef struct {
   unsigned int M1, S1, Q1, K1mS1;
   unsigned int J;
} Taus_param;

typedef struct {
   unsigned int ST1;
} Taus_state;

typedef struct {
   unsigned int M1, S1, Q1, K1mS1;
   unsigned int M2, S2, Q2, K2mS2;
} CombTaus2_param;

typedef struct {
   unsigned int ST1, ST2;
} CombTaus2_state;

typedef struct {
   unsigned int M1, S1, Q1, K1mS1;
   unsigned int M2, S2, Q2, K2mS2;
   unsigned int M3, S3, Q3, K3mS3;
} CombTaus3_param;

typedef struct {
   unsigned int ST1, ST2, ST3;
} CombTaus3_state;


#ifdef USE_LONGLONG
/* The 64 bits Tausworthe */
typedef struct {
   ulonglong M1, S1, Q1, K1mS1;
} LongTaus_param;

typedef struct {
   ulonglong ST1;
} LongTaus_state;
#endif


/*============================= Functions ===============================*/

static unsigned long Taus_Bits (void *vpar, void *vsta)
{
   Taus_param *param = vpar;
   Taus_state *state = vsta;
   unsigned int A, B;

   A = (state->ST1 & param->M1) << param->S1;
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> param->K1mS1;
   state->ST1 = A ^ B;
   return state->ST1;
}

static double Taus_U01 (void *vpar, void *vsta)
{
   return Taus_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static unsigned long TausJ_Bits (void *vpar, void *vsta)
/*
 * Tausworthe generator that goes forward J steps for each value returned.
 */
{
   Taus_param *param = vpar;
   Taus_state *state = vsta;
   unsigned int A, B;
   unsigned int i;

   for (i = 1; i <= param->J; i++) {
      A = (state->ST1 & param->M1) << param->S1;
      B = ((state->ST1 << param->Q1) ^ state->ST1) >> param->K1mS1;
      state->ST1 = A ^ B;
   }
   return state->ST1;
}

static double TausJ_U01 (void *vpar, void *vsta)
{
   return TausJ_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrTaus (void *vsta)
{
   Taus_state *state = vsta;
   printf (" S = %1u\n", state->ST1);
}

/*-------------------------------------------------------------------------*/

static unif01_Gen * CreateTaus_0 (char *na, unsigned int k, unsigned int q,
   unsigned int s, unsigned int Y)
{
   unif01_Gen *gen;
   Taus_param *param;
   Taus_state *state;
   size_t len;
   char name[LEN0 + 1];
   char str[LEN1 + 1];
   unsigned int B;

   strncpy (str, na, (size_t) LEN1);
   strncat (str, ":   Invalid Parameter", (size_t) LEN1 - 30);
   util_Assert ((k <= 32) && (k > 2 * q) && (s <= k - q) &&
                (s >= 1) && (q >= 1), str);
   /* These restrictions implies, for 32 bits int, */
   /* 0 < k <= 32 ; 0 < q < 16 ; 0 < s < 32 ; 0 < k - s < 32 */

   gen   = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Taus_param));
   state = util_Malloc (sizeof (Taus_state));

   strncpy (name, na, (size_t) LEN0);
   addstr_Uint (name, ":   k = ", k);
   addstr_Uint (name, ",  q = ", q);
   addstr_Uint (name, ",  s = ", s);
   addstr_Uint (name, ",  Y = ", Y);
   len = strlen (name);
   gen->name = util_Calloc (len + 1, sizeof (char));
   strncpy (gen->name, name, len);

   param->Q1 = q;
   param->K1mS1 = k - s;
   param->S1 = s;

  /* k most signif. bits at 1 */
   B = num_TwoExp[32 - k] - 1.0;
   param->M1 = ~B;

   util_Assert (param->M1 > 0, "CreateTaus_0:   M1 = 0");

   strncpy (str, na, (size_t) LEN1);
   strncat (str, ":   Y = 0", (size_t) LEN1 - 30);
   util_Assert (Y > 0, str);
   state->ST1 = Y & param->M1;

   /* make sure that the initial state is not 0 */
   while (state->ST1 == 0) {
      Y *= 2;
      state->ST1 = Y & param->M1;
   }
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> k;
   if (k >= 8 * sizeof (int))
      B = 0;                      /* B = B >> 32 does not work */
   state->ST1 ^= B;

   gen->GetBits = &Taus_Bits;
   gen->GetU01  = &Taus_U01;
   gen->Write   = &WrTaus;
   gen->param   = param;
   gen->state   = state;
   return gen;
}

/*-------------------------------------------------------------------------*/

unif01_Gen *utaus_CreateTaus (unsigned int k, unsigned int q, unsigned int s,
   unsigned int Y)
{
   return CreateTaus_0 ("utaus_CreateTaus", k, q, s, Y);
}

unif01_Gen *utaus_CreateTausJ (unsigned int k, unsigned int q,
   unsigned int s, unsigned int j, unsigned int Y)
{
   unif01_Gen *gen;
   Taus_param *param;
   unsigned int pos;
   int found;
   size_t len;
   char str[LEN1 + 1] = "";

   gen = CreateTaus_0 ("utaus_CreateTausJ", k, q, s, Y);
   param = gen->param;
   param->J = j;
   gen->GetBits = &TausJ_Bits;
   gen->GetU01  = &TausJ_U01;

   addstr_Uint (str, ",  j = ", j);
   len = strlen (gen->name) + strlen (str);
   gen->name = util_Realloc (gen->name, (len + 1) * sizeof (char));
   mystr_Position (",  Y =", gen->name, 0, &pos, &found);
   mystr_Insert (gen->name, str, pos);
   return gen;
}


/**************************************************************************/
#ifdef USE_LONGLONG

static unsigned long LongTaus_Bits (void *vpar, void *vsta)
/* 
 * Implementation of a Tausworthe generator for k <= 64. However, bits
 * 33 to 64 will be right-shifted 32 bits before being returned in an
 * unsigned int of 32 bits. This is to satisfy the return type
 * of unif01_Gen.GetBits.
 */
{
   LongTaus_param *param = vpar;
   LongTaus_state *state = vsta;
   ulonglong A, B;
   A = (state->ST1 & param->M1) << param->S1;
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> param->K1mS1;
   state->ST1 = A ^ B;
   return state->ST1 >> 32;
}

static double LongTaus_U01 (void *vpar, void *vsta)
{
   return LongTaus_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrLongTaus (void *vsta)
{
   LongTaus_state *state = vsta;
   printf (" S = %1" PRIuLEAST64 "\n", state->ST1);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utaus_CreateLongTaus (unsigned int k, unsigned int q,
                                   unsigned int s, ulonglong Y)
{
   unif01_Gen *gen;
   LongTaus_param *param;
   LongTaus_state *state;
   size_t len;
   char name[LEN0 + 1];
   ulonglong B;

   util_Assert ((k <= 64) && (k > 2 * q) && (s <= k - q) && (s >= 1) &&
      (q >= 1), "utaus_CreateLongTaus:   Invalid Parameter");

   gen   = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LongTaus_param));
   state = util_Malloc (sizeof (LongTaus_state));

   strncpy (name, "utaus_CreateLongTaus", (size_t) LEN0);
   addstr_Uint (name, ":   k = ", k);
   addstr_Uint (name, ",  q = ", q);
   addstr_Uint (name, ",  s = ", s);
   addstr_ULONG (name, ",  Y = ", Y);
   len = strlen (name);
   gen->name = util_Calloc (len + 1, sizeof (char));
   strncpy (gen->name, name, len);

   param->Q1 = q;
   param->K1mS1 = k - s;
   param->S1 = s;
   B = num_TwoExp[64 - k] - 1.0;        /* k most signif. bits at 1 */
   param->M1 = ~B;

   util_Assert (Y > 0, "utaus_CreateLongTaus:   Y = 0");
   state->ST1 = Y & param->M1;

   /* make sure that the initial state is not 0 */
   while (state->ST1 == 0) {
      Y *= 2;
      state->ST1 = Y & param->M1;
   }
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> k;
   state->ST1 ^= B;

   gen->GetBits = &LongTaus_Bits;
   gen->GetU01  = &LongTaus_U01;
   gen->Write   = &WrLongTaus;
   gen->param   = param;
   gen->state   = state;
   return gen;
}

#endif

/**************************************************************************/

static unsigned long CombTaus2_Bits (void *vpar, void *vsta)
{
   CombTaus2_param *param = vpar;
   CombTaus2_state *state = vsta;
   unsigned int A, B;

   A = (state->ST1 & param->M1) << param->S1;
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> param->K1mS1;
   state->ST1 = A ^ B;

   A = (state->ST2 & param->M2) << param->S2;
   B = ((state->ST2 << param->Q2) ^ state->ST2) >> param->K2mS2;
   state->ST2 = A ^ B;

   return state->ST1 ^ state->ST2;
}

static double CombTaus2_U01 (void *vpar, void *vsta)
{
   return CombTaus2_Bits (vpar, vsta) * unif01_INV32;
}

static void WrCombTaus2 (void *vsta)
{
   CombTaus2_state *state = vsta;
   printf (" S1 = %1u,   S2 = %1u\n", state->ST1, state->ST2);
}


/*-------------------------------------------------------------------------*/

unif01_Gen * utaus_CreateCombTaus2 (unsigned int k1, unsigned int k2,
   unsigned int q1, unsigned int q2, unsigned int s1, unsigned int s2,
   unsigned int Y1, unsigned int Y2)
{
   unif01_Gen *gen;
   CombTaus2_param *param;
   CombTaus2_state *state;
   size_t len;
   char name[LEN0 + 1];
   unsigned int B;

   util_Assert (
      (k1 <= 32) && (k1 > 2 * q1) && (s1 <= k1 - q1) && (s1 > 0) &&
      (k2 <= 32) && (k2 > 2 * q2) && (s2 <= k2 - q2) && (s2 > 0) &&
      (q1 > 0) && (q2 > 0) && (k1 >= k2),
      "utaus_CreateCombTaus2:   Invalid Parameter");

   strncpy (name, "utaus_CreateCombTaus2:", (size_t) LEN0);
   addstr_Uint (name, "   (k1, k2) = ", k1);
   addstr_Uint (name, ", ", k2);
   addstr_Uint (name, ",   (q1, q2) = ", q1);
   addstr_Uint (name, ", ", q2);
   addstr_Uint (name, ",   (s1, s2) = ", s1);
   addstr_Uint (name, ", ", s2);
   addstr_Uint (name, ",   (Y1, Y2) = ", Y1);
   addstr_Uint (name, ", ", Y2);

   gen = util_Malloc (sizeof (unif01_Gen));
   param = gen->param = util_Malloc (sizeof (CombTaus2_param));
   state = gen->state = util_Malloc (sizeof (CombTaus2_state));
   gen->GetU01  = &CombTaus2_U01;
   gen->GetBits = &CombTaus2_Bits;
   gen->Write   = &WrCombTaus2;

   len = strlen (name);
   gen->name = util_Calloc (len + 1, sizeof (char));
   strncpy (gen->name, name, len);

   param->Q1 = q1;
   param->S1 = s1;
   param->K1mS1 = k1 - s1;
   param->M1 = num_TwoExp[32 - k1] - 1; /* k1 most sig. bits at 1 */
   param->M1 = ~param->M1;

   param->Q2 = q2;
   param->S2 = s2;
   param->K2mS2 = k2 - s2;
   param->M2 = num_TwoExp[32 - k2] - 1; /* k2 most sig. bits at 1 */
   param->M2 = ~param->M2;

   /* Initialisation as in L'Ecuyer 1996. */
   util_Assert (Y1 > 0, "utaus_CreateCombTaus2:   seed1 = 0");
   state->ST1 = Y1 & param->M1;
   while (state->ST1 == 0) {      /* Make sure the initial state is not 0 */
      Y1 = Y1 * 2;
      state->ST1 = Y1 & param->M1;
   }
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> k1;
   if (k1 >= 8 * sizeof (int))
      B = 0;                      /* B = B >> 32 does not work */
   state->ST1 = state->ST1 ^ B;

   util_Assert (Y2 > 0, "utaus_CreateCombTaus2:   seed2 = 0");
   state->ST2 = Y2 & param->M2;
   while (state->ST2 == 0) {      /* Make sure the initial state is not 0 */
      Y2 = Y2 * 2;
      state->ST2 = Y2 & param->M2;
   }
   B = ((state->ST2 << param->Q2) ^ state->ST2) >> k2;
   if (k2 >= 8 * sizeof (int))
      B = 0;                      /* B = B >> 32 does not work */
   state->ST2 = state->ST2 ^ B;

   return gen;
}


/****************************************************************************/

static unsigned long CombTaus3_Bits (void *vpar, void *vsta)
{
   CombTaus3_param *param = vpar;
   CombTaus3_state *state = vsta;
   unsigned int A, B;

   A = (state->ST1 & param->M1) << param->S1;
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> param->K1mS1;
   state->ST1 = A ^ B;

   A = (state->ST2 & param->M2) << param->S2;
   B = ((state->ST2 << param->Q2) ^ state->ST2) >> param->K2mS2;
   state->ST2 = A ^ B;

   A = (state->ST3 & param->M3) << param->S3;
   B = ((state->ST3 << param->Q3) ^ state->ST3) >> param->K3mS3;
   state->ST3 = A ^ B;

   return state->ST1 ^ state->ST2 ^ state->ST3;
}

static double CombTaus3_U01 (void *vpar, void *vsta)
{
   return CombTaus3_Bits (vpar, vsta) * unif01_INV32;
}

static void WrCombTaus3 (void *vsta)
{
   CombTaus3_state *state = vsta;
   printf (" S1 = %1u,   S2 = %1u,   S3 =  %1u\n",
            state->ST1, state->ST2, state->ST3);
}


/*-------------------------------------------------------------------------*/

unif01_Gen * utaus_CreateCombTaus3 (unsigned int k1, unsigned int k2,
   unsigned int k3, unsigned int q1, unsigned int q2, unsigned int q3,
   unsigned int s1, unsigned int s2, unsigned int s3, unsigned int Y1,
   unsigned int Y2, unsigned int Y3)
{
   unif01_Gen *gen;
   CombTaus3_param *param;
   CombTaus3_state *state;
   size_t len;
   char name[LEN0 + 1];
   unsigned int B;

   util_Assert (
      (k1 <= 32) && (k1 > 2 * q1) && (s1 <= k1 - q1) && (s1 > 0) &&
      (k2 <= 32) && (k2 > 2 * q2) && (s2 <= k2 - q2) && (s2 > 0) &&
      (k3 <= 32) && (k3 > 2 * q3) && (s3 <= k3 - q3) && (s3 > 0) &&
      (q1 > 0) && (q2 > 0) && (q3 > 0) && (k1 >= k2) && (k2 >= k3),
      "utaus_CreateCombTaus3:   Invalid Parameter");

   strncpy (name, "utaus_CreateCombTaus3:", (size_t) LEN0);
   addstr_Uint (name, "   (k1, k2, k3) = ", k1);
   addstr_Uint (name, ", ", k2);
   addstr_Uint (name, ", ", k3);
   addstr_Uint (name, ",   (q1, q2, q3) = ", q1);
   addstr_Uint (name, ", ", q2);
   addstr_Uint (name, ", ", q3);
   addstr_Uint (name, ",   (s1, s2, s3) = ", s1);
   addstr_Uint (name, ", ", s2);
   addstr_Uint (name, ", ", s3);
   addstr_Uint (name, ",   (Y1, Y2, Y3) = ", Y1);
   addstr_Uint (name, ", ", Y2);
   addstr_Uint (name, ", ", Y3);

   gen = util_Malloc (sizeof (unif01_Gen));
   param = gen->param = util_Malloc (sizeof (CombTaus3_param));
   state = gen->state = util_Malloc (sizeof (CombTaus3_state));
   gen->GetU01  = &CombTaus3_U01;
   gen->GetBits = &CombTaus3_Bits;
   gen->Write   = &WrCombTaus3;

   len = strlen (name);
   gen->name = util_Calloc (len + 1, sizeof (char));
   strncpy (gen->name, name, len);

   param->Q1 = q1;
   param->S1 = s1;
   param->K1mS1 = k1 - s1;
   param->M1 = num_TwoExp[32 - k1] - 1; /* k1 most sig. bits at 1 */
   param->M1 = ~param->M1;

   param->Q2 = q2;
   param->S2 = s2;
   param->K2mS2 = k2 - s2;
   param->M2 = num_TwoExp[32 - k2] - 1; /* k2 most sig. bits at 1 */
   param->M2 = ~param->M2;

   param->Q3 = q3;
   param->S3 = s3;
   param->K3mS3 = k3 - s3;
   param->M3 = num_TwoExp[32 - k3] - 1; /* k3 most sig. bits at 1 */
   param->M3 = ~param->M3;

   /* Initialisation as in L'Ecuyer 1996. */
   util_Assert (Y1 > 0, "utaus_CreateCombTaus3:   seed1 = 0");
   state->ST1 = Y1 & param->M1;
   while (state->ST1 == 0) {      /* Make sure the initial state is not 0 */
      Y1 = Y1 * 2;
      state->ST1 = Y1 & param->M1;
   }
   B = ((state->ST1 << param->Q1) ^ state->ST1) >> k1;
   if (k1 >= 8 * sizeof (int))
      B = 0;                      /* B = B >> 32 does not work */
   state->ST1 = state->ST1 ^ B;

   util_Assert (Y2 > 0, "utaus_CreateCombTaus3:   seed2 = 0");
   state->ST2 = Y2 & param->M2;
   while (state->ST2 == 0) {      /* Make sure the initial state is not 0 */
      Y2 = Y2 * 2;
      state->ST2 = Y2 & param->M2;
   }
   B = ((state->ST2 << param->Q2) ^ state->ST2) >> k2;
   if (k2 >= 8 * sizeof (int))
      B = 0;                      /* B = B >> 32 does not work */
   state->ST2 = state->ST2 ^ B;

   util_Assert (Y3 > 0, "utaus_CreateCombTaus3:   seed3 = 0");
   state->ST3 = Y3 & param->M3;
   while (state->ST3 == 0) {      /* Make sure the initial state is not 0 */
      Y3 = Y3 * 2;
      state->ST3 = Y3 & param->M3;
   }
   B = ((state->ST3 << param->Q3) ^ state->ST3) >> k3;
   if (k3 >= 8 * sizeof (int))
      B = 0;                      /* B = B >> 32 does not work */
   state->ST3 = state->ST3 ^ B;

   return gen;
}


/****************************************************************************/

#define FAC32  4294967296.0            /* 2^32 */
#define FAC49  562949953421312.0       /* 2^49 */
#define FAC66  73786976294838206464.0  /* 2^66 */

static double CombTaus3T_U01 (void *vpar, void *vsta)
{
   double U;
   U =  CombTaus3_Bits (vpar, vsta) / FAC32;
   U += CombTaus3_Bits (vpar, vsta) / FAC49;
   U += CombTaus3_Bits (vpar, vsta) / FAC66;
   if (U < 1.0)
      return U;
   else
      return U - 1.0;
}

/*-------------------------------------------------------------------------*/

static unsigned long CombTaus3T_Bits (void *vpar, void *vsta)
{
  return (unsigned long) (CombTaus3T_U01 (vpar, vsta) * unif01_NORM32);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * utaus_CreateCombTaus3T (unsigned int k1, unsigned int k2,
   unsigned int k3, unsigned int q1, unsigned int q2, unsigned int q3,
   unsigned int s1, unsigned int s2, unsigned int s3, unsigned int Y1,
   unsigned int Y2, unsigned int Y3)
{
   unif01_Gen *gen;
   size_t j;

   gen = utaus_CreateCombTaus3 (
          k1, k2, k3, q1, q2, q3, s1, s2, s3, Y1, Y2, Y3);
   j = strlen (gen->name);
   gen->name = util_Realloc (gen->name, (j + 2) * sizeof (char));
   j = strcspn (gen->name, ":");
   mystr_Insert (gen->name, "T", (unsigned int) j);
   gen->GetU01  = &CombTaus3T_U01;
   gen->GetBits = &CombTaus3T_Bits;
 
   return gen;
}


/****************************************************************************/

void utaus_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
