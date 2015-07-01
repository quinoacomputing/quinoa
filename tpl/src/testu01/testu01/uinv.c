/*************************************************************************\
 *
 * Package:        TestU01
 * File:           uinv.c
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

#include "uinv.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>




/*============================= Constants ===============================*/

#define  LEN  300                  /* Max string length */
#define  MASK31 0x7fffffffUL       /* Mask of 31 bits */
#define  MASK32 0xffffffffUL       /* Mask of 32 bits */




/*============================== Types ==================================*/
/* 
 * Different styles for generators of same type: Small, Medium, Large
 */

typedef enum {
   StyleS, StyleM, StyleL
} GenStyle;

/*-------------------------------------------------------------------------*/

typedef struct {
   long A1, A2, M, Q, R;
   double Norm;
} InvImpl_param;

typedef struct {
   long Z;
   GenStyle Style;
} InvImpl_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long A1, A2, Mask, Shift, M;
} InvImpl2a_param;

typedef struct {
   unsigned long Z;
} InvImpl2a_state;

/*-------------------------------------------------------------------------*/

typedef InvImpl2a_param InvImpl2b_param;

typedef InvImpl2a_state InvImpl2b_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   long M, MmA, A2;
   double Norm;
} InvExpl_param;

typedef struct {
   long Z;
} InvExpl_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long A2, E, Mask, Shift;
} InvExpl2a_param;

typedef struct {
   unsigned long Z;
} InvExpl2a_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long A1, A2, E, Mask, Shift;
} InvExpl2b_param;

typedef struct {
   unsigned long N;
} InvExpl2b_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   long *A, *Q, *R;
   long M;
   double Norm;
} InvMRG_param;

typedef struct {
   long *S;
   int Order;
} InvMRG_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   double *A;
   double M;
   double Norm;
} InvMRGFloat_param;

typedef struct {
   double *S;
   int Order;
} InvMRGFloat_state;





/*============================= Functions ===============================*/


/**********************************************************************
 *
 * Non-linear Generator with inversion in the form
 * Z = (a1 + a2 * Zm1^(-1) ) % m
 * as mentionned in L'Ecuyer (1990).
 *
 * The values returned by the generators implemented below are
 * in the following set: { 0/m, 1/m, 2/m, ..., (m-1)/m }.
 *
 *********************************************************************/

static double SmallInvImpl_U01 (void *vpar, void *vsta)
/*
 * Direct Implementation used when (A1 + A2*(M - 1)) fits in a long int
 */
{
   InvImpl_param *param = vpar;
   InvImpl_state *state = vsta;

   if (state->Z == 0)
      state->Z = param->A1;
   else
      state->Z = (param->A1 + param->A2 * num_InvEuclid (param->M, state->Z))
         % param->M;
   return (state->Z * param->Norm);
}

static unsigned long SmallInvImpl_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SmallInvImpl_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double MediumInvImpl_U01 (void *vpar, void *vsta)
/*
 * Implementation used when a2 * (m % a2) < m
 */
{
   InvImpl_param *param = vpar;
   InvImpl_state *state = vsta;
   long Zinv, k;

   if (state->Z == 0)
      state->Z = param->A1;
   else {
      Zinv = num_InvEuclid (param->M, state->Z);
      k = Zinv / param->Q;
      state->Z = param->A2 * (Zinv - k * param->Q) - k * param->R;
      if (state->Z < 0)
         state->Z += param->A1;
      else
         state->Z = (state->Z - param->M) + param->A1;
      if (state->Z < 0)
         state->Z += param->M;
   }
   return (state->Z * param->Norm);
}

static unsigned long MediumInvImpl_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumInvImpl_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double LargeInvImpl_U01 (void *vpar, void *vsta)
/*
 * Implementation used when the last 2 above cannot be used
 */
{
   InvImpl_param *param = vpar;
   InvImpl_state *state = vsta;

   if (state->Z == 0)
      state->Z = param->A1;
   else
      state->Z = num_MultModL (param->A2,
         num_InvEuclid (param->M, state->Z), param->A1, param->M);
   return (state->Z * param->Norm);
}

static unsigned long LargeInvImpl_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LargeInvImpl_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrInvImpl (void *vsta)
{
   InvImpl_state *state = vsta;

   if (state->Style == StyleS)
      printf ("Small InvImpl,");
   else if (state->Style == StyleM)
      printf ("Medium InvImpl,");
   else if (state->Style == StyleL)
      printf ("Large InvImpl,");
   printf (" Z = %1ld\n", state->Z);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvImpl (long m, long a1, long a2, long z0)
{
   unif01_Gen *gen;
   InvImpl_param *param;
   InvImpl_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((z0 < 0) || (z0 >= m) || (a1 <= 0) || (a1 >= m) ||
      (a2 <= 0) || (a2 >= m) || (m < 2) || ((m % 2) == 0))
      util_Error ("uinv_CreateInvImpl:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvImpl_param));
   state = util_Malloc (sizeof (InvImpl_state));

   strcpy (name, "uinv_CreateInvImpl:");
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   z0 = ", z0);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Norm = 1.0 / m;
   param->M = m;
   param->A1 = a1;
   param->A2 = a2;
   state->Z = z0;

   if (m - 1 <= (LONG_MAX - a1) / a2) {  /* First case */
      state->Style = StyleS;
      gen->GetBits = &SmallInvImpl_Bits;
      gen->GetU01 = &SmallInvImpl_U01;

   } else {
      param->Q = m / a2;
      param->R = m % a2;
      if (param->R <= param->Q) {     /* Second case  */
         state->Style = StyleM;
         gen->GetBits = &MediumInvImpl_Bits;
         gen->GetU01 = &MediumInvImpl_U01;

      } else {                        /* Third case */
         state->Style = StyleL;
         gen->GetBits = &LargeInvImpl_Bits;
         gen->GetU01 = &LargeInvImpl_U01;
      }
   }
   gen->param = param;
   gen->state = state;
   gen->Write = &WrInvImpl;
   return gen;
}


/**************************************************************************/
/*
 * Non-linear Generator with inversion in the form
 * Z = (a1 + a2 * Zm1^(-1) ) % m   avec   m = 2^e
 * as mentionned in Eichenauer-Herrmann (1992).
 *
 * The values returned by the generators implemented below are
 * in the following set: { 1/m, 3/m, 5/m, ..., (m-1)/m }.
 *
 * N.B.  For e = 31 ou 32, the inverse is computed according to the formula:
 * X^(-1) = X^(m-1) = X^((m / 4) - 1)
 */


static unsigned long InvImpl2a_Bits (void *vpar, void *vsta)
/* 
 * Implementation used when M = 2^e  with  3 <= e <= 30.
 * Compute the inverse by the modified Euclide algorithm.
 */
{
   InvImpl2a_param *param = vpar;
   InvImpl2a_state *state = vsta;

   state->Z = (param->A1 + param->A2 *
      num_InvEuclid ((long) param->M, (long) state->Z)) & param->Mask;
   return (state->Z << param->Shift);
}

/*-------------------------------------------------------------------------*/

static double InvImpl2a_U01 (void *vpar, void *vsta)
{
   return InvImpl2a_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static unsigned long InvImpl2a31_Bits (void *vpar, void *vsta)
/*
 * Special Implementation used when M = 2^31.
 * The inverse must be computed by exponentiation even if it is slower.
 * Given the special value of the exponent ( 2^29 - 1 )
 * ( = 536870911 ), it is done in a very direct way ...
 */
{
   InvImpl2a_param *param = vpar;
   InvImpl2a_state *state = vsta;
   unsigned int res = state->Z;
   int fc;

   for (fc = 1; fc <= 28; fc++)
      res = res * res * state->Z;          /* res = Z^(-1) */

   res &= MASK31;
   state->Z = (param->A1 + param->A2 * res) & MASK31;
   return (state->Z << 1);
}

/*-------------------------------------------------------------------------*/

static double InvImpl2a31_U01 (void *vpar, void *vsta)
{
   return InvImpl2a31_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static unsigned long InvImpl2a32_Bits (void *vpar, void *vsta)
/*
 * Special Implementation used when M = 2^32.
 * The inverse must be computed by exponentiation even if it is slower.
 * Given the special value of the exponent ( 2^30 - 1 )
 * ( = 1073741823 ), it is done in a very direct way ...
 */
{
   InvImpl2a_param *param = vpar;
   InvImpl2a_state *state = vsta;
   unsigned int res = state->Z;
   int fc;

   for (fc = 1; fc <= 29; fc++)
      res = res * res * state->Z;              /* res = Z^(-1) */
   res &= MASK32;

   state->Z = (param->A1 + param->A2 * res) & MASK32;
   return state->Z;
}

/*-------------------------------------------------------------------------*/

static double InvImpl2a32_U01 (void *vpar, void *vsta)
{
   return InvImpl2a32_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrInvImpl2a (void *vsta)
{
   InvImpl2a_state *state = vsta;
   printf (" Z = %1lu\n", state->Z);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvImpl2a (int e, unsigned long a1, unsigned long a2,
   unsigned long z0)
{
   unif01_Gen *gen;
   InvImpl2a_param *param;
   InvImpl2a_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned long M;

   if (((z0 % 2) == 0) || ((a2 % 2) == 0) || ((a1 % 2) != 0) ||
      (e < 3) || (e > 32))
      util_Error ("uinv_CreateInvImpl2a:   Invalid parameter *");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvImpl2a_param));
   state = util_Malloc (sizeof (InvImpl2a_state));

   if (e < 32)
      M = num_TwoExp[e];
   if ((e <= 31) && ((z0 >= M) || (a1 >= M) || (a2 >= M)))
      util_Error ("uinv_CreateInvImpl2a:   Invalid parameter **");

   strcpy (name, "uinv_CreateInvImpl2a:");
   addstr_Long (name, "   e = ", e);
   addstr_Ulong (name, ",   a1 = ", a1);
   addstr_Ulong (name, ",   a2 = ", a2);
   addstr_Ulong (name, ",   z0 = ", z0);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Shift = 32 - e;
   param->M = M;
   param->A1 = a1;
   param->A2 = a2;
   state->Z = z0;
   if (e < 32)
      param->Mask = M - 1;
   else
      param->Mask = MASK32;

   if (e == 32) {
      gen->GetBits = &InvImpl2a32_Bits;
      gen->GetU01 = &InvImpl2a32_U01;

   } else if (e == 31) {
      gen->GetBits = &InvImpl2a31_Bits;
      gen->GetU01 = &InvImpl2a31_U01;

   } else {                       /* e <= 30 */
      gen->GetBits = &InvImpl2a_Bits;
      gen->GetU01 = &InvImpl2a_U01;
   }
   gen->param = param;
   gen->state = state;
   gen->Write = &WrInvImpl2a;
   return gen;
}


/**************************************************************************/

static unsigned long InvImpl2b_Bits (void *vpar, void *vsta)
{
   InvImpl2b_param *param = vpar;
   InvImpl2b_state *state = vsta;
   unsigned long div0;             /* Greatest common divisor of Z and 2^e */
   unsigned long mod0;             /* Remainder of division of Zc by div0 */

   mod0 = state->Z % 2;
   div0 = 1;
   while ((mod0 == 0) && (state->Z != 0)) {
      state->Z /= 2;
      div0 *= 2;
      mod0 = state->Z % 2;
   }
   /* Compute the next state */
   state->Z = param->A1 + div0 * param->A2 *
      num_InvEuclid ((long) param->M, (long) state->Z);
   state->Z &= param->Mask;
   return state->Z << param->Shift;
}

/*-------------------------------------------------------------------------*/

static double InvImpl2b_U01 (void *vpar, void *vsta)
{
   return InvImpl2b_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static unsigned long InvImpl2b31_Bits (void *vpar, void *vsta)
{
   InvImpl2b_param *param = vpar;
   InvImpl2b_state *state = vsta;
   unsigned long div0;             /* Greatest common divisor of Z and 2^e */
   unsigned long mod0;             /* Remainder of division of Zc by div0 */

   mod0 = state->Z % 2;
   div0 = 1;
   while ((mod0 == 0) && (state->Z != 0)) {
      state->Z /= 2;
      div0 *= 2;
      mod0 = state->Z % 2;
   }
   /* Compute the next state */
   state->Z = param->A1 + div0 * param->A2 * num_InvExpon (31, state->Z);
   state->Z &= MASK31;
   return state->Z << 1;
}

/*-------------------------------------------------------------------------*/

static double InvImpl2b31_U01 (void *vpar, void *vsta)
{
   return InvImpl2b31_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static unsigned long InvImpl2b32_Bits (void *vpar, void *vsta)
{
   InvImpl2b_param *param = vpar;
   InvImpl2b_state *state = vsta;
   unsigned long div0;             /* Greatest common divisor of Z and 2^e */
   unsigned long mod0;             /* Remainder of division of Zc by div0 */

   mod0 = state->Z % 2;
   div0 = 1;
   while ((mod0 == 0) && (state->Z != 0)) {
      state->Z /= 2;
      div0 *= 2;
      mod0 = state->Z % 2;
   }
   /* Compute the next state */
   state->Z = param->A1 + div0 * param->A2 * num_InvExpon (32, state->Z);
   state->Z &= MASK32;
   return state->Z;

}

/*-------------------------------------------------------------------------*/

static double InvImpl2b32_U01 (void *vpar, void *vsta)
{
   return InvImpl2b32_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrInvImpl2b (void *vsta)
{
   InvImpl2b_state *state = vsta;
   printf (" Z = %1lu\n", state->Z);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvImpl2b (int e, unsigned long a1, unsigned long a2,
   unsigned long z0)
{
   unif01_Gen *gen;
   InvImpl2b_param *param;
   InvImpl2b_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned long M;

   M = num_TwoExp[e];

   if ((((a1 >= M) || (a2 >= M) || (z0 >= M)) && (e < 32)) ||
      ((a1 % 2) != 1) || ((a2 % 2) != 1) || (e < 3) || (e > 32))
      util_Error ("uinv_CreateInvImpl2b:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvImpl2b_param));
   state = util_Malloc (sizeof (InvImpl2b_state));

   strcpy (name, "uinv_CreateInvImpl2b:");
   addstr_Long (name, "   e = ", e);
   addstr_Ulong (name, ",   a1 = ", a1);
   addstr_Ulong (name, ",   a2 = ", a2);
   addstr_Ulong (name, ",   z0 = ", z0);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Shift = 32 - e;
   param->M = M;
   param->A1 = a1;
   param->A2 = a2;
   state->Z = z0;
   param->Mask = M - 1;

   if (e == 32) {
      gen->GetBits = &InvImpl2b32_Bits;
      gen->GetU01 = &InvImpl2b32_U01;

   } else if (e == 31) {
      gen->GetBits = &InvImpl2b31_Bits;
      gen->GetU01 = &InvImpl2b31_U01;

   } else {                       /* e <= 30 */
      gen->GetBits = &InvImpl2b_Bits;
      gen->GetU01 = &InvImpl2b_U01;
   }
   gen->param = param;
   gen->state = state;
   gen->Write = &WrInvImpl2b;
   return gen;
}


/**************************************************************************/
/*
 * Non-linear Generator with "explicit" inversion in the form
 * Z = Xn^(-1)  when Xn # 0 ,  else 0 and
 * with Xn = a*n + c as mentionned in L'Ecuyer (1993). 
 *
 * The values returned by the generator implemented below are
 * in the following set: { 0/m, 1/m, 2/m, ..., (m-1)/m }.
 */

static double InvExpl_U01 (void *vpar, void *vsta)
{
   InvExpl_param *param = vpar;
   InvExpl_state *state = vsta;

   if (state->Z >= param->MmA)
      state->Z -= param->M;
   state->Z += param->A2;

   if (state->Z == 0)
      return 0.0;
   else
      return (num_InvEuclid (param->M, state->Z) * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long InvExpl_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * InvExpl_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrInvExpl (void *vsta)
{
   InvExpl_state *state = vsta;
   printf (" Z = %1ld\n", state->Z);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvExpl (long m, long a, long c)
{
   unif01_Gen *gen;
   InvExpl_param *param;
   InvExpl_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a <= 0) || (c < 0) || (a >= m) || (c >= m) || ((m % 2) == 0))
      util_Error ("uinv_CreateInvExpl:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvExpl_param));
   state = util_Malloc (sizeof (InvExpl_state));

   strcpy (name, "uinv_CreateInvExpl:");
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   c = ", c);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->Z = c;
   param->M = m;
   param->A2 = a;
   param->MmA = m - a;
   param->Norm = 1.0 / m;

   gen->param = param;
   gen->state = state;
   gen->Write = &WrInvExpl;
   gen->GetBits = &InvExpl_Bits;
   gen->GetU01 = &InvExpl_U01;

   return gen;
}


/**************************************************************************/

static unsigned long InvExpl2a_Bits (void *vpar, void *vsta)
{
   InvExpl2a_param *param = vpar;
   InvExpl2a_state *state = vsta;

   if (param->E < 31) {
      state->Z += param->A2;
      state->Z &= param->Mask;
      if (state->Z == 0)
         return 0;
      else
         return num_InvExpon (param->E, state->Z) << param->Shift;

   } else if (param->E == 31) {
      state->Z += param->A2;
      state->Z &= MASK31;
      if (state->Z == 0)
         return 0;
      else
         return num_InvExpon (31, state->Z) << 1;

   } else {  /* param->E == 32 */ 
      state->Z += param->A2;
      state->Z &= MASK32;
      if (state->Z == 0)
         return 0;
      else
         return num_InvExpon (32, state->Z);
   }
}

/*-------------------------------------------------------------------------*/

static double InvExpl2a_U01 (void *vpar, void *vsta)
{
   return InvExpl2a_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrInvExpl2a (void *vsta)
{
   InvExpl2a_state *state = vsta;
   printf (" Z = %1lu\n", state->Z);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvExpl2a (int e, long a, long c)
{
   unif01_Gen *gen;
   InvExpl2a_param *param;
   InvExpl2a_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned long M;

   M = num_TwoExp[e];
   if ((a <= 0) || (c <= 0) || (a % 4 != 2) || (c % 2 != 1) ||
      ((((unsigned long) a >= M) || ((unsigned long) c >= M))
         && (e < 32)) || (e < 3) || (e > 32))
      util_Error ("uinv_CreateInvExpl2a:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvExpl2a_param));
   state = util_Malloc (sizeof (InvExpl2a_state));

   state->Z = c;
   param->A2 = a;
   param->Mask = M - 1;
   param->E = e;
   param->Shift = 32 - e;

   strcpy (name, "uinv_CreateInvExpl2a:");
   addstr_Long (name, "   e = ", e);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   c = ", c);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &InvExpl2a_Bits;
   gen->GetU01 = &InvExpl2a_U01;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrInvExpl2a;
   return gen;
}


/**************************************************************************/

static unsigned long InvExpl2b_Bits (void *vpar, void *vsta)
{
   InvExpl2b_param *param = vpar;
   InvExpl2b_state *state = vsta;
   unsigned long Z;

   state->N++;
   if (param->E < 31) {
      Z = state->N * num_InvExpon (param->E,
                      (param->A1 + param->A2 * state->N) & param->Mask);
#ifndef IS_ULONG32
      Z &= param->Mask;
#endif
      return Z << param->Shift;

   } else if (param->E == 31) {
      Z = state->N * num_InvExpon (31,
                     (param->A1 + param->A2 * state->N) & MASK31);
#ifndef IS_ULONG32
      Z &= MASK31;
#endif
      return Z << 1;

   } else {     /* E == 32 */
      Z = state->N * num_InvExpon (32, param->A1 + param->A2 * state->N);
#ifndef IS_ULONG32
      Z &= MASK32;
#endif
      return Z;
   } 
}

/*-------------------------------------------------------------------------*/

static double InvExpl2b_U01 (void *vpar, void *vsta)
{
   return InvExpl2b_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrInvExpl2b (void *vsta)
{
   InvExpl2b_state *state = vsta;
   printf (" N = %1lu\n", state->N);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvExpl2b (int e, long a, long c)
{
   unif01_Gen *gen;
   InvExpl2b_param *param;
   InvExpl2b_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned long M;

   M = num_TwoExp[e];

   if ((a <= 0) || (c <= 0) || (a % 4 != 2) || (c % 2 != 1) ||
      ((((unsigned long) a >= M) || ((unsigned long) c >= M))
         && (e < 32)) || (e < 3) || (e > 32))
      util_Error ("uinv_CreateInvExpl2b:   parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvExpl2b_param));
   state = util_Malloc (sizeof (InvExpl2b_state));

   strcpy (name, "uinv_CreateInvExpl2b:");
   addstr_Long (name, "   e = ", e);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   c = ", c);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->A1 = c;
   param->A2 = a;
   param->E = e;
   param->Shift = 32 - e;
   param->Mask = M - 1;
   state->N = 0;

   gen->GetBits = &InvExpl2b_Bits;
   gen->GetU01 = &InvExpl2b_U01;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrInvExpl2b;
   return gen;
}


/**************************************************************************/
/*
 * INVERSIVE Multiple recursive generator (MRG):   Z = Xn^(-1)
 *
 * The values returned by the generator InvMRG are
 * in the following set: { 1/(m+1), 2/(m+1), ... , m/(m+1) }.
 */

static double InvMRG_U01 (void *vpar, void *vsta)
{
   InvMRG_param *param = vpar;
   InvMRG_state *state = vsta;
   int i;
   long h, t, p = 0;

   for (i = state->Order; i > 0; i--) {
      if (param->A[i] != 0) {
         h = state->S[i] / param->Q[i];
         t = labs (param->A[i]) * (state->S[i] - h * param->Q[i]) -
              h * param->R[i];
         if (t < 0)
            t += param->M;
         if (param->A[i] < 0)
            p -= t;
         else
            p += (t - param->M);
         if (p < 0)
            p += param->M;
      }
      if (i > 1)
         state->S[i] = state->S[i - 1];
      else
         state->S[i] = p;
   }
   if (p == 0)
      return param->M * param->Norm;
   else
      return num_InvEuclid (param->M, p) * param->Norm;
}

/*-------------------------------------------------------------------------*/

static unsigned long InvMRG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * InvMRG_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrInvMRG (void *vsta)
{
   InvMRG_state *state = vsta;
   int i;

   for (i = 1; i <= state->Order; i++) {
      printf ("   S[%1d] = %10ld  ", i, state->S[i]);
      if (i % 3 == 0)
         printf ("\n");
   }
   if ((state->Order % 3) != 0)
      printf ("\n");
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvMRG (long m, int k, long A[], long S[])
{
   unif01_Gen *gen;
   InvMRG_param *param;
   InvMRG_state *state;
   size_t leng;
   char name[LEN + 1];
   int i, n;
   long *aa, *ss, *qq, *rr;

   if ((k < 2) || (m < 2) || (m % 2 == 0))
      util_Error ("uinv_CreateInvMRG:   Invalid parameter *");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvMRG_param));
   state = util_Malloc (sizeof (InvMRG_state));

   strcpy (name, "uinv_CreateInvMRG:");
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   k = ", k);
   addstr_ArrayLong (name, ",   A = ", k, A);
   addstr_ArrayLong (name, ",   S = ", k, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   aa = util_Calloc ((size_t) k + 1, sizeof (long));
   rr = util_Calloc ((size_t) k + 1, sizeof (long));
   qq = util_Calloc ((size_t) k + 1, sizeof (long));
   ss = util_Calloc ((size_t) k + 1, sizeof (long));

   n = 0;
   for (i = 1; i <= k; i++) {
      aa[i] = A[i - 1];
      ss[i] = S[i - 1];
      if ((labs (aa[i]) >= m) || (ss[i] >= m) || (ss[i] < 0))
         util_Error ("uinv_CreateInvMRG:   Invalid parameter **");
      if (aa[i] != 0) {
         rr[i] = m % labs (aa[i]);
         qq[i] = m / labs (aa[i]);
         if (rr[i] > qq[i])
            util_Error ("uinv_CreateInvMRG:   Invalid parameter ***");
      }
      if (ss[i] != 0)
         n++;
   }
   if (n == 0)
      util_Error ("uinv_CreateInvMRG:   Invalid parameter ****");

   param->M = m;
   param->Norm = 1.0 / (m + 1.0);
   param->A = aa;
   param->R = rr;
   param->Q = qq;
   state->Order = k;
   state->S = ss;

   gen->param = param;
   gen->state = state;
   gen->GetBits = &InvMRG_Bits;
   gen->GetU01 = &InvMRG_U01;
   gen->Write = &WrInvMRG;

   return gen;
}

/*-------------------------------------------------------------------------*/

void uinv_DeleteInvMRG (unif01_Gen * gen)
{
   InvMRG_param *param;
   InvMRG_state *state;

   if (NULL == gen)
      return;
   param = gen->param;
   state = gen->state;
   util_Free (state->S);
   util_Free (param->A);
   util_Free (param->Q);
   util_Free (param->R);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/**************************************************************************/

static double InvMRGFloat_U01 (void *vpar, void *vsta)
/*
 * Generator InvMRG of order k. The implementation uses floating-point
 * arithmetic. Similar to umrg_MultRecKFloat
 */
{
   InvMRGFloat_param *param = vpar;
   InvMRGFloat_state *state = vsta;
   long k;
   int i;
   double p = 0.0;

   p = 0.0;
   for (i = state->Order; i > 0; i--) {
      if (param->A[i] != 0.0)
         p += param->A[i] * state->S[i];
      if (i > 1)
         state->S[i] = state->S[i - 1];
   }
   k = p / param->M;
   if (p >= 0.0)
      p -= k * param->M;               /* p = sum ...  % m */
   else
      p += (1 - k) * param->M;

   k = state->S[1] = p;
   if (k == 0)
      return (param->M * param->Norm);
   else
      return (num_InvEuclid ((long) param->M, k) * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long InvMRGFloat_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * InvMRGFloat_U01 (vpar, vsta));
}

/*------------------------------------------------------------------------*/

static void WrInvMRGFloat (void *vsta)
{
   InvMRGFloat_state *state = vsta;
   int i;

   for (i = 1; i <= state->Order; i++) {
      printf ("   S[%1d] = %10ld  ", i, (long) state->S[i]);
      if ((i % 3) == 0)
         printf ("\n");
   }
   if (state->Order % 3 != 0)
      printf ("\n");
}

/*------------------------------------------------------------------------*/

unif01_Gen *uinv_CreateInvMRGFloat (long m, int k, long A[], long S[])
{
   unif01_Gen *gen;
   InvMRGFloat_param *param;
   InvMRGFloat_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;
   unsigned int n;
   double pr1, pr2;
   double *ar, *sr;

   if ((k < 2) || (m < 2) || (m % 2 == 0))
      util_Error
         ("uinv.CreateInvMRGFloat:   k or m invalid");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (InvMRGFloat_param));
   state = util_Malloc (sizeof (InvMRGFloat_state));

   ar = util_Calloc ((size_t) k + 1, sizeof (double));
   sr = util_Calloc ((size_t) k + 1, sizeof (double));

   n = 0;
   pr2 = pr1 = 0.0;
   for (i = 1; i <= k; i++) {
      ar[i] = A[i - 1];
      sr[i] = S[i - 1];
      if ((A[i - 1] >= m) || (-A[i - 1] >= m))
         util_Error ("uinv.CreateInvMRGFloat:   |a_i| >= m");
      else if ((S[i - 1] >= m) || (S[i - 1] < 0))
         util_Error ("uinv.CreateInvMRGFloat:    S_i >= m   or   S_i < 0");
      else if (A[i - 1] < 0)
         pr2 -= ar[i];
      else
         pr1 += ar[i];

      if (S[i - 1] != 0)
         n++;
   }
   if (n == 0)
      util_Error ("uinv.CreateInvMRGFloat:   all S[i] = 0");
   if ((pr1 * m >= num_TwoExp[53]) || (pr2 * m >= num_TwoExp[53]))
      util_Error ("uinv.CreateInvMRGFloat:   invalid a_i");

   strcpy (name, "uinv_CreateInvMRGFloat:");
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   k = ", k);
   addstr_ArrayLong (name, ",   A = ", k, A);
   addstr_ArrayLong (name, ",   S = ", k, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->M = m;
   param->Norm = 1.0 / (m + 1.0);
   param->A = ar;
   state->Order = k;
   state->S = sr;

   gen->param = param;
   gen->state = state;
   gen->GetBits = &InvMRGFloat_Bits;
   gen->GetU01 = &InvMRGFloat_U01;
   gen->Write = &WrInvMRGFloat;

   return gen;
}

/*------------------------------------------------------------------------*/

void uinv_DeleteInvMRGFloat (unif01_Gen * gen)
{
   InvMRGFloat_param *param;
   InvMRGFloat_state *state;

   if (NULL == gen)
      return;
   param = gen->param;
   state = gen->state;
   util_Free (state->S);
   util_Free (param->A);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/***************************************************************************/

void uinv_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}

