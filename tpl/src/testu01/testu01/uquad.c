/*************************************************************************\
 *
 * Package:        TestU01
 * File:           uquad.c
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

#include "uquad.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>




/*============================= Constants =================================*/


#if LONG_MAX <= 2147483647L
#define  MASK64  0
#define  RacLONG_MAX  46340           /* Trunc (sqrt (LONG_MAX)) */
#else
#define  MASK64  0xffffffffffffffffUL /* Mask 64 bits */
#define  RacLONG_MAX  3037000499
#endif

#define  MASK32  0xffffffffUL        /* Mask 32 bits */
#define  LEN 200                     /* Max length of strings */




/*============================== Types ===================================*/

/* Different styles of generator of the same type: */

typedef enum {
   XXQ, SSQ, SMQ, SLQ, MSQ, MMQ,
   MLQ, LSQ, LMQ, LLQ
} GenStyle1;

typedef enum {
   Qu2, Qu2e32
} GenStyle2;

/*-------------------------------------------------------------------------*/

typedef struct {
   long C,                    /* Constant */
      A, B,                   /* Multipliers */
      Qa, Ra,                 /* Qa = M / A,  Ra = M % A */
      Qb, Rb,                 /* Qb = M / B,  Rb = M % B */
      M;                      /* Modulus */
   double Norm;
} Quad_param;

typedef struct {
   long S;
   GenStyle1 style;
} Quad_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long A, B, C, Mask, Shift;
   lebool Flag;
} Quad2_param;

typedef struct {
   unsigned long S;
   GenStyle2 style;
} Quad2_state;





/*============================= Functions ===============================*/

static double LLQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used in the general case.  " Very slow "
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2;

   W1 = num_MultModL (state->S, state->S, 0L, param->M);
   W1 = num_MultModL (param->A, W1, 0L, param->M);
   W2 = num_MultModL (param->B, state->S, param->C, param->M);
   state->S = (W1 - param->M) + W2;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long LLQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LLQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double LMQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when  b * (m % b) < m
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2, h;

   W1 = num_MultModL (state->S, state->S, 0L, param->M);
   W1 = num_MultModL (param->A, W1, param->C, param->M);
   h = state->S / param->Qb;
   W2 = param->B * (state->S - h * param->Qb) - h * param->Rb;
   if (W2 < 0)
      state->S = W2 + W1;
   else
      state->S = (W2 - param->M) + W1;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long LMQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LMQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double LSQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when b*(m-1) holds in a long int.
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2;

   W1 = num_MultModL (state->S, state->S, 0L, param->M);
   W1 = num_MultModL (param->A, W1, param->C, param->M);
   W2 = (param->B * state->S) % param->M;
   state->S = (W1 - param->M) + W2;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long LSQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LSQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double MLQuad_U01 (void *vpar, void *vsta)
/* 
 * Implementation used when  a * (m % a) < m
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2, k;

   W1 = num_MultModL (param->B, state->S, param->C, param->M);
   W2 = num_MultModL (state->S, state->S, 0L, param->M);
   k = W2 / param->Qa;
   W2 = param->A * (W2 - k * param->Qa) - k * param->Ra;
   if (W2 < 0)
      state->S = W2 + W1;
   else
      state->S = (W2 - param->M) + W1;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MLQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MLQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double MMQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when  a * (m % a) < m   and when  b * (m % b) < m
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2, k, h;

   W1 = num_MultModL (state->S, state->S, 0L, param->M);
   k = W1 / param->Qa;
   W1 = param->A * (W1 - k * param->Qa) - k * param->Ra;
   if (W1 < 0)
      W1 += param->M;
   h = state->S / param->Qb;
   W2 = param->B * (state->S - h * param->Qb) - h * param->Rb;
   if (W2 < 0)
      state->S = W2 + W1;
   else
      state->S = (W2 - param->M) + W1;
   if (state->S < 0)
      state->S += param->C;
   else
      state->S = (state->S - param->M) + param->C;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MMQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MMQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double MSQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when  a * (m % a) < m  and when  b*(m-1) holds in a
 * long int
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2, k;

   W1 = (param->B * state->S) % param->M;
   W2 = num_MultModL (state->S, state->S, 0L, param->M);
   k = W2 / param->Qa;
   W2 = param->A * (W2 - k * param->Qa) - k * param->Ra;
   if (W2 < 0)
      state->S = W2 + W1;
   else
      state->S = (W2 - param->M) + W1;
   if (state->S < 0)
      state->S += param->C;
   else
      state->S = (state->S - param->M) + param->C;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MSQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MSQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double SLQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when a*(m-1) holds in a long int
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2;

   W2 = num_MultModL (state->S, state->S, 0L, param->M);
   W2 = (param->A * W2) % param->M;
   W1 = num_MultModL (param->B, state->S, param->C, param->M);
   state->S = (W2 - param->M) + W1;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long SLQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SLQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double SMQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when  b * (m % b) < m  and when  a*(m-1) holds in
 * a long int
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2, h;

   W1 = num_MultModL (state->S, state->S, 0L, param->M);
   W1 = (param->A * W1) % param->M;
   h = state->S / param->Qb;
   W2 = param->B * (state->S - h * param->Qb) - h * param->Rb;
   if (W2 < 0)
      state->S = W2 + W1;
   else
      state->S = (W2 - param->M) + W1;
   if (state->S < 0)
      state->S += param->C;
   else
      state->S = (state->S - param->M) + param->C;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long SMQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SMQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double SSQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when a*(m-1) and b*(m-1) hold in a long int
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2;

   W1 = num_MultModL (state->S, state->S, 0L, param->M);
   W1 = (param->A * W1) % param->M;
   W2 = (param->B * state->S) % param->M;
   state->S = (W1 - param->M) + W2;
   if (state->S < 0)
      state->S += param->C;
   else
      state->S = (state->S - param->M) + param->C;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long SSQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SSQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double XXQuad_U01 (void *vpar, void *vsta)
/*
 * Implementation used when all intermediary results hold in a long int.
 * It is the case when  M < Rac2MaxInt.
 */
{
   Quad_param *param = vpar;
   Quad_state *state = vsta;
   long W1, W2;

   W1 = (param->A * ((state->S * state->S) % param->M)) % param->M;
   W2 = ((param->B * state->S) + param->C) % param->M;
   state->S = (W1 - param->M) + W2;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long XXQuad_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * XXQuad_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrQuad (void *vsta)
{
   Quad_state *state = vsta;

   switch (state->style) {
   case XXQ:
      printf (" XXQuad");
      break;
   case SSQ:
      printf (" SSQuad");
      break;
   case SMQ:
      printf (" SMQuad");
      break;
   case SLQ:
      printf (" SLQuad");
      break;
   case MSQ:
      printf (" MSQuad");
      break;
   case MMQ:
      printf (" MMQuad");
      break;
   case MLQ:
      printf (" MLQuad");
      break;
   case LSQ:
      printf (" LSQuad");
      break;
   case LMQ:
      printf (" LMQuad");
      break;
   case LLQ:
      printf (" LLQuad");
      break;
   default:
      util_Error ("WrQuad:   impossible case");
   }
   printf (",   S = %1ld\n", state->S);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * uquad_CreateQuadratic (long m, long a, long b, long c, long s)
{
   unif01_Gen *gen;
   Quad_param *param;
   Quad_state *state;
   size_t leng;
   char name[LEN + 1];
   GenStyle1 style;

   util_Assert ((a >= 0) && (b >= 0) && (c >= 0) && (s >= 0) &&
      (a < m) && (b < m) && (c < m) && (s < m) && (0 < m) &&
      ((a != 0) || (b != 0)) && ((c != 0) || (s != 0)),
      "uquad_CreateQuadratic:   Invalid Parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Quad_param));
   state = util_Malloc (sizeof (Quad_state));

   strcpy (name, "uquad_CreateQuadratic:");
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   b = ", b);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->A = a;
   param->B = b;
   param->C = c;
   param->M = m;
   param->Qa = m / a;
   param->Ra = m % a;
   param->Qb = m / b;
   param->Rb = m % b;
   param->Norm = 1.0 / m;
   state->S = s;

   if (m - 1 <= LONG_MAX / a)
      style = XXQ;
   else if (param->Ra <= param->Qa)
      style = SLQ;
   else
      style = MLQ;

   if (m - 1 <= LONG_MAX / b)
      style += 1;
   else if (param->Rb <= param->Qb)
      style += 2;
   else
      style += 3;

   if (m - 1 <= RacLONG_MAX)
      style = XXQ;
   state->style = style;

   switch (style) {
   case XXQ:
      gen->GetBits = &XXQuad_Bits;
      gen->GetU01  = &XXQuad_U01;
      break;
   case SSQ:
      gen->GetBits = &SSQuad_Bits;
      gen->GetU01  = &SSQuad_U01;
      break;
   case SMQ:
      gen->GetBits = &SMQuad_Bits;
      gen->GetU01  = &SMQuad_U01;
      break;
   case SLQ:
      gen->GetBits = &SLQuad_Bits;
      gen->GetU01  = &SLQuad_U01;
      break;
   case MSQ:
      gen->GetBits = &MSQuad_Bits;
      gen->GetU01  = &MSQuad_U01;
      break;
   case MMQ:
      gen->GetBits = &MMQuad_Bits;
      gen->GetU01  = &MMQuad_U01;
      break;
   case MLQ:
      gen->GetBits = &MLQuad_Bits;
      gen->GetU01  = &MLQuad_U01;
      break;
   case LSQ:
      gen->GetBits = &LSQuad_Bits;
      gen->GetU01  = &LSQuad_U01;
      break;
   case LMQ:
      gen->GetBits = &LMQuad_Bits;
      gen->GetU01  = &LMQuad_U01;
      break;
   case LLQ:
      gen->GetBits = &LLQuad_Bits;
      gen->GetU01  = &LLQuad_U01;
      break;
   default:
      util_Error ("uquad_CreateQuadratic:   impossible case");
   }

   gen->Write = &WrQuad;
   gen->param = param;
   gen->state = state;
   return gen;
}


/**************************************************************************/

static unsigned long Quad2_Bits (void *vpar, void *vsta)
{
   Quad2_param *param = vpar;
   Quad2_state *state = vsta;

   state->S = (param->A * state->S * state->S + param->B * state->S +
               param->C) & param->Mask;
   if (param->Flag)
      return state->S << param->Shift;
   else
      return state->S >> param->Shift;
}

/*-------------------------------------------------------------------------*/

static double Quad2_U01 (void *vpar, void *vsta)
{
   return Quad2_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static unsigned long Quad2e32_Bits (void *vpar, void *vsta)
/*
 * Special case when M = 2^32
 */
{
   Quad2_param *param = vpar;
   Quad2_state *state = vsta;

   state->S = param->A * state->S * state->S + param->B * state->S + param->C;
#ifndef IS_ULONG32
   state->S &= MASK32;
#endif
   return state->S;
}

/*-------------------------------------------------------------------------*/

static double Quad2e32_U01 (void *vpar, void *vsta)
{
   return Quad2e32_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrQuad2 (void *vsta)
{
   Quad2_state *state = vsta;

   if (state->style == Qu2)
      printf (" Quad2");
   else if (state->style == Qu2e32)
      printf (" Quad2e32");
   printf (":   S = %1lu\n", state->S);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *uquad_CreateQuadratic2 (int e, unsigned long a, unsigned long b,
   unsigned long c, unsigned long s)
{
   unif01_Gen *gen;
   Quad2_param *param;
   Quad2_state *state;
   size_t leng;
   char name[LEN + 1];
   double M;

   util_Assert (((a != 0) || (b != 0)) && ((s != 0) || (c != 0)) &&
      (e >= 2), "uquad_CreateQuadratic2:   Invalid Parameter *");
#if ULONG_MAX <= 4294967295UL
   util_Assert (e <= 32, "uquad_CreateQuadratic2:   e > 32");
#else
   util_Assert (e <= 64, "uquad_CreateQuadratic2:   e > 64");
#endif

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Quad2_param));
   state = util_Malloc (sizeof (Quad2_state));

   strcpy (name, "uquad_CreateQuadratic2: ");
   addstr_Uint (name, "   e = ", (unsigned) e);
   addstr_Ulong (name, ",   a = ", a);
   addstr_Ulong (name, ",   b = ", b);
   addstr_Ulong (name, ",   c = ", c);
   addstr_Ulong (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   M = num_TwoExp[e];
   util_Assert ((s < M) && (a < M) && (b < M) && (c < M),
         "uquad_CreateQuadratic2:   Invalid Parameter **");

   param->A = a;
   param->B = b;
   param->C = c;
   if (e == 64)
      param->Mask = MASK64;
   else if (32 == e)
      param->Mask = MASK32;
   else
      param->Mask = M - 1;
   if (e <= 32) {
      param->Shift = 32 - e;
      param->Flag = TRUE;
   } else {
      param->Shift = e - 32;
      param->Flag = FALSE;
   }
   state->S = s;

   if (e == 32) {
      state->style = Qu2e32;
      gen->GetBits = &Quad2e32_Bits;
      gen->GetU01  = &Quad2e32_U01;

   } else {
      state->style = Qu2;
      gen->GetBits = &Quad2_Bits;
      gen->GetU01  = &Quad2_U01;
   }

   gen->Write = &WrQuad2;
   gen->param = param;
   gen->state = state;
   return gen;
}

/**************************************************************************/

void uquad_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
