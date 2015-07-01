/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ulcg.c
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
#include "addstr.h"

#include "ulcg.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>

#ifdef USE_GMP
#include <gmp.h>
#endif


/*============================== constants ================================*/

#define  LEN  300                 /* Max length of strings */

#define DeuxExp32m1	4294967295UL           /* 2^32 - 1 */
#define DeuxExp31m1	2147483647             /* 2^31 - 1 */
#define UnSur2e31m1	4.656612875245797E-10  /* 1 / (2^31 - 1) */
#define DeuxExp53	9007199254740992.0     /* 2^53 */
#define DeuxExp48       281474976710656ULL     /* 2^48 */
#define MASK32 0xffffffffUL                    /* 2^32 - 1 */
#define MASK31 0x7fffffffUL                    /* 2^31 - 1 */




/*================================= Types =================================*/

typedef struct {
   long M, A, C, q, r;
   double Norm;
} LCG_param;

typedef struct {
   long S;
} LCG_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   double M, A, C;
   double Norm;
} LCGFloat_param;

typedef struct {
   double S;
} LCGFloat_state;

/*-------------------------------------------------------------------------*/
#ifdef USE_GMP

typedef struct {
   mpz_t M, A, C;
   mpq_t U;                       /* U = S / M */
   int cflag;                     /* cflag = 0 means C = 0 */
} BigLCG_param;

typedef struct {
   mpz_t S;
} BigLCG_state;

#endif
/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long A, C, Mask, Shift;
} Pow2LCG_param;

typedef struct {
   unsigned long S;
} Pow2LCG_state;

/*-------------------------------------------------------------------------*/
#ifdef USE_GMP

typedef struct {
   mpz_t M, A, C;
   mpq_t U;                       /* U = S / M */
   unsigned long E;               /* M = 2^E */
   int cflag;                     /* cflag = 0 means C = 0 */
} BigPow2LCG_param;

typedef struct {
   mpz_t S;
} BigPow2LCG_state;

#endif
/*-------------------------------------------------------------------------*/

typedef Pow2LCG_param LCG2e31_param;
typedef Pow2LCG_state LCG2e31_state;

typedef Pow2LCG_param LCG2e32_param;
typedef Pow2LCG_state LCG2e32_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long ahi, alo, Dalo;
} LCG2e31m1HD_param;

typedef struct {
   unsigned long S;
} LCG2e31m1HD_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long a, c;
} LCGPayne_param;

typedef struct {
   unsigned long S;
} LCGPayne_state;

/*-------------------------------------------------------------------------*/

#ifdef USE_LONGLONG

typedef struct {
   ulonglong A, C;
   ulonglong Mask, Shift;
} Pow2LCGL_param;

typedef struct {
   ulonglong S;
} Pow2LCGL_state;

typedef Pow2LCGL_param LCG2e48L_param;
typedef Pow2LCGL_state LCG2e48L_state;

#endif

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long M, H, Q, R, mask1, mask2, emq, emr;
   double norm;
} Wu2_param;

typedef struct {
   unsigned long S;
} Wu2_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   long A1, A2, C1, C2, M1, M2, M1m1, q1, q2, r1, r2;
   double Norm;
} CombLEC2_param;

typedef struct {
   long S1, S2;
} CombLEC2_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   double A1, A2, C1, C2, M1, M2, M1m1;
   double Norm;
} CombLEC2Float_param;

typedef struct {
   double S1, S2;
} CombLEC2Float_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   long A1, A2, A3, C1, C2, C3, M1, M2, M3, M1m1, M1mM3;
   long q1, q2, q3, r1, r2, r3;
   double Norm;
} CombLEC3_param;

typedef struct {
   long S1, S2, S3;
} CombLEC3_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   long A1, A2, C1, C2, M1, M2, q1, q2, r1, r2;
   double Norm1, Norm2;
} CombWH2_param;

typedef struct {
   long S1, S2;
} CombWH2_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   double A1, A2, C1, C2, M1, M2;
   double Norm1, Norm2;
} CombWH2Float_param;

typedef struct {
   double S1, S2;
} CombWH2Float_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   long A1, A2, A3, C1, C2, C3, M1, M2, M3, q1, q2, q3, r1, r2, r3;
   double Norm1, Norm2, Norm3;
} CombWH3_param;

typedef struct {
   long S1, S2, S3;
} CombWH3_state;







/*============================ fonctions ===================================*/

static double SmallLCG_U01 (void *vpar, void *vsta)
/*
 * Implementation used when (a*(m-1) + c) holds in a long int
 */
{
   LCG_param *param = vpar;
   LCG_state *state = vsta;

   state->S = (param->A * state->S + param->C) % param->M;
   return (state->S * param->Norm);
}

static double MediumLCG_U01 (void *vpar, void *vsta)
/*
 * Implementation used when a * (m % a) < m and c != 0
 */
{
   LCG_param *param = vpar;
   LCG_state *state = vsta;
   long k;

   k = state->S / param->q;
   state->S = param->A * (state->S - k * param->q) - k * param->r;
   if (state->S < 0)
      state->S += param->C;
   else
      state->S = (state->S - param->M) + param->C;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

static double MediumMLCG_U01 (void *vpar, void *vsta)
/*
 * Implementation used when a * (m % a) < m and c = 0
 */
{
   LCG_param *param = vpar;
   LCG_state *state = vsta;
   long k;

   k = state->S / param->q;
   state->S = param->A * (state->S - k * param->q) - k * param->r;
   if (state->S < 0)
      state->S += param->M;
   return (state->S * param->Norm);
}

static double LargeLCG_U01 (void *vpar, void *vsta)
/*
 * Implementation used when none of the three above (SmallLCG, MediumMLCG,
 * MediumLCG) can be used
 */
{
   LCG_param *param = vpar;
   LCG_state *state = vsta;

   state->S = num_MultModL (param->A, state->S, param->C, param->M);
   return (state->S * param->Norm);
}


static unsigned long SmallLCG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SmallLCG_U01 (vpar, vsta));
}

static unsigned long MediumLCG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumLCG_U01 (vpar, vsta));
}

static unsigned long MediumMLCG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumMLCG_U01 (vpar, vsta));
}

static unsigned long LargeLCG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LargeLCG_U01 (vpar, vsta));
}

static void WrLCG (void *vsta)
{
   LCG_state *state = vsta;
   printf (" s = %ld\n", state->S);
}


unif01_Gen *ulcg_CreateLCG (long m, long a, long c, long s)
{
   unif01_Gen *gen;
   LCG_param *param;
   LCG_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a < 0) || (c < 0) || (s < 0) || (a >= m) || (c >= m) ||
       (s >= m) || (m <= 0))
      util_Error ("ulcg_CreateLCG:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LCG_param));
   state = util_Malloc (sizeof (LCG_state));

   strncpy (name, "ulcg_CreateLCG:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Norm = 1.0 / m;
   param->M = m;
   param->A = a;
   param->C = c;
   state->S = s;

   gen->param = param;
   gen->state = state;
   gen->Write = &WrLCG;

   if (m - 1 <= (LONG_MAX - c) / a) {
      gen->GetBits = &SmallLCG_Bits;
      gen->GetU01 = &SmallLCG_U01;
   } else {
      param->q = m / a;
      param->r = m % a;
      if (param->r <= param->q) {
         if (c != 0) {
            gen->GetBits = &MediumLCG_Bits;
            gen->GetU01 = &MediumLCG_U01;
         } else {
            gen->GetBits = &MediumMLCG_Bits;
            gen->GetU01 = &MediumMLCG_U01;
         }
      } else {
         gen->GetBits = &LargeLCG_Bits;
         gen->GetU01 = &LargeLCG_U01;
      }
   }
   return gen;
}


/**************************************************************************/
#ifdef USE_GMP

static double BigLCG_U01 (void *vpar, void *vsta)
{
   BigLCG_param *param = vpar;
   BigLCG_state *state = vsta;

   mpz_mul (state->S, param->A, state->S);     /* S = A * S */
   if (param->cflag)
      mpz_add (state->S, state->S, param->C);  /* S = S + C */
   mpz_mod (state->S, state->S, param->M);     /* S = S % M */
   mpq_set_num (param->U, state->S);           /* Numerator of U */
   mpq_set_den (param->U, param->M);           /* Denominator of U */
   return mpq_get_d (param->U);                /* U = S /  M */
}

static unsigned long BigLCG_Bits (void *vpar, void *vsta)
{
   return BigLCG_U01 (vpar, vsta) * unif01_NORM32;
}

static void WrBigLCG (void *vsta)
{
   BigLCG_state *state = vsta;
   printf (" s = ");
   mpz_out_str (NULL, 10, state->S); 
   printf ("\n");
}

unif01_Gen *ulcg_CreateBigLCG (char *m, char *a, char *c, char *s)
{
   unif01_Gen *gen;
   BigLCG_param *param;
   BigLCG_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (BigLCG_param));
   state = util_Malloc (sizeof (BigLCG_state));

   strncpy (name, "ulcg_CreateBigLCG:", (size_t) LEN);
   leng = 36 + strlen (name) + strlen (m) + strlen (a) + strlen (c) +
      strlen (s);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncat (gen->name, name, leng);
   strcat (gen->name, "   m = ");
   strcat (gen->name, m);
   strcat (gen->name, ",   a = ");
   strcat (gen->name, a);
   strcat (gen->name, ",   c = ");
   strcat (gen->name, c);
   strcat (gen->name, ",   s = ");
   strcat (gen->name, s);

   mpz_init (state->S);
   if (mpz_set_str (state->S, s, 10))
      util_Error ("ulcg_CreateBigLCG:   s is not a valid decimal number");
   if (mpz_sgn (state->S) < 0)
      util_Error ("ulcg_CreateBigLCG:   s < 0");

   mpz_init (param->M);
   if (mpz_set_str (param->M, m, 10))
      util_Error ("ulcg_CreateBigLCG:   m is not a valid decimal number");

   mpz_init (param->A);
   if (mpz_set_str (param->A, a, 10))
      util_Error ("ulcg_CreateBigLCG:   a is not a valid decimal number");
   if (mpz_sgn (param->A) < 0)
      util_Error ("ulcg_CreateBigLCG:   a < 0");

   /* Checks if c = 0 */
   if (strcmp (c, "0") == 0)
      param->cflag = 0;
   else {
      param->cflag = 1;
      mpz_init (param->C);
      if (mpz_set_str (param->C, c, 10))
         util_Error ("ulcg_CreateBigLCG:   c is not a valid decimal number");
      if (mpz_sgn (param->C) < 0)
         util_Error ("ulcg_CreateBigLCG:   c < 0");
   }

   mpq_init (param->U);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrBigLCG;
   gen->GetBits = &BigLCG_Bits;
   gen->GetU01 = &BigLCG_U01;
   return gen;
}

void ulcg_DeleteBigLCG (unif01_Gen * gen)
{
   BigLCG_param *param;
   BigLCG_state *state;

   if (NULL == gen)
      return;
   param = gen->param;
   state = gen->state;
   mpz_clear (state->S);
   mpz_clear (param->M);
   mpz_clear (param->A);
   if (param->cflag)
      mpz_clear (param->C);
   mpq_clear (param->U);
   unif01_DeleteGen (gen);
}

#endif
/**************************************************************************/

static double LCGFloat_U01 (void *vpar, void *vsta)
{
   LCGFloat_param *param = vpar;
   LCGFloat_state *state = vsta;
   long k;

   state->S = param->A * state->S + param->C;
   k = state->S / param->M;
   state->S -= k * param->M;
   return (state->S * param->Norm);
}

static double LCGFloatNeg_U01 (void *vpar, void *vsta)
/*
 * When a < 0 and c = 0.
 */
{
   LCGFloat_param *param = vpar;
   LCGFloat_state *state = vsta;
   long k;

   state->S = param->A * state->S;
   k = state->S / param->M;
   state->S += (1 - k) * param->M;
   return (state->S * param->Norm);
}

static unsigned long LCGFloat_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LCGFloat_U01 (vpar, vsta));
}

static unsigned long LCGFloatNeg_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LCGFloatNeg_U01 (vpar, vsta));
}


static void WrLCGFloat (void *vsta)
{
   LCGFloat_state *state = vsta;
   printf (" s = %ld\n", (long) state->S);
}

unif01_Gen *ulcg_CreateLCGFloat (long m, long a, long c, long s)
{
   unif01_Gen *gen;
   LCGFloat_param *param;
   LCGFloat_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((c < 0) || ((a < 0) && (c != 0)) || (a >= m) || (c >= m) || (s >= m))
      util_Error ("ulcg_CreateLCGFloat:   Invalid parameter");
   if (((double) a * m + c >= DeuxExp53) || (-a * (double) m >= DeuxExp53))
      util_Error ("ulcg_CreateLCGFloat:   |am| + c >= 2^{53}");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LCGFloat_param));
   state = util_Malloc (sizeof (LCGFloat_state));

   strncpy (name, "ulcg_CreateLCGFloat:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrLCGFloat;

   param->Norm = 1.0 / m;
   param->M = m;
   param->A = a;
   param->C = c;
   state->S = s;

   if (a < 0) {
      gen->GetBits = &LCGFloatNeg_Bits;
      gen->GetU01 = &LCGFloatNeg_U01;
   } else {
      gen->GetBits = &LCGFloat_Bits;
      gen->GetU01 = &LCGFloat_U01;
   }
   return gen;
}


/**************************************************************************/
#ifdef USE_GMP

static double BigPow2LCG_U01 (void *vpar, void *vsta)
{
   BigPow2LCG_param *param = vpar;
   BigPow2LCG_state *state = vsta;

   mpz_mul (state->S, param->A, state->S);     /* S = A * S */
   if (param->cflag)
      mpz_add (state->S, state->S, param->C);  /* S = S + C */
   mpz_tdiv_r_2exp (state->S, state->S, param->E);     /* S = S % 2^E */
   mpq_set_num (param->U, state->S);           /* Numerator of U */
   mpq_set_den (param->U, param->M);           /* Denominator of U */
   return mpq_get_d (param->U);                /* U = S /  M */
}

static unsigned long BigPow2LCG_Bits (void *vpar, void *vsta)
{
   return (BigPow2LCG_U01 (vpar, vsta) * unif01_NORM32);
}

static void WrBigPow2LCG (void *vsta)
{
   BigPow2LCG_state *state = vsta;
   printf (" s = ");
   mpz_out_str (stdout, 10, state->S); 
   printf ("\n");
}

unif01_Gen *ulcg_CreateBigPow2LCG (long e, char *a, char *c, char *s)
{
   unif01_Gen *gen;
   BigPow2LCG_param *param;
   BigPow2LCG_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (e > 0, "ulcg_CreateBigPow2LCG:   e < 1");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (BigPow2LCG_param));
   state = util_Malloc (sizeof (BigPow2LCG_state));

   strncpy (name, "ulcg_CreateBigPow2LCG:   ", (size_t) LEN);
   addstr_Long (name, "e = ", e);
   leng = 30 + strlen (name) + strlen (a) + strlen (c) + strlen (s);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);
   strcat (gen->name, ",   a = ");
   strcat (gen->name, a);
   strcat (gen->name, ",   c = ");
   strcat (gen->name, c);
   strcat (gen->name, ",   s = ");
   strcat (gen->name, s);

   param->E = e;
   mpz_init (param->M);
   mpz_ui_pow_ui (param->M, 2, (unsigned long) e);           /* M = 2^E */

   mpz_init (state->S);
   if (mpz_set_str (state->S, s, 10))
      util_Error ("ulcg_CreateBigPow2LCG:   s is not a valid decimal number");
   if (mpz_sgn (state->S) < 0)
      util_Error ("ulcg_CreateBigPow2LCG:   s < 0");

   mpz_init (param->A);
   if (mpz_set_str (param->A, a, 10))
      util_Error ("ulcg_CreateBigPow2LCG:   a is not a valid decimal number");
   if (mpz_sgn (param->A) < 0)
      util_Error ("ulcg_CreateBigPow2LCG:   a < 0");

   /* Checks if c = 0 */
   if (strcmp (c, "0") == 0)
      param->cflag = 0;
   else {
      param->cflag = 1;
      mpz_init (param->C);
      if (mpz_set_str (param->C, c, 10))
         util_Error ("ulcg_CreateBigPow2LCG:   c is not a valid decimal number");
      if (mpz_sgn (param->C) < 0)
         util_Error ("ulcg_CreateBigPow2LCG:   c < 0");
   }

   mpq_init (param->U);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrBigPow2LCG;
   gen->GetBits = &BigPow2LCG_Bits;
   gen->GetU01 = &BigPow2LCG_U01;
   return gen;
}

void ulcg_DeleteBigPow2LCG (unif01_Gen * gen)
{
   BigPow2LCG_param *param;
   BigPow2LCG_state *state;

   if (NULL == gen)
      return;
   param = gen->param;
   state = gen->state;
   mpz_clear (state->S);
   mpz_clear (param->M);
   mpz_clear (param->A);
   if (param->cflag)
      mpz_clear (param->C);
   mpq_clear (param->U);
   unif01_DeleteGen (gen);
}

#endif
/**************************************************************************/

static unsigned long Pow2LCG_Bits (void *vpar, void *vsta)
{
   Pow2LCG_param *param = vpar;
   Pow2LCG_state *state = vsta;
   state->S = (param->A * state->S + param->C) & param->Mask;
   return state->S << param->Shift;
}

static double Pow2LCG_U01 (void *vpar, void *vsta)
{
   return (Pow2LCG_Bits (vpar, vsta) * unif01_INV32);
}

static void WrPow2LCG (void *vsta)
{
   Pow2LCG_state *state = vsta;
   printf (" s = %1lu\n", state->S);
}

unif01_Gen *ulcg_CreatePow2LCG (int e, long a, long c, long s)
{
   unif01_Gen *gen;
   Pow2LCG_param *param;
   Pow2LCG_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (e <= 31, "ulcg_CreatePow2LCG:   e > 31");
   if (((a <= 0 || c < 0) || s < 0) || e < 0) {
      util_Error ("ulcg_CreatePow2LCG:   parameter < 0");
   }
   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Pow2LCG_param));
   state = util_Malloc (sizeof (Pow2LCG_state));

   strncpy (name, "ulcg_CreatePow2LCG: ", (size_t) LEN);
   addstr_Int (name, "  e = ", e);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Mask = (unsigned long) (num_TwoExp[e] - 1.0);
   param->Shift = 32 - e;
   param->A = a;
   param->C = c;
   state->S = s;

   gen->param = param;
   gen->state = state;
   gen->Write = &WrPow2LCG;
   gen->GetBits = &Pow2LCG_Bits;
   gen->GetU01 = &Pow2LCG_U01;
   return gen;
}


/**************************************************************************/

static unsigned long LCG2e31_Bits (void *vpar, void *vsta)
{
   LCG2e31_param *param = vpar;
   LCG2e31_state *state = vsta;

   state->S = (param->A * state->S + param->C) & MASK31;
   return state->S << 1;
}

static double LCG2e31_U01 (void *vpar, void *vsta)
{
   return (LCG2e31_Bits (vpar, vsta) * unif01_INV32);
}

static void WrLCG2e31 (void *vsta)
{
   LCG2e31_state *state = vsta;
   printf (" s = %1lu\n", state->S);
}

unif01_Gen *ulcg_CreateLCG2e31 (long a, long c, long s)
{
   unif01_Gen *gen;
   LCG2e31_param *param;
   LCG2e31_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a <= 0) || (c < 0) || (s <= 0) || (s >= DeuxExp31m1)
      || (c >= DeuxExp31m1) || (a >= DeuxExp31m1))
      util_Error ("ulcg_CreateLCG2e31:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LCG2e31_param));
   state = util_Malloc (sizeof (LCG2e31_state));

   strncpy (name, "ulcg_CreateLCG2e31: ", (size_t) LEN);
   addstr_Long (name, "  a = ", a);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->A = a;
   param->C = c;
   state->S = s;

   gen->param = param;
   gen->state = state;
   gen->Write = &WrLCG2e31;
   gen->GetBits = &LCG2e31_Bits;
   gen->GetU01 = &LCG2e31_U01;
   return gen;
}


/**************************************************************************/

static unsigned long LCG2e32_Bits (void *vpar, void *vsta)
{
   LCG2e32_param *param = vpar;
   LCG2e32_state *state = vsta;
#ifdef IS_ULONG32
   state->S = (param->A * state->S + param->C);
#else
   state->S = (param->A * state->S + param->C) & MASK32;
#endif
   return state->S;
}

static double LCG2e32_U01 (void *vpar, void *vsta)
{
   return (LCG2e32_Bits (vpar, vsta) * unif01_INV32);
}

static void WrLCG2e32 (void *vsta)
{
   LCG2e32_state *state = vsta;
   printf (" s = %1lu\n", state->S);
}


unif01_Gen *ulcg_CreateLCG2e32 (unsigned long a, unsigned long c,
   unsigned long s)
{
   unif01_Gen *gen;
   LCG2e32_param *param;
   LCG2e32_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a > DeuxExp32m1) || (c > DeuxExp32m1) || (s > DeuxExp32m1))
      util_Error ("ulcg_CreateLCG2e32:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LCG2e32_param));
   state = util_Malloc (sizeof (LCG2e32_state));

   strncpy (name, "ulcg_CreateLCG2e32: ", (size_t) LEN);
   addstr_Ulong (name, "  a = ", a);
   addstr_Ulong (name, ",   c = ", c);
   addstr_Ulong (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->A = a;
   param->C = c;
   state->S = s;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrLCG2e32;
   gen->GetBits = &LCG2e32_Bits;
   gen->GetU01 = &LCG2e32_U01;
   return gen;
}


/**************************************************************************/

static double LCGPayne_U01 (void *vpar, void *vsta)
{
   LCGPayne_param *param = vpar;
   LCGPayne_state *state = vsta;
   unsigned long q;
#ifdef USE_LONGLONG
   ulonglong res;
   res = state->S * (ulonglong) param->a + param->c;
#else
   unsigned long res;
   res = state->S * param->a + param->c;
#endif
   q = (res & MASK31) + (res >> 31);
   if (q >= DeuxExp31m1)
      q -= DeuxExp31m1;
   state->S = q;
   return (q * UnSur2e31m1);
}


static unsigned long LCGPayne_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LCGPayne_U01 (vpar, vsta));
}


static void WrLCGPayne (void *vsta)
{
   LCGPayne_state *state = vsta;
   printf (" s = %1lu\n", state->S);
}


unif01_Gen *ulcg_CreateLCGPayne (long a, long c, long s)
{
   unif01_Gen *gen;
   LCGPayne_param *param;
   LCGPayne_state *state;
   size_t leng;
   char name[LEN + 1];

#ifndef USE_LONGLONG
#ifdef IS_ULONG32
    util_Error ("ulcg_CreateLCGPayne will not work");
#endif
#endif
   if ((a < 1) || (s < 0) || (s >= DeuxExp31m1))
      util_Error ("ulcg_CreateLCGPayne:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LCGPayne_param));
   state = util_Malloc (sizeof (LCGPayne_state));

   strncpy (name, "ulcg_CreateLCGPayne:", (size_t) LEN);
   addstr_Long (name, "   a = ", a);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->a = a;
   param->c = c;
   state->S = s;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrLCGPayne;
   gen->GetBits = &LCGPayne_Bits;
   gen->GetU01 = &LCGPayne_U01;
   return gen;
}


/**************************************************************************/

static double LCG2e31m1HD_U01 (void *vpar, void *vsta)
{
   LCG2e31m1HD_param *param = vpar;
   LCG2e31m1HD_state *state = vsta;
   unsigned long xhi, xlo, mid;

   xhi = state->S >> 16;
   xlo = state->S & 0xffff;
   mid = param->ahi * xlo + param->Dalo * xhi;
   state->S = param->ahi * xhi + mid / 65536 + param->alo * xlo;
   if (state->S > DeuxExp31m1)
      state->S -= DeuxExp31m1;
   state->S += (32768 * (mid & 0xffff));
   if (state->S > DeuxExp31m1)
      state->S -= DeuxExp31m1;
   return (state->S * UnSur2e31m1);
}


static unsigned long LCG2e31m1HD_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LCG2e31m1HD_U01 (vpar, vsta));
}


static void WrLCG2e31m1HD (void *vsta)
{
   LCG2e31m1HD_state *state = vsta;
   printf (" s = %1lu\n", state->S);
}


unif01_Gen *ulcg_CreateLCG2e31m1HD (long a, long s)
{
   unif01_Gen *gen;
   LCG2e31m1HD_param *param;
   LCG2e31m1HD_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a <= 1) || (s <= 0) || (s >= DeuxExp31m1)
      || (a >= 1073741824))       /* must have a < 2^30 */
      util_Error ("ulcg_CreateLCG2e31m1HD:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LCG2e31m1HD_param));
   state = util_Malloc (sizeof (LCG2e31m1HD_state));

   strncpy (name, "ulcg_CreateLCG2e31m1HD: ", (size_t) LEN);
   addstr_Long (name, "  a = ", a);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->ahi = a >> 15;
   param->alo = a & 0x7fff;
   param->Dalo = 2 * param->alo;
   state->S = s;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrLCG2e31m1HD;
   gen->GetBits = &LCG2e31m1HD_Bits;
   gen->GetU01 = &LCG2e31m1HD_U01;
   return gen;
}


/**************************************************************************/
#ifdef USE_LONGLONG

static unsigned long LCG2e48L_Bits (void *vpar, void *vsta)
{
   LCG2e48L_param *param = vpar;
   LCG2e48L_state *state = vsta;

   state->S = (param->A * state->S + param->C) & 0xffffffffffffULL;
   /* We return the bits in a 32 bits int; thus shift right. */
   return (unsigned long) (state->S >> 16);
}

static double LCG2e48L_U01 (void *vpar, void *vsta)
{
   return (LCG2e48L_Bits (vpar, vsta) * unif01_INV32);
}

static void WrLCG2e48L (void *vsta)
{
   LCG2e48L_state *state = vsta;
   printf (" s = %" PRIuLEAST64 "\n\n", state->S);
}

unif01_Gen *ulcg_CreateLCG2e48L (ulonglong a, ulonglong c, ulonglong s)
{
   unif01_Gen *gen;
   LCG2e48L_param *param;
   LCG2e48L_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((s >= DeuxExp48) || (a >= DeuxExp48) || (c >= DeuxExp48))
      util_Error ("ulcg_CreateLCG2e48L:   parameter >= 281474976710656");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LCG2e48L_param));
   state = util_Malloc (sizeof (LCG2e48L_state));

   strncpy (name, "ulcg_CreateLCG2e48L:", (size_t) LEN);
   addstr_ULONG (name, "   a = ", a);
   addstr_ULONG (name, ",   c = ", c);
   addstr_ULONG (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->A = a;
   param->C = c;
   state->S = s;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrLCG2e48L;
   gen->GetBits = &LCG2e48L_Bits;
   gen->GetU01 = &LCG2e48L_U01;
   return gen;
}


/**************************************************************************/

static unsigned long Pow2LCGLB_Bits (void *vpar, void *vsta)
{
   /* e >= 32 */
   Pow2LCGL_param *param = vpar;
   Pow2LCGL_state *state = vsta;
   state->S = (param->A * state->S + param->C) & param->Mask;
   return state->S >> param->Shift;
}

static unsigned long Pow2LCGLA_Bits (void *vpar, void *vsta)
{
   /* e < 32 */
   Pow2LCGL_param *param = vpar;
   Pow2LCGL_state *state = vsta;
   state->S = (param->A * state->S + param->C) & param->Mask;
   return state->S << param->Shift;
}

static double Pow2LCGLB_U01 (void *vpar, void *vsta)
{
   return (Pow2LCGLB_Bits (vpar, vsta) * unif01_INV32);
}

static double Pow2LCGLA_U01 (void *vpar, void *vsta)
{
   return (Pow2LCGLA_Bits (vpar, vsta) * unif01_INV32);
}

static void WrPow2LCGL (void *vsta)
{
   Pow2LCGL_state *state = vsta;
   printf (" s = %1" PRIuLEAST64 "\n", state->S);
}

unif01_Gen *ulcg_CreatePow2LCGL (int e, ulonglong a, ulonglong c, ulonglong s)
{
   unif01_Gen *gen;
   Pow2LCGL_param *param;
   Pow2LCGL_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (e <= 64, "ulcg_CreatePow2LCGL:   e > 64");
   util_Assert (e >  0,  "ulcg_CreatePow2LCGL:   e <= 0");
   util_Assert  (a != 0, "ulcg_CreatePow2LCGL:   a = 0");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Pow2LCGL_param));
   state = util_Malloc (sizeof (Pow2LCGL_state));

   strncpy (name, "ulcg_CreatePow2LCGL: ", (size_t) LEN);
   addstr_Int (name, "  e = ", e);
   addstr_ULONG (name, ",   a = ", a);
   addstr_ULONG (name, ",   c = ", c);
   addstr_ULONG (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (e < 64)
      param->Mask = (1ULL << e) - 1;
   else
      param->Mask = 0xffffffffffffffffULL;
   if (e <= 32) {
      param->Shift = 32 - e;
      gen->GetBits = &Pow2LCGLA_Bits;
      gen->GetU01 = &Pow2LCGLA_U01;
   } else {
      param->Shift = e - 32;
      gen->GetBits = &Pow2LCGLB_Bits;
      gen->GetU01 = &Pow2LCGLB_U01;
   }
   param->A = a;
   param->C = c;
   state->S = s;

   gen->param = param;
   gen->state = state;
   gen->Write = &WrPow2LCGL;
   return gen;
}

#endif


/**************************************************************************/

static double Wu2pp_U01 (void *vpar, void *vsta)
{
   Wu2_param *param = vpar;
   Wu2_state *state = vsta;
   unsigned long k, k1, x0, x1;

   x0 = state->S & param->mask1;
   x1 = state->S >> param->emq;
   k = (x0 << param->Q) + param->H * x1;
   if (k >= param->M)
      k -= param->M;
   x0 = state->S & param->mask2;
   x1 = state->S >> param->emr;
   k1 = (x0 << param->R) + param->H * x1;
   if (k1 >= param->M)
      k1 -= param->M;
   state->S = k + k1;
   if (state->S >= param->M)
      state->S -= param->M;
   return state->S * param->norm;
}

static double Wu2mm_U01 (void *vpar, void *vsta)
{
   Wu2_param *param = vpar;
   Wu2_state *state = vsta;
   unsigned long k, k1, x0, x1;

   x0 = state->S & param->mask1;
   x1 = state->S >> param->emq;
   k = (x0 << param->Q) + param->H * x1;
   if (k >= param->M)
      k -= param->M;
   x0 = state->S & param->mask2;
   x1 = state->S >> param->emr;
   k1 = (x0 << param->R) + param->H * x1;
   if (k1 >= param->M)
      k1 -= param->M;
   k += k1;
   if (k < param->M)
      state->S = param->M - k;
   else
      state->S = 2 * param->M - k;
   return state->S * param->norm;
}

static double Wu2pm_U01 (void *vpar, void *vsta)
{
   Wu2_param *param = vpar;
   Wu2_state *state = vsta;
   unsigned long k, k1, x0, x1;

   x0 = state->S & param->mask1;
   x1 = state->S >> param->emq;
   k = (x0 << param->Q) + param->H * x1;
   if (k >= param->M)
      k -= param->M;
   x0 = state->S & param->mask2;
   x1 = state->S >> param->emr;
   k1 = (x0 << param->R) + param->H * x1;
   if (k1 >= param->M)
      k1 -= param->M;
   if (k >= k1)
      state->S = k - k1;
   else
      state->S = param->M + k - k1;
   return state->S * param->norm;
}

static unsigned long Wu2pp_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Wu2pp_U01 (vpar, vsta));
}

static unsigned long Wu2mm_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Wu2mm_U01 (vpar, vsta));
}

static unsigned long Wu2pm_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Wu2pm_U01 (vpar, vsta));
}


static void WrWu2 (void *vsta)
{
   Wu2_state *state = vsta;
   printf (" s = %lu\n", state->S);
}

unif01_Gen *ulcg_CreateLCGWu2 (long m, char o1, unsigned int q, char o2,
   unsigned int r, long s)
{
   unif01_Gen *gen;
   Wu2_param *param;
   Wu2_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int E = 1;
   double w;

   util_Assert (s < m, "ulcg_CreateLCGWu2:   s >= m");
   util_Assert ((unsigned long) m <= MASK31, "ulcg_CreateLCGWu2:   m > 2^31 - 1");
   util_Assert ((o1 == '+') || (o1 == '-'),
      "ulcg_CreateLCGWu2:   o1 must be '+' or '-'");
   util_Assert ((o2 == '+') || (o2 == '-'),
      "ulcg_CreateLCGWu2:   o2 must be '+' or '-'");

   /* find the smallest power of 2 > m */
   while (num_TwoExp[E] < m)
      E++;
   util_Assert (q <= E, "ulcg_CreateLCGWu2:   q > E");
   util_Assert (r <= E, "ulcg_CreateLCGWu2:   r > E");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Wu2_param));
   state = util_Malloc (sizeof (Wu2_state));

   strncpy (name, "ulcg_CreateLCGWu2: ", (size_t) LEN);
   addstr_Long (name, "  m = ", m);
   addstr_Char (name, ",   o1 = ", o1);
   addstr_Long (name, ",   q = ", (long) q);
   addstr_Char (name, ",   o2 = ", o2);
   addstr_Long (name, ",   r = ", (long) r);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->H = num_TwoExp[E] - m;
   param->R = r;
   param->Q = q;

   util_Assert (param->H < num_TwoExp[q], "ulcg_CreateLCGWu2:   h >= 2^q");
   w = param->H*(num_TwoExp[q] - (param->H + 1)/num_TwoExp[E-q]);
   util_Assert (w < m, "ulcg_CreateLCGWu2:   parameters (q)");
   util_Assert (param->H < num_TwoExp[r], "ulcg_CreateLCGWu2:   h >= 2^r");
   w = param->H*(num_TwoExp[r] - (param->H + 1)/num_TwoExp[E-r]);
   util_Assert (w < m, "ulcg_CreateLCGWu2:   parameters (r)");

   w = num_TwoExp[E] - num_TwoExp[q] + param->H * ((m - 1) >> (E - q));
   util_Assert (w < 2.0 * m, "ulcg_CreateLCGWu2:   parameters (Q)");
   w = num_TwoExp[E] - num_TwoExp[r] + param->H * ((m - 1) >> (E - r));
   util_Assert (w < 2.0 * m, "ulcg_CreateLCGWu2:   parameters (R)");
 
   if (o1 == '-') {
      if (o2 == '-') {
         gen->GetBits = &Wu2mm_Bits;
         gen->GetU01 = &Wu2mm_U01;
      } else {
         param->R = q;
         param->Q = r;
         gen->GetBits = &Wu2pm_Bits;
         gen->GetU01 = &Wu2pm_U01;
      }
   } else {
      if (o2 == '-') {
         gen->GetBits = &Wu2pm_Bits;
         gen->GetU01 = &Wu2pm_U01;
      } else {
         gen->GetBits = &Wu2pp_Bits;
         gen->GetU01 = &Wu2pp_U01;
      }
   }

   param->emq = E - q;
   param->emr = E - r;
   param->mask1 = num_TwoExp[param->emq] - 1.0;
   param->mask2 = num_TwoExp[param->emr] - 1.0;
   param->M = m;
   param->norm = 1.0 / m;
   state->S = s % m;

   gen->param = param;
   gen->state = state;
   gen->Write = &WrWu2;
   return gen;
}

/**************************************************************************/

/**********************************************************************
 *
 * Combined generators LCG following L'Ecuyer's method (1988).
 * The state of each generator is in {1, 2, ..., Mi-1}.
 * The returned values are in {1/M1, 2/M1, ..., (M1-1)/M1}.
 * 
 **********************************************************************/

static double SmallCombLEC2_U01 (void *vpar, void *vsta)
{
   CombLEC2_param *param = vpar;
   CombLEC2_state *state = vsta;
   long z;

   state->S1 = ((param->A1 * state->S1) + param->C1) % param->M1;
   state->S2 = ((param->A2 * state->S2) + param->C2) % param->M2;
   z = state->S1 - state->S2;
   if (z < 1)
      z += param->M1m1;
   return (z * param->Norm);
}


static double MediumMCombLEC2_U01 (void *vpar, void *vsta)
{
   CombLEC2_param *param = vpar;
   CombLEC2_state *state = vsta;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->M1;
   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->M2;
   k = state->S1 - state->S2;
   if (k < 1)
      k += param->M1m1;
   return (k * param->Norm);
}


static double MediumCombLEC2_U01 (void *vpar, void *vsta)
{
   CombLEC2_param *param = vpar;
   CombLEC2_state *state = vsta;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->C1;
   else
      state->S1 = (state->S1 - param->M1) + param->C1;
   if (state->S1 < 0)
      state->S1 += param->M1;

   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->C2;
   else
      state->S2 = (state->S2 - param->M2) + param->C2;
   if (state->S2 < 0)
      state->S2 += param->M2;

   k = state->S1 - state->S2;
   if (k < 1)
      k += param->M1m1;

   return (k * param->Norm);
}


static double LargeCombLEC2_U01 (void *vpar, void *vsta)
{
   CombLEC2_param *param = vpar;
   CombLEC2_state *state = vsta;
   long z;

   state->S1 = num_MultModL (param->A1, state->S1, param->C1, param->M1);
   state->S2 = num_MultModL (param->A2, state->S2, param->C2, param->M2);
   z = state->S1 - state->S2;
   if (z < 1)
      z += param->M1m1;
   return (z * param->Norm);
}

static unsigned long SmallCombLEC2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SmallCombLEC2_U01 (vpar, vsta));
}

static unsigned long MediumCombLEC2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumCombLEC2_U01 (vpar, vsta));
}

static unsigned long MediumMCombLEC2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumMCombLEC2_U01 (vpar, vsta));
}

static unsigned long LargeCombLEC2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LargeCombLEC2_U01 (vpar, vsta));
}

static void WrCombLEC2 (void *vsta)
{
   CombLEC2_state *state = vsta;
   printf (" s1 = %1ld,   s2 = %1ld\n", state->S1, state->S2);
}


unif01_Gen *ulcg_CreateCombLEC2 (long m1, long m2, long a1, long a2,
   long c1, long c2, long s1, long s2)
{
   unif01_Gen *gen;
   CombLEC2_param *param;
   CombLEC2_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a1 < 0) || (c1 < 0) || (s1 < 0) || (a1 >= m1) || (c1 >= m1) ||
      (s1 >= m1) || (a2 < 0) || (c2 < 0) || (s2 < 0) || (a2 >= m2) ||
      (c2 >= m2) || (s2 >= m2) || (m2 > m1) || (m2 <= 0) || (m1 <= 0))
      util_Error ("ulcg_CreateCombLEC2:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CombLEC2_param));
   state = util_Malloc (sizeof (CombLEC2_state));

   strncpy (name, "ulcg_CreateCombLEC2:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   c1 = ", c1);
   addstr_Long (name, ",   s1 = ", s1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   c2 = ", c2);
   addstr_Long (name, ",   s2 = ", s2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrCombLEC2;

   if ((m1 - 1 <= (LONG_MAX - c1) / a1) && (m2 - 1 <= (LONG_MAX - c2) / a2)) {
      gen->GetBits = &SmallCombLEC2_Bits;
      gen->GetU01 = &SmallCombLEC2_U01;
   } else {
      param->q1 = m1 / a1;
      param->r1 = m1 % a1;
      param->q2 = m2 / a2;
      param->r2 = m2 % a2;
      if ((param->r1 <= param->q1) && (param->r2 <= param->q2)) {
         if ((c1 == 0) && (c2 == 0)) {
            gen->GetBits = &MediumMCombLEC2_Bits;
            gen->GetU01 = &MediumMCombLEC2_U01;
         } else {
            gen->GetBits = &MediumCombLEC2_Bits;
            gen->GetU01 = &MediumCombLEC2_U01;
         }
      } else {
         gen->GetBits = &LargeCombLEC2_Bits;
         gen->GetU01 = &LargeCombLEC2_U01;
      }
   }
   param->M1 = m1;
   param->M2 = m2;
   param->A1 = a1;
   param->A2 = a2;
   param->C1 = c1;
   param->C2 = c2;
   state->S1 = s1 % m1;
   state->S2 = s2 % m2;
   param->M1m1 = m1 - 1;
   param->Norm = 1.0 / m1;

   return gen;
}


/**************************************************************************/

static double CombLEC2Float_U01 (void *vpar, void *vsta)
{
   CombLEC2Float_param *param = vpar;
   CombLEC2Float_state *state = vsta;
   double z;
   long k;

   state->S1 = param->A1 * state->S1 + param->C1;
   k = state->S1 / param->M1;
   state->S1 -= k * param->M1;
   state->S2 = param->A2 * state->S2 + param->C2;
   k = state->S2 / param->M2;
   state->S2 -= k * param->M2;
   z = state->S1 - state->S2;
   if (z < 1.0)
      z = z + param->M1m1;
   return (z * param->Norm);
}

static unsigned long CombLEC2Float_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombLEC2Float_U01 (vpar, vsta));
}

static void WrCombLEC2Float (void *vsta)
{
   CombLEC2Float_state *state = vsta;
   printf (" s1 = %ld,   s2 =  %ld\n", (long) state->S1, (long) state->S2);
}

unif01_Gen *ulcg_CreateCombLEC2Float (long m1, long m2, long a1, long a2,
   long c1, long c2, long s1, long s2)
{
   unif01_Gen *gen;
   CombLEC2Float_param *param;
   CombLEC2Float_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a1 < 0) || (c1 < 0) || (s1 < 0) || (a1 >= m1) || (c1 >= m1) ||
      (s1 >= m1) || (a2 < 0) || (c2 < 0) || (s2 < 0) || (a2 >= m2) ||
      (c2 >= m2) || (s2 >= m2) || (m2 > m1))
      util_Error ("ulcg_CreateCombLEC2Float:   Invalid parameter");

   if ((a1 * (m1 - 1.0) + c1) >= num_TwoExp[53])
      util_Error ("ulcg_CreateCombLEC2Float:   a1m1 + c1 - a1 >= 2^{53}");
   if ((a2 * (m2 - 1.0) + c2) >= num_TwoExp[53])
      util_Error ("ulcg_CreateCombLEC2Float:   a2m2 + c2 - a2 >= 2^{53}");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CombLEC2Float_param));
   state = util_Malloc (sizeof (CombLEC2Float_state));

   strncpy (name, "ulcg_CreateCombLEC2Float:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   c1 = ", c1);
   addstr_Long (name, ",   s1 = ", s1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   c2 = ", c2);
   addstr_Long (name, ",   s2 = ", s2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrCombLEC2Float;
   gen->GetBits = &CombLEC2Float_Bits;
   gen->GetU01 = &CombLEC2Float_U01;

   param->M1 = m1;
   param->A1 = a1;
   param->C1 = c1;
   state->S1 = s1 % m1;
   param->M2 = m2;
   param->A2 = a2;
   param->C2 = c2;
   state->S2 = s2 % m2;
   param->M1m1 = m1 - 1;
   param->Norm = 1.0 / m1;

   return gen;
}


/**************************************************************************/

static double SmallCombLEC3_U01 (void *vpar, void *vsta)
{
   CombLEC3_param *param = vpar;
   CombLEC3_state *state = vsta;
   long z;

   state->S1 = ((param->A1 * state->S1) + param->C1) % param->M1;
   state->S2 = ((param->A2 * state->S2) + param->C2) % param->M2;
   state->S3 = ((param->A3 * state->S3) + param->C3) % param->M3;

   z = state->S1 - state->S2;
   if (z > param->M1mM3)
      z -= param->M1m1;
   z += state->S3;
   if (z < 1)
      z += param->M1m1;

   return (z * param->Norm);
}


static double MediumCombLEC3_U01 (void *vpar, void *vsta)
{
   CombLEC3_param *param = vpar;
   CombLEC3_state *state = vsta;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->C1;
   else
      state->S1 = (state->S1 - param->M1) + param->C1;
   if (state->S1 < 0)
      state->S1 += param->M1;

   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->C2;
   else
      state->S2 = (state->S2 - param->M2) + param->C2;
   if (state->S2 < 0)
      state->S2 += param->M2;

   k = state->S3 / param->q3;
   state->S3 = param->A3 * (state->S3 - k * param->q3) - k * param->r3;
   if (state->S3 < 0)
      state->S3 += param->C3;
   else
      state->S3 = (state->S3 - param->M3) + param->C3;
   if (state->S3 < 0)
      state->S3 += param->M3;

   k = state->S1 - state->S2;
   if (k > param->M1mM3)
      k -= param->M1m1;
   k += state->S3;
   if (k < 1)
      k += param->M1m1;

   return (k * param->Norm);
}


static double MediumMCombLEC3_U01 (void *vpar, void *vsta)
{
   CombLEC3_param *param = vpar;
   CombLEC3_state *state = vsta;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->M1;

   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->M2;

   k = state->S3 / param->q3;
   state->S3 = param->A3 * (state->S3 - k * param->q3) - k * param->r3;
   if (state->S3 < 0)
      state->S3 += param->M3;

   k = state->S1 - state->S2;
   if (k > param->M1mM3)
      k -= param->M1m1;
   k += state->S3;
   if (k < 1)
      k += param->M1m1;

   return (k * param->Norm);
}


static double LargeCombLEC3_U01 (void *vpar, void *vsta)
{
   CombLEC3_param *param = vpar;
   CombLEC3_state *state = vsta;
   long z;

   state->S1 = num_MultModL (param->A1, state->S1, param->C1, param->M1);
   state->S2 = num_MultModL (param->A2, state->S2, param->C2, param->M2);
   state->S3 = num_MultModL (param->A3, state->S3, param->C3, param->M3);

   z = state->S1 - state->S2;
   if (z > param->M1mM3)
      z -= param->M1m1;
   z += state->S3;
   if (z < 1)
      z += param->M1m1;

   return (z * param->Norm);
}


static unsigned long SmallCombLEC3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SmallCombLEC3_U01 (vpar, vsta));
}

static unsigned long MediumCombLEC3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumCombLEC3_U01 (vpar, vsta));
}

static unsigned long MediumMCombLEC3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumMCombLEC3_U01 (vpar, vsta));
}

static unsigned long LargeCombLEC3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LargeCombLEC3_U01 (vpar, vsta));
}

static void WrCombLEC3 (void *vsta)
{
   CombLEC3_state *state = vsta;
   printf (" s1 = %1ld,   s2 = %1ld,   s3 = %1ld\n",
      state->S1, state->S2, state->S3);
}


unif01_Gen *ulcg_CreateCombLEC3 (long m1, long m2, long m3, long a1,
   long a2, long a3, long c1, long c2, long c3, long s1, long s2, long s3)
{
   unif01_Gen *gen;
   CombLEC3_param *param;
   CombLEC3_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a1 < 0) || (c1 < 0) || (s1 < 0) || (a1 >= m1) || (c1 >= m1) ||
      (s1 >= m1) || (a2 < 0) || (c2 < 0) || (s2 < 0) || (a2 >= m2) ||
      (c2 >= m2) || (s2 >= m2) || (a3 < 0) || (c3 < 0) || (s3 < 0) ||
      (a3 >= m3) || (c3 >= m3) || (s3 >= m3) || (m2 > m1) || (m3 > m2) ||
      (m3 <= 0) || (m2 <= 0) || (m1 <= 0))
      util_Error ("ulcg_CreateCombLEC3:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CombLEC3_param));
   state = util_Malloc (sizeof (CombLEC3_state));

   strncpy (name, "ulcg_CreateCombLEC3:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   c1 = ", c1);
   addstr_Long (name, ",   s1 = ", s1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   c2 = ", c2);
   addstr_Long (name, ",   s2 = ", s2);
   addstr_Long (name, ",   m3 = ", m3);
   addstr_Long (name, ",   a3 = ", a3);
   addstr_Long (name, ",   c3 = ", c3);
   addstr_Long (name, ",   s3 = ", s3);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrCombLEC3;

     if ((m1 - 1 <= (LONG_MAX - c1) / a1) && (m2 - 1 <= (LONG_MAX - c2) / a2) &&
      (m3 - 1 <= (LONG_MAX - c3) / a3)) {
      gen->GetBits = &SmallCombLEC3_Bits;
      gen->GetU01 = &SmallCombLEC3_U01;
   } else {
      param->q1 = m1 / a1;
      param->r1 = m1 % a1;
      param->q2 = m2 / a2;
      param->r2 = m2 % a2;
      param->q3 = m3 / a3;
      param->r3 = m3 % a3;
      if ((param->r1 <= param->q1) && (param->r2 <= param->q2) &&
         (param->r3 <= param->q3)) {
         if ((c1 == 0) && (c2 == 0) && (c3 == 0)) {
            gen->GetBits = &MediumMCombLEC3_Bits;
            gen->GetU01 = &MediumMCombLEC3_U01;
         } else {
            gen->GetBits = &MediumCombLEC3_Bits;
            gen->GetU01 = &MediumCombLEC3_U01;
         }
      } else {
         gen->GetBits = &LargeCombLEC3_Bits;
         gen->GetU01 = &LargeCombLEC3_U01;
      }
   }
   param->M1 = m1;
   param->M2 = m2;
   param->M3 = m3;
   param->A1 = a1;
   param->A2 = a2;
   param->A3 = a3;
   param->C1 = c1;
   param->C2 = c2;
   param->C3 = c3;
   state->S1 = s1 % m1;
   state->S2 = s2 % m2;
   state->S3 = s3 % m3;
   param->M1mM3 = m1 - m3;
   param->M1m1 = m1 - 1;
   param->Norm = 1.0 / m1;

   return gen;
}


/**************************************************************************/

static double SmallCombWH2_U01 (void *vpar, void *vsta)
{
   CombWH2_param *param = vpar;
   CombWH2_state *state = vsta;
   double Sum;

   state->S1 = (param->A1 * state->S1 + param->C1) % param->M1;
   state->S2 = (param->A2 * state->S2 + param->C2) % param->M2;
   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2;
   if (Sum >= 1.0)
      return Sum - 1.0;
   else
      return Sum;
}


static double MediumMCombWH2_U01 (void *vpar, void *vsta)
{
   CombWH2_param *param = vpar;
   CombWH2_state *state = vsta;
   double Sum;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->M1;

   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->M2;

   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2;
   if (Sum >= 1.0)
      return Sum - 1.0;
   else
      return Sum;
}


static double MediumCombWH2_U01 (void *vpar, void *vsta)
{
   CombWH2_param *param = vpar;
   CombWH2_state *state = vsta;
   double Sum;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->C1;
   else
      state->S1 = (state->S1 - param->M1) + param->C1;
   if (state->S1 < 0)
      state->S1 += param->M1;

   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->C2;
   else
      state->S2 = (state->S2 - param->M2) + param->C2;
   if (state->S2 < 0)
      state->S2 += param->M2;

   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2;
   if (Sum >= 1.0)
      return Sum - 1.0;
   else
      return Sum;
}


static double LargeCombWH2_U01 (void *vpar, void *vsta)
{
   CombWH2_param *param = vpar;
   CombWH2_state *state = vsta;
   double Sum;

   state->S1 = num_MultModL (param->A1, state->S1, param->C1, param->M1);
   state->S2 = num_MultModL (param->A2, state->S2, param->C2, param->M2);
   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2;
   if (Sum >= 1.0)
      return Sum - 1.0;
   else
      return Sum;
}


static unsigned long SmallCombWH2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SmallCombWH2_U01 (vpar, vsta));
}

static unsigned long MediumCombWH2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumCombWH2_U01 (vpar, vsta));
}

static unsigned long MediumMCombWH2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumMCombWH2_U01 (vpar, vsta));
}

static unsigned long LargeCombWH2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LargeCombWH2_U01 (vpar, vsta));
}

static void WrCombWH2 (void *vsta)
{
   CombWH2_state *state = vsta;
   printf (" s1 = %1ld,   s2 = %1ld\n", state->S1, state->S2);
}


unif01_Gen *ulcg_CreateCombWH2 (long m1, long m2, long a1, long a2,
   long c1, long c2, long s1, long s2)
{
   unif01_Gen *gen;
   CombWH2_param *param;
   CombWH2_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a1 < 0) || (c1 < 0) || (s1 < 0) || (a1 >= m1) || (c1 >= m1) ||
      (s1 >= m1) || (a2 < 0) || (c2 < 0) || (s2 < 0) || (a2 >= m2) ||
      (c2 >= m2) || (s2 >= m2) || (m2 <= 0) || (m1 <= 0))
      util_Error ("ulcg_CreateCombWH2:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CombWH2_param));
   state = util_Malloc (sizeof (CombWH2_state));

   strncpy (name, "ulcg_CreateCombWH2:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   c1 = ", c1);
   addstr_Long (name, ",   s1 = ", s1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   c2 = ", c2);
   addstr_Long (name, ",   s2 = ", s2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrCombWH2;

   if ((m1 - 1 <= (LONG_MAX - c1) / a1) && (m2 - 1 <= (LONG_MAX - c2) / a2)) {
      gen->GetBits = &SmallCombWH2_Bits;
      gen->GetU01 = &SmallCombWH2_U01;
   } else {
      param->q1 = m1 / a1;
      param->r1 = m1 % a1;
      param->q2 = m2 / a2;
      param->r2 = m2 % a2;
      if ((param->r1 <= param->q1) && (param->r2 <= param->q2)) {
         if ((c1 == 0) && (c2 == 0)) {
            gen->GetBits = &MediumMCombWH2_Bits;
            gen->GetU01 = &MediumMCombWH2_U01;
         } else {
            gen->GetBits = &MediumCombWH2_Bits;
            gen->GetU01 = &MediumCombWH2_U01;
         }
      } else {
         gen->GetBits = &LargeCombWH2_Bits;
         gen->GetU01 = &LargeCombWH2_U01;
      }
   }
   param->M1 = m1;
   param->M2 = m2;
   param->A1 = a1;
   param->A2 = a2;
   param->C1 = c1;
   param->C2 = c2;
   state->S1 = s1 % m1;
   state->S2 = s2 % m2;
   param->Norm1 = 1.0 / m1;
   param->Norm2 = 1.0 / m2;
   return gen;
}


/**************************************************************************/

static double CombWH2Float_U01 (void *vpar, void *vsta)
{
   CombWH2Float_param *param = vpar;
   CombWH2Float_state *state = vsta;
   double Sum;
   long k;

   state->S1 = param->A1 * state->S1 + param->C1;
   k = state->S1 / param->M1;
   state->S1 -= k * param->M1;
   state->S2 = param->A2 * state->S2 + param->C2;
   k = state->S2 / param->M2;
   state->S2 -= k * param->M2;

   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2;
   if (Sum >= 1.0)
      return Sum - 1.0;
   else
      return Sum;
}

static unsigned long CombWH2Float_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombWH2Float_U01 (vpar, vsta));
}

static void WrCombWH2Float (void *vsta)
{
   CombWH2Float_state *state = vsta;
   printf (" s1 = %ld,   s2 = %ld\n", (long) state->S1, (long) state->S2);
}

unif01_Gen *ulcg_CreateCombWH2Float (long m1, long m2, long a1, long a2,
   long c1, long c2, long s1, long s2)
{
   unif01_Gen *gen;
   CombWH2Float_param *param;
   CombWH2Float_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a1 < 0) || (c1 < 0) || (s1 < 0) || (a1 >= m1) || (c1 >= m1) ||
      (s1 >= m1) || (a2 < 0) || (c2 < 0) || (s2 < 0) || (a2 >= m2) ||
      (c2 >= m2) || (s2 >= m2) || (m2 > m1))
      util_Error ("ulcg_CreateCombWH2Float:   Invalid parameter");

   if ((a1 * (m1 - 1.0) + c1) >= num_TwoExp[53])
      util_Error ("ulcg_CreateCombWH2Float:   a1m1 + c1 - a1 >= 2^{53}");

   if ((a2 * (m2 - 1.0) + c2) >= num_TwoExp[53])
      util_Error ("ulcg_CreateCombWH2Float:   a2m2 + c2 - a2 >= 2^{53}");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CombWH2Float_param));
   state = util_Malloc (sizeof (CombWH2Float_state));

   strncpy (name, "ulcg_CreateCombWH2Float:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   c1 = ", c1);
   addstr_Long (name, ",   s1 = ", s1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   c2 = ", c2);
   addstr_Long (name, ",   s2 = ", s2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrCombWH2Float;

   param->M1 = m1;
   param->A1 = a1;
   param->C1 = c1;
   state->S1 = s1 % m1;
   param->M2 = m2;
   param->A2 = a2;
   param->C2 = c2;
   state->S2 = s2 % m2;
   param->Norm1 = 1.0 / m1;
   param->Norm2 = 1.0 / m2;

   gen->GetBits = &CombWH2Float_Bits;
   gen->GetU01 = &CombWH2Float_U01;
   return gen;
}

/**************************************************************************/

static double SmallCombWH3_U01 (void *vpar, void *vsta)
{
   CombWH3_param *param = vpar;
   CombWH3_state *state = vsta;
   double Sum;

   state->S1 = (param->A1 * state->S1 + param->C1) % param->M1;
   state->S2 = (param->A2 * state->S2 + param->C2) % param->M2;
   state->S3 = (param->A3 * state->S3 + param->C3) % param->M3;

   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2 +
         state->S3 * param->Norm3;
   if (Sum >= 2.0)
      return Sum - 2.0;
   if (Sum >= 1.0)
      return Sum - 1.0;
   return Sum;
}

static double MediumMCombWH3_U01 (void *vpar, void *vsta)
{
   CombWH3_param *param = vpar;
   CombWH3_state *state = vsta;
   double Sum;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->M1;

   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->M2;

   k = state->S3 / param->q3;
   state->S3 = param->A3 * (state->S3 - k * param->q3) - k * param->r3;
   if (state->S3 < 0)
      state->S3 += param->M3;

   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2 +
         state->S3 * param->Norm3;
   if (Sum >= 2.0)
      return Sum - 2.0;
   if (Sum >= 1.0)
      return Sum - 1.0;
   return Sum;
}

static double MediumCombWH3_U01 (void *vpar, void *vsta)
{
   CombWH3_param *param = vpar;
   CombWH3_state *state = vsta;
   double Sum;
   long k;

   k = state->S1 / param->q1;
   state->S1 = param->A1 * (state->S1 - k * param->q1) - k * param->r1;
   if (state->S1 < 0)
      state->S1 += param->C1;
   else
      state->S1 = (state->S1 - param->M1) + param->C1;
   if (state->S1 < 0)
      state->S1 += param->M1;

   k = state->S2 / param->q2;
   state->S2 = param->A2 * (state->S2 - k * param->q2) - k * param->r2;
   if (state->S2 < 0)
      state->S2 += param->C2;
   else
      state->S2 = (state->S2 - param->M2) + param->C2;
   if (state->S2 < 0)
      state->S2 += param->M2;

   k = state->S3 / param->q3;
   state->S3 = param->A3 * (state->S3 - k * param->q3) - k * param->r3;
   if (state->S3 < 0)
      state->S3 += param->C3;
   else
      state->S3 = (state->S3 - param->M3) + param->C3;
   if (state->S3 < 0)
      state->S3 += param->M3;

   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2 +
         state->S3 * param->Norm3;
   if (Sum >= 2.0)
      return Sum - 2.0;
   if (Sum >= 1.0)
      return Sum - 1.0;
   return Sum;
}


static double LargeCombWH3_U01 (void *vpar, void *vsta)
{
   CombWH3_param *param = vpar;
   CombWH3_state *state = vsta;
   double Sum;

   state->S1 = num_MultModL (param->A1, state->S1, param->C1, param->M1);
   state->S2 = num_MultModL (param->A2, state->S2, param->C2, param->M2);
   state->S3 = num_MultModL (param->A3, state->S3, param->C3, param->M3);

   Sum = state->S1 * param->Norm1 + state->S2 * param->Norm2 +
         state->S3 * param->Norm3;
   if (Sum >= 2.0)
      return Sum - 2.0;
   if (Sum >= 1.0)
      return Sum - 1.0;
   return Sum;
}

static unsigned long SmallCombWH3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SmallCombWH3_U01 (vpar, vsta));
}

static unsigned long MediumCombWH3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumCombWH3_U01 (vpar, vsta));
}

static unsigned long MediumMCombWH3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MediumMCombWH3_U01 (vpar, vsta));
}

static unsigned long LargeCombWH3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LargeCombWH3_U01 (vpar, vsta));
}

static void WrCombWH3 (void *vsta)
{
   CombWH3_state *state = vsta;
   printf (" s1 = %1ld,   s2 = %1ld,   s3 = %1ld\n",
      state->S1, state->S2, state->S3);
}


unif01_Gen *ulcg_CreateCombWH3 (long m1, long m2, long m3, long a1,
   long a2, long a3, long c1, long c2, long c3, long s1, long s2, long s3)
{
   unif01_Gen *gen;
   CombWH3_param *param;
   CombWH3_state *state;
   size_t leng;
   char name[LEN + 1];

   if ((a1 < 0) || (c1 < 0) || (s1 < 0) || (a1 >= m1) || (c1 >= m1) ||
      (s1 >= m1) || (a2 < 0) || (c2 < 0) || (s2 < 0) || (a2 >= m2) ||
      (c2 >= m2) || (s2 >= m2) || (a3 < 0) || (c3 < 0) || (s3 < 0) ||
      (a3 >= m3) || (c3 >= m3) || (s3 >= m3) ||
      (m2 <= 0) || (m1 <= 0) || (m3 <= 0))
      util_Error ("ulcg_CreateCombWH3:   Invalid parameter");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CombWH3_param));
   state = util_Malloc (sizeof (CombWH3_state));

   strncpy (name, "ulcg_CreateCombWH3:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   c1 = ", c1);
   addstr_Long (name, ",   s1 = ", s1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   c2 = ", c2);
   addstr_Long (name, ",   s2 = ", s2);
   addstr_Long (name, ",   m3 = ", m3);
   addstr_Long (name, ",   a3 = ", a3);
   addstr_Long (name, ",   c3 = ", c3);
   addstr_Long (name, ",   s3 = ", s3);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param = param;
   gen->state = state;
   gen->Write = &WrCombWH3;

   if ((m1 - 1 <= (LONG_MAX - c1) / a1) && (m2 - 1 <= (LONG_MAX - c2) / a2) &&
      (m3 - 1 <= (LONG_MAX - c3) / a3)) {
      gen->GetBits = &SmallCombWH3_Bits;
      gen->GetU01 = &SmallCombWH3_U01;
   } else {
      param->q1 = m1 / a1;
      param->r1 = m1 % a1;
      param->q2 = m2 / a2;
      param->r2 = m2 % a2;
      param->q3 = m3 / a3;
      param->r3 = m3 % a3;
      if ((param->r1 <= param->q1) && (param->r2 <= param->q2) &&
         (param->r3 <= param->q3)) {
         if ((c1 == 0) && (c2 == 0) && (c3 == 0)) {
            gen->GetBits = &MediumMCombWH3_Bits;
            gen->GetU01 = &MediumMCombWH3_U01;
         } else {
            gen->GetBits = &MediumCombWH3_Bits;
            gen->GetU01 = &MediumCombWH3_U01;
         }
      } else {
         gen->GetBits = &LargeCombWH3_Bits;
         gen->GetU01 = &LargeCombWH3_U01;
      }
   }
   param->M1 = m1;
   param->M2 = m2;
   param->M3 = m3;
   param->A1 = a1;
   param->A2 = a2;
   param->A3 = a3;
   param->C1 = c1;
   param->C2 = c2;
   param->C3 = c3;
   state->S1 = s1 % m1;
   state->S2 = s2 % m2;
   state->S3 = s3 % m3;
   param->Norm1 = 1.0 / m1;
   param->Norm2 = 1.0 / m2;
   param->Norm3 = 1.0 / m3;
   return gen;
}


/**************************************************************************/

void ulcg_DeleteGen (unif01_Gen * gen)
{
   unif01_DeleteGen (gen);
}
