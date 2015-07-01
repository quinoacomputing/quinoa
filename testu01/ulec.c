/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ulec.c
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
#include "mystr.h"

#include "ulec.h"
#include "utaus.h"
#include "ulcg.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>

#define LEN 255



/*================================ Types ================================*/

typedef struct {
   long S1, S2;
} CombLec88_state;


typedef struct {
   double S1, S2;
} CombLec88Float_state;


typedef struct {
   long S1, S2, S3, S4;
} CLCG4_state;


typedef struct {
   double S1, S2, S3, S4;
} CLCG4Float_state;


typedef struct {
   long S1, S2, S3, S4, S5;
} MRG93_state;


typedef struct {
   long s10, s11, s12, s20, s21, s22;
} CombMRG96_state;


typedef struct {
   double s10, s11, s12, s20, s21, s22;
} CombMRG96Float_state;


typedef struct {
   double S10, S11, S12, S20, S21, S22;
} MRG32k3_state;


typedef struct {
   long S10, S11, S12, S20, S21, S22;
} MRG32k3_L_state;


typedef struct {
   double S10, S11, S12, S13, S14, S20, S21, S22, S23, S24;
} MRG32k5_state;


#ifdef USE_LONGLONG
typedef struct {
   longlong S10, S11, S12, S20, S21, S22;
} MRG63k3_state;


typedef struct {
   ulonglong y1, y2, y3, y4, y5;
} lfsr258_state;
#endif


typedef struct {
   unsigned int z1, z2, z3;
} lfsr88_state;


typedef struct {
   unsigned int z1, z2, z3, z4;
} lfsr113_state;


typedef struct {
   unsigned long x10, x11, x12, x20, x21, x22;
} MRG31k3p_state;






/*============================== Functions ================================*/

#define  UnSurM1  4.65661305739176919E-10

static double CombLec88_U01 (void *junk, void *vsta)
{
   CombLec88_state *state = vsta;
   long k;
#if LONG_MAX >= 87385394472108L
   state->S1 = (40014 * state->S1) % 2147483563;
   state->S2 = (40692 * state->S2) % 2147483399;
#else
   k = state->S1 / 53668;
   state->S1 = 40014 * (state->S1 - k * 53668) - k * 12211;
   if (state->S1 < 0)
      state->S1 += 2147483563;
   k = state->S2 / 52774;
   state->S2 = 40692 * (state->S2 - k * 52774) - k * 3791;
   if (state->S2 < 0)
      state->S2 += 2147483399;
#endif
   k = state->S1 - state->S2;
   if (k < 1)
      k += 2147483562;
   return (k * UnSurM1);
}

/*-------------------------------------------------------------------------*/

static unsigned long CombLec88_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombLec88_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrCombLec88 (void *vsta)
{
   CombLec88_state *state = vsta;
   printf (" S1 = %1ld", state->S1);
   printf (",   S2 = %1ld\n\n", state->S2);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * ulec_CreateCombLec88 (long S1, long S2)
{
   unif01_Gen *gen;
   CombLec88_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert ((0 < S1) && (S1 < 2147483563),
      "ulec_CreateCombLec88:   S1 must be in [1, 2147483562]");
   util_Assert ((0 < S2) && (S2 < 2147483399),
      "ulec_CreateCombLec88:   S2 must be in [1, 2147483399]");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CombLec88_state));

   strncpy (name, "ulec_CreateCombLec88:", (size_t) LEN);
   addstr_Long (name, "   S1 = ", S1);
   addstr_Long (name, ",   S2 = ", S2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = S1;
   state->S2 = S2;

   gen->GetBits = &CombLec88_Bits;
   gen->GetU01 = &CombLec88_U01;
   gen->Write = &WrCombLec88;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*************************************************************************/

static double CombLec88Float_U01 (void *junk, void *vsta)
{
   CombLec88Float_state *state = vsta;
   long k;
   double Z;

   state->S1 *= 40014.0;
   k = state->S1 / 2147483563.0;
   state->S1 -= k * 2147483563.0;

   state->S2 *= 40692.0;
   k = state->S2 / 2147483399.0;
   state->S2 -= k * 2147483399.0;

   Z = state->S1 - state->S2;
   if (Z < 1.0)
      Z += 2147483562.0;
   return Z * UnSurM1;
}

/*-------------------------------------------------------------------------*/

static unsigned long CombLec88Float_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombLec88Float_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrCombLec88Float (void *vsta)
{
   CombLec88Float_state *state = vsta;
   printf (" S1 = %1ld", (long) state->S1);
   printf (",   S2 = %1ld\n\n", (long) state->S2);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * ulec_CreateCombLec88Float (long S1, long S2)
{
   unif01_Gen *gen;
   CombLec88Float_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert ((0 < S1) && (S1 < 2147483563),
      "ulec_CreateCombLec88Float:   S1 must be in [1, 2147483562]");
   util_Assert ((0 < S2) && (S2 < 2147483399),
      "ulec_CreateCombLec88Float:   S2 must be in [1, 2147483399]");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CombLec88Float_state));

   strncpy (name, "ulec_CreateCombLec88Float:", (size_t) LEN);
   addstr_Long (name, "   S1 = ", S1);
   addstr_Long (name, ",   S2 = ", S2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = S1;
   state->S2 = S2;

   gen->GetBits = &CombLec88Float_Bits;
   gen->GetU01 = &CombLec88Float_U01;
   gen->Write = &WrCombLec88Float;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*************************************************************************/

static double CLCG4_U01 (void *junk, void *vsta)
{
   CLCG4_state *state = vsta;
   double U;

#if LONG_MAX >= 446047364265901L
   state->S1 = (45991 * state->S1) % 2147483647;
#else
   long k;
   k = state->S1 / 46693;
   state->S1 = 45991 * (state->S1 - k * 46693) - k * 25884;
#endif
   if (state->S1 < 0)
      state->S1 += 2147483647;
   U = 4.65661287524579692E-10 * state->S1;

#if LONG_MAX >= 446047364265901L
   state->S2 = (207707 * state->S2) % 2147483543;
#else
   k = state->S2 / 10339;
   state->S2 = 207707 * (state->S2 - k * 10339) - k * 870;
#endif
   if (state->S2 < 0)
      state->S2 += 2147483543;
   U -= 4.65661310075985993E-10 * state->S2;
   if (U < 0.0)
      U += 1.0;

#if LONG_MAX >= 446047364265901L
   state->S3 = (138556 * state->S3) % 2147483423;
#else
   k = state->S3 / 15499;
   state->S3 = 138556 * (state->S3 - k * 15499) - k * 3979;
#endif
   if (state->S3 < 0)
      state->S3 += 2147483423;
   U += 4.65661336096842131E-10 * state->S3;
   if (U > 1.0)
      U -= 1.0;

#if LONG_MAX >= 446047364265901L
   state->S4 = (49689 * state->S4) % 2147483323;
#else
   k = state->S4 / 43218;
   state->S4 = 49689 * (state->S4 - k * 43218) - k * 24121;
#endif
   if (state->S4 < 0)
      state->S4 += 2147483323;
   U -= 4.65661357780891134E-10 * state->S4;
   if (U < 0.0)
      U += 1.0;
   return U;
}

/*-------------------------------------------------------------------------*/

static unsigned long CLCG4_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CLCG4_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrCLCG4 (void *vsta)
{
   CLCG4_state *state = vsta;
   printf (" S1 = %1ld", state->S1);
   printf (",   S2 = %1ld", state->S2);
   printf (",   S3 = %1ld", state->S3);
   printf (",   S4 = %1ld\n\n", state->S4);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * ulec_CreateCLCG4 (long S1, long S2, long S3, long S4)
{
   unif01_Gen *gen;
   CLCG4_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CLCG4_state));

   strncpy (name, "ulec_CreateCLCG4:", (size_t) LEN);
   addstr_Long (name, "   S1 = ", S1);
   addstr_Long (name, ",   S2 = ", S2);
   addstr_Long (name, ",   S3 = ", S3);
   addstr_Long (name, ",   S4 = ", S4);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = S1;
   state->S2 = S2;
   state->S3 = S3;
   state->S4 = S4;

   gen->GetBits = &CLCG4_Bits;
   gen->GetU01 = &CLCG4_U01;
   gen->Write = &WrCLCG4;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

/*******************************************************************/

#undef  UnSurM1

#define  m1r  2147483647.0
#define  m2r  2147483543.0
#define  m3r  2147483423.0
#define  m4r  2147483323.0
#define  UnSurM1  (1.0 / m1r)
#define  UnSurM2  (1.0 / m2r)
#define  UnSurM3  (1.0 / m3r)
#define  UnSurM4  (1.0 / m4r)
#define  a1r  45991.0
#define  a2r  207707.0
#define  a3r  138556.0
#define  a4r  49689.0

static double CLCG4Float_U01 (void *junk, void *vsta)
{
   CLCG4Float_state *state = vsta;
   double U;
   long k;

   state->S1 = a1r * state->S1;
   k = state->S1 * UnSurM1;
   state->S1 -= k * m1r;
   U = state->S1 * UnSurM1;

   state->S2 *= a2r;
   k = state->S2 * UnSurM2;
   state->S2 -= k * m2r;
   U -= state->S2 * UnSurM2;
   if (U < 0.0)
      U += 1.0;

   state->S3 *= a3r;
   k = state->S3 * UnSurM3;
   state->S3 -= k * m3r;
   U += state->S3 * UnSurM3;
   if (U > 1.0)
      U -= 1.0;

   state->S4 *= a4r;
   k = state->S4 * UnSurM4;
   state->S4 -= k * m4r;
   U -= state->S4 * UnSurM4;
   if (U < 0.0)
      U += 1.0;

   return U;
}

/*-------------------------------------------------------------------------*/

static unsigned long CLCG4Float_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CLCG4Float_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrCLCG4Float (void *vsta)
{
   CLCG4Float_state *state = vsta;
   printf (" S1 = %1ld", (long) state->S1);
   printf (",   S2 = %1ld", (long) state->S2);
   printf (",   S3 = %1ld", (long) state->S3);
   printf (",   S4 = %1ld\n\n", (long) state->S4);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * ulec_CreateCLCG4Float (long S1, long S2, long S3, long S4)
{
   unif01_Gen *gen;
   CLCG4Float_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert ((0 <= S1) && (S1 < 2147483647),
      "ulec_CreateCLCG4Float:   S1 must be in [0, 2147483646]");
   util_Assert ((0 <= S2) && (S2 < 2147483543),
      "ulec_CreateCLCG4Float:   S2 must be in [0, 2147483542]");
   util_Assert ((0 <= S3) && (S3 < 2147483423),
      "ulec_CreateCLCG4Float:   S3 must be in [0, 2147483422]");
   util_Assert ((0 <= S4) && (S4 < 2147483323),
      "ulec_CreateCLCG4Float:   S4 must be in [0, 2147483322]");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CLCG4Float_state));

   strncpy (name, "ulec_CreateCLCG4Float:", (size_t) LEN);
   addstr_Long (name, "   S1 = ", S1);
   addstr_Long (name, ",   S2 = ", S2);
   addstr_Long (name, ",   S3 = ", S3);
   addstr_Long (name, ",   S4 = ", S4);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = S1;
   state->S2 = S2;
   state->S3 = S3;
   state->S4 = S4;

   gen->GetBits = &CLCG4Float_Bits;
   gen->GetU01 = &CLCG4Float_U01;
   gen->Write = &WrCLCG4Float;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*******************************************************************/

#undef   norm
#define  norm 4.656612873077393E-10

static double MRG93_U01 (void *junk, void *vsta)
{
   MRG93_state *state = vsta;
   long p1, p5;
#if LONG_MAX < 230584299847627572L
   long h;
#endif

   if (state->S1 == 2147483647)
      state->S1 = 0;

#if LONG_MAX >= 230584299847627572L
   p5 = (104480 * state->S5) % 2147483647;
#else
   h = state->S5 / 20554;
   p5 = 104480 * (state->S5 - h * 20554) - h * 1727;
#endif
   state->S5 = state->S4;
   state->S4 = state->S3;
   state->S3 = state->S2;
   state->S2 = state->S1;

#if LONG_MAX >= 230584299847627572L
   p1 = (107374182 * state->S1) % 2147483647;
#else
   h = state->S1 / 20;
   p1 = 107374182 * (state->S1 - h * 20) - h * 7;
#endif
   if (p1 < 0)
      p1 += 2147483647;
   if (p5 > 0)
      p5 -= 2147483647;
   state->S1 = p1 + p5;
   if (state->S1 <= 0)
      state->S1 += 2147483647;
   return (state->S1 * norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG93_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG93_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG93 (void *vsta)
{
   MRG93_state *state = vsta;
   printf (" S1 = %10ld,   S2 = %10ld,   S3 = %10ld,\n",
      state->S1, state->S2, state->S3);
   printf (" S4 = %10ld,   S5 = %10ld\n\n", state->S4, state->S5);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * ulec_CreateMRG93 (long S1, long S2, long S3, long S4, long S5)
{
   unif01_Gen *gen;
   MRG93_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG93_state));

   strncpy (name, "ulec_CreateMRG93:", (size_t) LEN);
   addstr_Long (name, "   S1 = ", S1);
   addstr_Long (name, ",   S2 = ", S2);
   addstr_Long (name, ",   S3 = ", S3);
   addstr_Long (name, ",   S4 = ", S4);
   addstr_Long (name, ",   S5 = ", S5);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S1 = S1;
   state->S2 = S2;
   state->S3 = S3;
   state->S4 = S4;
   state->S5 = S5;

   gen->GetBits = &MRG93_Bits;
   gen->GetU01  = &MRG93_U01;
   gen->Write   = &WrMRG93;
   gen->param   = NULL;
   gen->state   = state;
   return gen;
}

/*******************************************************************/
/*  Combined MRGs  */

#undef  norm
#undef  m2r

#define   m1    2147483647
#define   m1r   2147483647.0
#define   m2    2145483479
#define   m2r   2145483479.0
#define   a12   63308
#define   q12   33921
#define   r12   12979
#define   a13  -183326
#define   q13   11714
#define   r13   2883
#define   a21   86098
#define   q21   24919
#define   r21   7417
#define   a23  -539608
#define   q23   3976
#define   r23   2071
#define   a12r  63308.0
#define   a13r -183326.0
#define   a21r  86098.0
#define   a23r -539608.0
#define   norm 4.656612873077393e-10


static double CombMRG96_U01 (void *junk, void *vsta)
{
   CombMRG96_state *state = vsta;
   long p12, p21;
#if LONG_MAX < 1342441885711174L
   long h, p13, p23;
#endif

#if LONG_MAX >= 1342441885711174L
   p12 = (a12 * state->s11 + a13 * state->s10) % m1;
   if (p12 < 0)
      p12 += m1;
#else
   h = state->s10 / q13;
   p13 = -a13 * (state->s10 - h * q13) - h * r13;
   h = state->s11 / q12;
   p12 = a12 * (state->s11 - h * q12) - h * r12;
   if (p13 < 0)
      p13 += m1;
   if (p12 < 0)
      p12 += m1;
#endif
   state->s10 = state->s11;
   state->s11 = state->s12;
#if LONG_MAX >= 1342441885711174L
   state->s12 = p12;
#else
   state->s12 = p12 - p13;
   if (state->s12 < 0)
      state->s12 += m1;
#endif

#if LONG_MAX >= 1342441885711174L
   p21 = (a21 * state->s22 + a23 * state->s20) % m2;
   if (p21 < 0)
      p21 += m2;
#else
   h = state->s20 / q23;
   p23 = -a23 * (state->s20 - h * q23) - h * r23;
   h = state->s22 / q21;
   p21 = a21 * (state->s22 - h * q21) - h * r21;
   if (p23 < 0)
      p23 += m2;
   if (p21 < 0)
      p21 += m2;
#endif
   state->s20 = state->s21;
   state->s21 = state->s22;

#if LONG_MAX >= 1342441885711174L
   state->s22 = p21;
#else
   state->s22 = p21 - p23;
   if (state->s22 < 0)
      state->s22 += m2;
#endif

   if (state->s12 <= state->s22)
      return ((state->s12 - state->s22 + m1) * norm);
   else
      return ((state->s12 - state->s22) * norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long CombMRG96_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombMRG96_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrCombMRG96 (void *vsta)
{
   CombMRG96_state *state = vsta;
   printf (" s11 = %10ld", state->s10);
   printf (",   s12 = %10ld", state->s11);
   printf (",   s13 = %10ld,\n", state->s12);
   printf (" s21 = %10ld", state->s20);
   printf (",   s22 = %10ld", state->s21);
   printf (",   s23 = %10ld\n\n", state->s22);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateCombMRG96 (long S11, long S12, long S13,
   long S21, long S22, long S23)
{
   unif01_Gen *gen;
   CombMRG96_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CombMRG96_state));

   strncpy (name, "ulec_CreateCombMRG96:", (size_t) LEN);
   addstr_Long (name, " (S11, ..., S23) = (", S11);
   addstr_Long (name, ", ", S12);
   addstr_Long (name, ", ", S13);
   addstr_Long (name, ", ", S21);
   addstr_Long (name, ", ", S22);
   addstr_Long (name, ", ", S23);
   addstr_Char (name, ")", ' ');
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->s10 = S11;
   state->s11 = S12;
   state->s12 = S13;
   state->s20 = S21;
   state->s21 = S22;
   state->s22 = S23;

   gen->param   = NULL;
   gen->state   = state;
   gen->GetBits = &CombMRG96_Bits;
   gen->GetU01  = &CombMRG96_U01;
   gen->Write   = &WrCombMRG96;
   return gen;
}


/***************************************************************************/

#undef  m1
#undef  m2
#undef  a12
#undef  a21

#define m1  2147483647.0
#define m2  2145483479.0
#define a12      63308.0
#define a13n    183326.0
#define a21      86098.0
#define a23n    539608.0

static double CombMRG96Float_U01 (void *junk, void *vsta)
{
   CombMRG96Float_state *state = vsta;
   long k;
   double p1, p2;

   /* Component 1 */
   p1 = a12 * state->s11 - a13n * state->s10;
   k = p1 / m1;
   p1 -= k * m1;
   if (p1 < 0.0)
      p1 += m1;
   state->s10 = state->s11;
   state->s11 = state->s12;
   state->s12 = p1;

   /* Component 2 */
   p2 = a21 * state->s22 - a23n * state->s20;
   k = p2 / m2;
   p2 -= k * m2;
   if (p2 < 0.0)
      p2 += m2;
   state->s20 = state->s21;
   state->s21 = state->s22;
   state->s22 = p2;

   /* Combination */
   if (p1 <= p2)
      return ((p1 - p2 + m1) * norm);
   else
      return ((p1 - p2) * norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long CombMRG96Float_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombMRG96Float_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrCombMRG96Float (void *vsta)
{
   CombMRG96Float_state *state = vsta;
   printf (" s11 = %10ld", (long) state->s10);
   printf (",   s12 = %10ld", (long) state->s11);
   printf (",   s13 = %10ld,\n", (long) state->s12);
   printf (" s21 = %10ld", (long) state->s20);
   printf (",   s22 = %10ld", (long) state->s21);
   printf (",   s23 = %10ld\n\n", (long) state->s22);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateCombMRG96Float (long S11, long S12, long S13,
   long S21, long S22, long S23)
{
   unif01_Gen *gen;
   CombMRG96Float_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (CombMRG96Float_state));

   strncpy (name, "ulec_CreateCombMRG96Float:", (size_t) LEN);
   addstr_Long (name, " (S11, ..., S23) = (", S11);
   addstr_Long (name, ", ", S12);
   addstr_Long (name, ", ", S13);
   addstr_Long (name, ", ", S21);
   addstr_Long (name, ", ", S22);
   addstr_Long (name, ", ", S23);
   addstr_Char (name, ")", ' ');
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->s10 = S11;
   state->s11 = S12;
   state->s12 = S13;
   state->s20 = S21;
   state->s21 = S22;
   state->s22 = S23;

   gen->param   = NULL;
   gen->state   = state;
   gen->GetBits = &CombMRG96Float_Bits;
   gen->GetU01  = &CombMRG96Float_U01;
   gen->Write   = &WrCombMRG96Float;
   return gen;
}


/***************************************************************************/
              
#define FAC24  0.59604644775390625e-7         /* 1 / 2^24 */

static double CombMRG96D_U01 (void *vpar, void *vsta)
{
   double U;
   U = CombMRG96_U01 (vpar, vsta);
   U += CombMRG96_U01 (vpar, vsta) * FAC24;
   if (U < 1.0)
      return U;
   else
      return U - 1.0;
}

/*-------------------------------------------------------------------------*/

static unsigned long CombMRG96D_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (CombMRG96D_U01 (vpar, vsta) * unif01_NORM32);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateCombMRG96D (long S11, long S12, long S13,
   long S21, long S22, long S23)
{
   unif01_Gen *gen;
   size_t j;

   gen = ulec_CreateCombMRG96 (S11, S12, S13, S21, S22, S23);
   j = strlen (gen->name);
   gen->name = util_Realloc (gen->name, (j + 2) * sizeof (char));
   j = strcspn (gen->name, ":");
   mystr_Insert (gen->name, "D", (unsigned int) j);
   gen->GetU01  = &CombMRG96D_U01;
   gen->GetBits = &CombMRG96D_Bits;

   return gen;
}


/***************************************************************************/
#undef FAC24
#define FAC24  0.59604644775390625e-7         /* 1 / 2^24 */

static double CombMRG96FloatD_U01 (void *vpar, void *vsta)
{
   double U;
   U = CombMRG96Float_U01 (vpar, vsta);
   U += CombMRG96Float_U01 (vpar, vsta) * FAC24;
   if (U < 1.0)
      return U;
   else
      return U - 1.0;
}

/*-------------------------------------------------------------------------*/

static unsigned long CombMRG96FloatD_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (CombMRG96FloatD_U01 (vpar, vsta) * unif01_NORM32);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateCombMRG96FloatD (long S11, long S12, long S13,
   long S21, long S22, long S23)
{
   unif01_Gen *gen;
   size_t j;

   gen = ulec_CreateCombMRG96Float (S11, S12, S13, S21, S22, S23);
   j = strlen (gen->name);
   gen->name = util_Realloc (gen->name, (j + 2) * sizeof (char));
   j = strcspn (gen->name, ":");
   mystr_Insert (gen->name, "D", (unsigned int) j);
   gen->GetU01  = &CombMRG96FloatD_U01;
   gen->GetBits = &CombMRG96FloatD_Bits;
   return gen;
}


/***************************************************************************/

#undef m1
#undef m2
#undef a12
#undef a13n
#undef a21
#undef a23n

#define norm1p1 2.328306549295727688e-10
#define norm1   2.32830654983782883e-10
#define norm2   2.328318825240739e-10

#if LONG_MAX >= 9007199254740992L
#define m1   4294967087L
#define m2   4294944443L
#define a12     1403580L
#define a13n     810728L
#define a21      527612L
#define a23n    1370589L
#else
#define m1    1
#define m2    1
#define a12   0
#define a13n  0
#define a21   0
#define a23n  0
#endif



/*-------------------------------------------------------------------------*/

static void WrMRG32k3_L (void *vsta)
{
   MRG32k3_L_state *state = vsta;
   printf (" (s12, s11, s10, s22, s21, s20) = \n");
   printf (" ( %12ld,  %12ld,  %12ld,\n   %12ld,  %12ld,  %12ld )\n",
      state->S10, state->S11, state->S12, state->S20,
      state->S21, state->S22);
}

/*-------------------------------------------------------------------------*/

static unif01_Gen *CreateMRG32k3ii_L (long x10, long x11, long x12, 
   long x20, long x21, long x22, char *nom)
{
   unif01_Gen *gen;
   MRG32k3_L_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG32k3_L_state));

   strcpy (name, nom);
   addstr_Long (name, " (s10, s11, s12, s20, s21, s22) = (", x10);
   addstr_Long (name, ", ", x11);
   addstr_Long (name, ", ", x12);
   addstr_Long (name, ", ", x20);
   addstr_Long (name, ", ", x21);
   addstr_Long (name, ", ", x22);
   addstr_Char (name, "", ')');
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   util_Assert (x10 < m1 && x10 >= 0,
       "ulec_CreateMRG32k3aL:   x10 not in [0, m1)");
   util_Assert (x11 < m1 && x11 >= 0,
       "ulec_CreateMRG32k3aL:   x11 not in [0, m1)");
   util_Assert (x12 < m1 && x12 >= 0,
       "ulec_CreateMRG32k3aL:   x12 not in [0, m1)");
   util_Assert (x20 < m2 && x20 >= 0,
       "ulec_CreateMRG32k3aL:   x20 not in [0, m2)");
   util_Assert (x21 < m2 && x21 >= 0,
       "ulec_CreateMRG32k3aL:   x21 not in [0, m2)");
   util_Assert (x22 < m2 && x22 >= 0,
       "ulec_CreateMRG32k3aL:   x22 not in [0, m2)");
   state->S10 = x10;
   state->S11 = x11;
   state->S12 = x12;
   state->S20 = x20;
   state->S21 = x21;
   state->S22 = x22;

   gen->param = NULL;
   gen->state = state;
   gen->Write = &WrMRG32k3_L;
   return gen;
}

/*-------------------------------------------------------------------------*/

static double MRG32k3a_L_U01 (void *junk, void *vsta)
{
   MRG32k3_L_state *state = vsta;
   long p;

   /* Component 1 */
   p = (a12 * state->S11 - a13n * state->S10) % m1;
   if (p < 0)
      p += m1;
   state->S10 = state->S11;
   state->S11 = state->S12;
   state->S12 = p;

   /* Component 2 */
   p = (a21 * state->S22 - a23n * state->S20) % m2;
   if (p < 0)
      p += m2;
   state->S20 = state->S21;
   state->S21 = state->S22;
   state->S22 = p;

   /* Combination */
   if (state->S12 <= state->S22)
      return (state->S12 - state->S22 + m1) * norm1p1;
   else
      return (state->S12 - state->S22) * norm1p1;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG32k3a_L_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG32k3a_L_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateMRG32k3aL (long x10, long x11, long x12,
   long x20, long x21, long x22)
{
   unif01_Gen *gen;
   util_Assert (sizeof (long) > 6,
      "ulec_CreateMRG32k3aL:   LONG_MAX is too small for this implementation");
   gen = CreateMRG32k3ii_L (x10, x11, x12, x20, x21, x22,
      "ulec_CreateMRG32k3aL:");
   gen->GetBits = &MRG32k3a_L_Bits;
   gen->GetU01  = &MRG32k3a_L_U01;
   return gen;
}


/***************************************************************************/

#undef m1
#undef m2
#undef a12
#undef a13n
#undef a21
#undef a23n

#define norm1p1 2.328306549295727688e-10
#define norm1   2.32830654983782883e-10
#define norm2   2.328318825240739e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0


/*-------------------------------------------------------------------------*/

static void WrMRG32k3 (void *vsta)
{
   MRG32k3_state *state = vsta;
   printf (" (s12, s11, s10, s22, s21, s20) = \n");
   printf (" ( %12.0f,  %12.0f,  %12.0f,\n   %12.0f,  %12.0f,  %12.0f )\n",
      state->S10, state->S11, state->S12, state->S20,
      state->S21, state->S22);
}

/*-------------------------------------------------------------------------*/

static unif01_Gen *CreateMRG32k3ii (double x10, double x11, double x12,
   double x20, double x21, double x22, char *nom)
{
   unif01_Gen *gen;
   MRG32k3_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG32k3_state));

   strcpy (name, nom);
   addstr_Double (name, " (s10, s11, s12, s20, s21, s22) = (", x10);
   addstr_Double (name, ", ", x11);
   addstr_Double (name, ", ", x12);
   addstr_Double (name, ", ", x20);
   addstr_Double (name, ", ", x21);
   addstr_Double (name, ", ", x22);
   addstr_Char (name, "", ')');
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   util_Assert (x10 < m1 && x10 >= 0, "ulec_CreateMRG32k3:   x10 not in [0, m1)");
   util_Assert (x11 < m1 && x11 >= 0, "ulec_CreateMRG32k3:   x11 not in [0, m1)");
   util_Assert (x12 < m1 && x12 >= 0, "ulec_CreateMRG32k3:   x12 not in [0, m1)");
   util_Assert (x20 < m2 && x20 >= 0, "ulec_CreateMRG32k3:   x20 not in [0, m2)");
   util_Assert (x21 < m2 && x21 >= 0, "ulec_CreateMRG32k3:   x21 not in [0, m2)");
   util_Assert (x22 < m2 && x22 >= 0, "ulec_CreateMRG32k3:   x22 not in [0, m2)");
   state->S10 = x10;
   state->S11 = x11;
   state->S12 = x12;
   state->S20 = x20;
   state->S21 = x21;
   state->S22 = x22;

   gen->param = NULL;
   gen->state = state;
   gen->Write = &WrMRG32k3;
   return gen;
}

/*-------------------------------------------------------------------------*/

static double MRG32k3a_U01 (void *junk, void *vsta)
{
   MRG32k3_state *state = vsta;
   long k;
   double p;

   /* Component 1 */
   p = a12 * state->S11 - a13n * state->S10;
   k = p / m1;
   p -= k * m1;
   if (p < 0.0)
      p += m1;
   state->S10 = state->S11;
   state->S11 = state->S12;
   state->S12 = p;

   /* Component 2 */
   p = a21 * state->S22 - a23n * state->S20;
   k = p / m2;
   p -= k * m2;
   if (p < 0.0)
      p += m2;
   state->S20 = state->S21;
   state->S21 = state->S22;
   state->S22 = p;

   /* Combination */
   if (state->S12 <= state->S22)
      return (state->S12 - state->S22 + m1) * norm1p1;
   else
      return (state->S12 - state->S22) * norm1p1;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG32k3a_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG32k3a_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateMRG32k3a (double x10, double x11, double x12,
   double x20, double x21, double x22)
{
   unif01_Gen *gen;
   gen = CreateMRG32k3ii (x10, x11, x12, x20, x21, x22,
      "ulec_CreateMRG32k3a:");
   gen->GetBits = &MRG32k3a_Bits;
   gen->GetU01  = &MRG32k3a_U01;
   return gen;
}


/***************************************************************************/

static double MRG32k3b_U01 (void *junk, void *vsta)
{
   MRG32k3_state *state = vsta;
   long k;
   double p;

   /* Component 1 */
   p = a12 * state->S11 - a13n * state->S10;
   k = p / m1;
   p -= k * m1;
   if (p < 0.0)
      p += m1;
   state->S10 = state->S11;
   state->S11 = state->S12;
   state->S12 = p;

   /* Component 2 */
   p = a21 * state->S22 - a23n * state->S20;
   k = p / m2;
   p -= k * m2;
   if (p < 0.0)
      p += m2;
   state->S20 = state->S21;
   state->S21 = state->S22;
   state->S22 = p;

   /* Combination */
   p = state->S12 * norm1 - state->S22 * norm2;
   if (p < 0.0)
      return (p + 1.0);
   else
      return p;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG32k3b_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG32k3b_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateMRG32k3b (double x10, double x11, double x12,
   double x20, double x21, double x22)
{
   unif01_Gen *gen;
   gen = CreateMRG32k3ii (x10, x11, x12, x20, x21, x22,
      "ulec_CreateMRG32k3b:");
   gen->GetBits = &MRG32k3b_Bits;
   gen->GetU01  = &MRG32k3b_U01;
   return gen;
}


/***************************************************************************/

#undef norm1
#undef norm2
#undef m1
#undef m2
#undef a12
#undef a21
#undef a23

#define norm1 2.3283163396834613e-10
#define norm2 2.3283243092066027e-10
#define m1   4294949027.0
#define m2   4294934327.0
#define a12     1154721.0
#define a14     1739991.0
#define a15n    1108499.0
#define a21     1776413.0
#define a23      865203.0
#define a25n    1641052.0


static double MRG32k5a_U01 (void *junk, void *vsta)
{
   MRG32k5_state *state = vsta;
   long k;
   double p;

   /* Component 1 */
   p = a12 * state->S13 - a15n * state->S10;
   if (p > 0.0)
      p -= a14 * m1;
   p += a14 * state->S11;
   k = p / m1;
   p -= k * m1;
   if (p < 0.0)
      p += m1;
   state->S10 = state->S11;
   state->S11 = state->S12;
   state->S12 = state->S13;
   state->S13 = state->S14;
   state->S14 = p;

   /* Component 2 */
   p = a21 * state->S24 - a25n * state->S20;
   if (p > 0.0)
      p -= a23 * m2;
   p += a23 * state->S22;
   k = p / m2;
   p -= k * m2;
   if (p < 0.0)
      p += m2;
   state->S20 = state->S21;
   state->S21 = state->S22;
   state->S22 = state->S23;
   state->S23 = state->S24;
   state->S24 = p;

   /* Combination */
   if (state->S14 <= state->S24)
      return ((state->S14 - state->S24 + m1) * norm1);
   else
      return ((state->S14 - state->S24) * norm1);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG32k5a_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG32k5a_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double MRG32k5b_U01 (void *junk, void *vsta)
{
   MRG32k5_state *state = vsta;
   long k;
   double p;

   /* Component 1 */
   p = a12 * state->S13 - a15n * state->S10;
   if (p > 0.0)
      p -= a14 * m1;
   p += a14 * state->S11;
   k = p / m1;
   p -= k * m1;
   if (p < 0.0)
      p += m1;
   state->S10 = state->S11;
   state->S11 = state->S12;
   state->S12 = state->S13;
   state->S13 = state->S14;
   state->S14 = p;

   /* Component 2 */
   p = a21 * state->S24 - a25n * state->S20;
   if (p > 0.0)
      p -= a23 * m2;
   p += a23 * state->S22;
   k = p / m2;
   p -= k * m2;
   if (p < 0.0)
      p += m2;
   state->S20 = state->S21;
   state->S21 = state->S22;
   state->S22 = state->S23;
   state->S23 = state->S24;
   state->S24 = p;

   /* Combination */
   p = state->S14 * norm1 - state->S24 * norm2;
   if (p < 0.0)
      return (p + 1.0);
   else
      return p;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG32k5b_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG32k5b_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG32k5 (void *vsta)
{
   MRG32k5_state *state = vsta;
   printf ("(s14, s13, s12, s11, s10) = \n");
   printf ("( %12.0f, %12.0f, %12.0f, %12.0f, %12.0f )\n",
      state->S14, state->S13, state->S12, state->S11, state->S10);
   printf ("\n(s24, s23, s22, s21, s20) = \n");
   printf ("( %12.0f, %12.0f, %12.0f, %12.0f, %12.0f )\n",
      state->S24, state->S23, state->S22, state->S21, state->S20);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateMRG32k5a (double x10, double x11, double x12,
   double x13, double x14, double x20, double x21, double x22,
   double x23, double x24)
{
   unif01_Gen *gen;
   MRG32k5_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG32k5_state));

   strcpy (name, "ulec_CreateMRG32k5a:");
   addstr_Double (name, "   x10 = ", x10);
   addstr_Double (name, ",   x11 = ", x11);
   addstr_Double (name, ",   x12 = ", x12);
   addstr_Double (name, ",   x13 = ", x13);
   addstr_Double (name, ",   x14 = ", x14);
   addstr_Double (name, ",   x20 = ", x20);
   addstr_Double (name, ",   x21 = ", x21);
   addstr_Double (name, ",   x22 = ", x22);
   addstr_Double (name, ",   x23 = ", x23);
   addstr_Double (name, ",   x24 = ", x24);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S10 = x10;
   state->S11 = x11;
   state->S12 = x12;
   state->S13 = x13;
   state->S14 = x14;
   state->S20 = x20;
   state->S21 = x21;
   state->S22 = x22;
   state->S23 = x23;
   state->S24 = x24;

   gen->param = NULL;
   gen->state = state;
   gen->Write = &WrMRG32k5;
   gen->GetBits = &MRG32k5a_Bits;
   gen->GetU01  = &MRG32k5a_U01;
   return gen;
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateMRG32k5b (double x10, double x11, double x12,
   double x13, double x14, double x20, double x21, double x22,
   double x23, double x24)
{
   unif01_Gen *gen;
   MRG32k5_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG32k5_state));

   strcpy (name, "ulec_CreateMRG32k5b:");
   addstr_Double (name, "   x10 = ", x10);
   addstr_Double (name, ",   x11 = ", x11);
   addstr_Double (name, ",   x12 = ", x12);
   addstr_Double (name, ",   x13 = ", x13);
   addstr_Double (name, ",   x14 = ", x14);
   addstr_Double (name, ",   x20 = ", x20);
   addstr_Double (name, ",   x21 = ", x21);
   addstr_Double (name, ",   x22 = ", x22);
   addstr_Double (name, ",   x23 = ", x23);
   addstr_Double (name, ",   x24 = ", x24);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S10 = x10;
   state->S11 = x11;
   state->S12 = x12;
   state->S13 = x13;
   state->S14 = x14;
   state->S20 = x20;
   state->S21 = x21;
   state->S22 = x22;
   state->S23 = x23;
   state->S24 = x24;

   gen->param = NULL;
   gen->state = state;
   gen->Write = &WrMRG32k5;
   gen->GetBits = &MRG32k5b_Bits;
   gen->GetU01  = &MRG32k5b_U01;
   return gen;
}


/***************************************************************************/
#ifdef USE_LONGLONG

#undef norm1
#undef norm2
#undef m1
#undef m2
#undef a12
#undef q12
#undef r12
#undef a13n
#undef q13
#undef r13
#undef a21
#undef q21
#undef r21
#undef a23n
#undef q23
#undef r23

#define norm1 1.0842021724855052e-19
#define norm2 1.0842021724855071e-19
#define m1    9223372036854769163LL
#define m2    9223372036854754679LL
#define a12   1754669720LL
#define q12   5256471877LL
#define r12   251304723LL
#define a13n  3182104042ULL
#define q13   2898513661ULL
#define r13   394451401LL
#define a21   31387477935LL
#define q21   293855150LL
#define r21   143639429LL
#define a23n  6199136374LL
#define q23   1487847900LL
#define r23   985240079LL


static double MRG63k3a_U01 (void *junk, void *vsta)
{
   MRG63k3_state *state = vsta;
   longlong h, p12, p13, p21, p23;

   /* Component 1 */
   h = state->S10 / q13;
   p13 = a13n * (state->S10 - h * q13) - h * r13;
   h = state->S11 / q12;
   p12 = a12 * (state->S11 - h * q12) - h * r12;
   if (p13 < 0)
      p13 += m1;
   if (p12 < 0)
      p12 += m1;
   state->S10 = state->S11;
   state->S11 = state->S12;
   state->S12 = p12 - p13;
   if (state->S12 < 0)
      state->S12 += m1;

   /* Component 2 */
   h = state->S20 / q23;
   p23 = a23n * (state->S20 - h * q23) - h * r23;
   h = state->S22 / q21;
   p21 = a21 * (state->S22 - h * q21) - h * r21;
   if (p23 < 0)
      p23 += m2;
   if (p21 < 0)
      p21 += m2;
   state->S20 = state->S21;
   state->S21 = state->S22;
   state->S22 = p21 - p23;
   if (state->S22 < 0)
      state->S22 += m2;

   /* Combination */
   if (state->S12 - state->S22 > 0)
      return ((state->S12 - state->S22) * norm1);
   else
      return ((state->S12 - state->S22 + m1) * norm1);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG63k3a_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG63k3a_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static double MRG63k3b_U01 (void *junk, void *vsta)
{
   MRG63k3_state *state = vsta;
   longlong h, p12, p13, p21, p23;
   double p;

   /* Component 1 */
   h = state->S10 / q13;
   p13 = a13n * (state->S10 - h * q13) - h * r13;
   h = state->S11 / q12;
   p12 = a12 * (state->S11 - h * q12) - h * r12;
   if (p13 < 0)
      p13 += m1;
   if (p12 < 0)
      p12 += m1;
   state->S10 = state->S11;
   state->S11 = state->S12;
   state->S12 = p12 - p13;
   if (state->S12 < 0)
      state->S12 += m1;

   /* Component 2 */
   h = state->S20 / q23;
   p23 = a23n * (state->S20 - h * q23) - h * r23;
   h = state->S22 / q21;
   p21 = a21 * (state->S22 - h * q21) - h * r21;
   if (p23 < 0)
      p23 += m2;
   if (p21 < 0)
      p21 += m2;
   state->S20 = state->S21;
   state->S21 = state->S22;
   state->S22 = p21 - p23;
   if (state->S22 < 0)
      state->S22 += m2;

   /* Combination */
   p = state->S12 * norm1 - state->S22 * norm2;
   if (p < 0.0)
      return (p + 1.0);
   else
      return p;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG63k3b_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG63k3b_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG63k3 (void *vsta)
{
   MRG63k3_state *state = vsta;
   printf
      ("(s12, s11, s10, s22, s21, s20) = \n( %20" PRIdLEAST64
       ",  %20" PRIdLEAST64 ",  %20" PRIdLEAST64 ",\n  %20" PRIdLEAST64
       ",  %20" PRIdLEAST64 ",  %20" PRIdLEAST64 " )\n",
      state->S12, state->S11, state->S10, state->S22, state->S21,
      state->S20);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_CreateMRG63k3a (longlong s10, longlong s11, longlong s12,
   longlong s20, longlong s21, longlong s22)
{
   unif01_Gen *gen;
   MRG63k3_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG63k3_state));

   strcpy (name, "ulec_CreateMRG63k3a:");
   addstr_LONG (name, "   s10 = ", s10);
   addstr_LONG (name, ",   s11 = ", s11);
   addstr_LONG (name, ",   s12 = ", s12);
   addstr_LONG (name, ",   s20 = ", s20);
   addstr_LONG (name, ",   s21 = ", s21);
   addstr_LONG (name, ",   s22 = ", s22);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S10 = s10;
   state->S11 = s11;
   state->S12 = s12;
   state->S20 = s20;
   state->S21 = s21;
   state->S22 = s22;

   gen->param   = NULL;
   gen->state   = state;
   gen->Write   = &WrMRG63k3;
   gen->GetBits = &MRG63k3a_Bits;
   gen->GetU01  = &MRG63k3a_U01;
   return gen;
}

unif01_Gen *ulec_CreateMRG63k3b (longlong s10, longlong s11, longlong s12,
   longlong s20, longlong s21, longlong s22)
{
   unif01_Gen *gen;
   MRG63k3_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG63k3_state));

   strcpy (name, "ulec_CreateMRG63k3b:");
   addstr_LONG (name, "   s10 = ", s10);
   addstr_LONG (name, ",   s11 = ", s11);
   addstr_LONG (name, ",   s12 = ", s12);
   addstr_LONG (name, ",   s20 = ", s20);
   addstr_LONG (name, ",   s21 = ", s21);
   addstr_LONG (name, ",   s22 = ", s22);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->S10 = s10;
   state->S11 = s11;
   state->S12 = s12;
   state->S20 = s20;
   state->S21 = s21;
   state->S22 = s22;

   gen->param = NULL;
   gen->state = state;
   gen->Write = &WrMRG63k3;
   gen->GetBits = &MRG63k3b_Bits;
   gen->GetU01  = &MRG63k3b_U01;
   return gen;
}


#endif                            /* USE_LONGLONG */

/**************************************************************************/
/*  LFSR Generators.  */

/*  d_j = k_j - s_j.  */

#define q1   13
#define q2    2
#define q3    3
#define d1   19
#define d2   25
#define d3   11
#define s1   12
#define s2    4
#define s3   17
#define msk1 4294967294U            /* (2^32) - 2  */
#define msk2 4294967288U            /* (2^32) - 8  */
#define msk3 4294967280U            /* (2^32) - 16 */


static unsigned long lfsr88_Bits (void *junk, void *vsta)
{
   lfsr88_state *state = vsta;
   unsigned int b;

   b = ((state->z1 << q1) ^ state->z1) >> d1;
   state->z1 = ((state->z1 & msk1) << s1) ^ b;
   b = ((state->z2 << q2) ^ state->z2) >> d2;
   state->z2 = ((state->z2 & msk2) << s2) ^ b;
   b = ((state->z3 << q3) ^ state->z3) >> d3;
   state->z3 = ((state->z3 & msk3) << s3) ^ b;
   return state->z1 ^ state->z2 ^ state->z3;
}

/*-------------------------------------------------------------------------*/

static double lfsr88_U01 (void *junk, void *vsta)
{
   return lfsr88_Bits (junk, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void Wrlfsr88 (void *vsta)
{
   lfsr88_state *state = vsta;
   printf (" z1 = %1u", state->z1);
   printf (",   z2 = %1u", state->z2);
   printf (",   z3 = %1u\n\n", state->z3);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_Createlfsr88 (unsigned int us1, unsigned int us2,
   unsigned int us3)
{
   unif01_Gen *gen;
   lfsr88_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (us1 >= 2, "ulec_Createlfsr88:   s1 < 2");
   util_Assert (us2 >= 8, "ulec_Createlfsr88:   s2 < 8");
   util_Assert (us3 >= 16, "ulec_Createlfsr88:   s3 < 16");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (lfsr88_state));

   strcpy (name, "ulec_Createlfsr88:");
   addstr_Uint (name, "   s1 = ", us1);
   addstr_Uint (name, ",   s2 = ", us2);
   addstr_Uint (name, ",   s3 = ", us3);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->z1 = us1;
   state->z2 = us2;
   state->z3 = us3;

   gen->GetBits = &lfsr88_Bits;
   gen->GetU01 = &lfsr88_U01;
   gen->Write = &Wrlfsr88;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*======================================================================= */

#undef q1
#undef q3
#undef d1
#undef d2
#undef d3
#undef s1
#undef s2
#undef s3

#define q1    6
#define q2    2
#define q3   13
#define q4    3
#define d1   13
#define d2   27
#define d3   21
#define d4   12
#define s1   18
#define s2    2
#define s3    7
#define s4   13
#define msk4 4294967168U            /* 2^(32) - 128 */


static unsigned long lfsr113_Bits (void *junk, void *vsta)
{
   lfsr113_state *state = vsta;
   unsigned int b;

   b = ((state->z1 << q1) ^ state->z1) >> d1;
   state->z1 = ((state->z1 & msk1) << s1) ^ b;
   b = ((state->z2 << q2) ^ state->z2) >> d2;
   state->z2 = ((state->z2 & msk2) << s2) ^ b;
   b = ((state->z3 << q3) ^ state->z3) >> d3;
   state->z3 = ((state->z3 & msk3) << s3) ^ b;
   b = ((state->z4 << q4) ^ state->z4) >> d4;
   state->z4 = ((state->z4 & msk4) << s4) ^ b;
   return state->z1 ^ state->z2 ^ state->z3 ^ state->z4;
}

/*-------------------------------------------------------------------------*/

static double lfsr113_U01 (void *junk, void *vsta)
{
   return lfsr113_Bits (junk, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void Wrlfsr113 (void *vsta)
{
   lfsr113_state *state = vsta;
   printf (" z1 = %1u", state->z1);
   printf (",   z2 = %1u", state->z2);
   printf (",   z3 = %1u", state->z3);
   printf (",   z4 = %1u\n\n", state->z4);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_Createlfsr113 (unsigned int us1, unsigned int us2,
   unsigned int us3, unsigned int us4)
{
   unif01_Gen *gen;
   lfsr113_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (us1 >= 2, "ulec_Createlfsr113:   s1 < 2");
   util_Assert (us2 >= 8, "ulec_Createlfsr113:   s2 < 8");
   util_Assert (us3 >= 16, "ulec_Createlfsr113:   s3 < 16");
   util_Assert (us4 >= 128, "ulec_Createlfsr113:   s4 < 128");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (lfsr113_state));

   strcpy (name, "ulec_Createlfsr113:");
   addstr_Uint (name, "   s1 = ", us1);
   addstr_Uint (name, ",   s2 = ", us2);
   addstr_Uint (name, ",   s3 = ", us3);
   addstr_Uint (name, ",   s4 = ", us4);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->z1 = us1;
   state->z2 = us2;
   state->z3 = us3;
   state->z4 = us4;

   gen->GetBits = &lfsr113_Bits;
   gen->GetU01 = &lfsr113_U01;
   gen->Write = &Wrlfsr113;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/*======================================================================= */
#ifdef USE_LONGLONG

static unsigned long lfsr258_Bits (void *junk, void *vsta)
{
   lfsr258_state *state = vsta;
   ulonglong b;

   b = ((state->y1 << 1) ^ state->y1) >> 53;
   state->y1 = ((state->y1 & 18446744073709551614ULL) << 10) ^ b;
   b = ((state->y2 << 24) ^ state->y2) >> 50;
   state->y2 = ((state->y2 & 18446744073709551104ULL) << 5) ^ b;
   b = ((state->y3 << 3) ^ state->y3) >> 23;
   state->y3 = ((state->y3 & 18446744073709547520ULL) << 29) ^ b;
   b = ((state->y4 << 5) ^ state->y4) >> 24;
   state->y4 = ((state->y4 & 18446744073709420544ULL) << 23) ^ b;
   b = ((state->y5 << 3) ^ state->y5) >> 33;
   state->y5 = ((state->y5 & 18446744073701163008ULL) << 8) ^ b;
   return (state->y1 ^ state->y2 ^ state->y3 ^ state->y4 ^ state->y5) >> 32;
}

/*-------------------------------------------------------------------------*/

static double lfsr258_U01 (void *junk, void *vsta)
{
   return lfsr258_Bits (junk, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void Wrlfsr258 (void *vsta)
{
   lfsr258_state *state = vsta;
   printf (" y1 = %21" PRIuLEAST64 "", state->y1);
   printf (",   y2 = %21" PRIuLEAST64 ",\n", state->y2);
   printf (" y3 = %21" PRIuLEAST64 "", state->y3);
   printf (",   y4 = %21" PRIuLEAST64 ",\n", state->y4);
   printf (" y5 = %21" PRIuLEAST64 "\n\n", state->y5);
}

/*-------------------------------------------------------------------------*/

unif01_Gen *ulec_Createlfsr258 (ulonglong us1, ulonglong us2,
   ulonglong us3, ulonglong us4, ulonglong us5)
{
   unif01_Gen *gen;
   lfsr258_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (us1 >= 2, "ulec_Createlfsr258:   s1 < 2");
   util_Assert (us2 >= 512, "ulec_Createlfsr258:   s2 < 512");
   util_Assert (us3 >= 4096, "ulec_Createlfsr258:   s3 < 4096");
   util_Assert (us4 >= 131072, "ulec_Createlfsr258:   s4 < 131072");
   util_Assert (us5 >= 8388608, "ulec_Createlfsr258:   s5 < 8388608");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (lfsr258_state));

   strcpy (name, "ulec_Createlfsr258:");
   addstr_ULONG (name, "   s1 = ", us1);
   addstr_ULONG (name, ",   s2 = ", us2);
   addstr_ULONG (name, ",   s3 = ", us3);
   addstr_ULONG (name, ",   s4 = ", us4);
   addstr_ULONG (name, ",   s5 = ", us5);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->y1 = us1;
   state->y2 = us2;
   state->y3 = us3;
   state->y4 = us4;
   state->y5 = us5;

   gen->GetBits = &lfsr258_Bits;
   gen->GetU01 = &lfsr258_U01;
   gen->Write = &Wrlfsr258;
   gen->param = NULL;
   gen->state = state;
   return gen;
}

#endif                            /* USE_LONGLONG */


/***************************************************************************/

unif01_Gen *ulec_CreateCombTausLCG11 (unsigned int k, unsigned int q,
   unsigned int s, unsigned int Y1, long m, long a, long c, long Y2)
{
   unif01_Gen *gen1, *gen2;

   gen1 = utaus_CreateTaus (k, q, s, Y1);
   gen2 = ulcg_CreateLCG (m, a, c, Y2);
   return unif01_CreateCombAdd2 (gen1, gen2, "ulec_CreateCombTausLCG11:");
}


void ulec_DeleteCombTausLCG11 (unif01_Gen *gen)
{
   unif01_Comb2_Param *param = gen->param;

   ulcg_DeleteGen (param->gen2);
   utaus_DeleteGen (param->gen1);
   unif01_DeleteCombGen (gen);
}


/***************************************************************************/

#undef q1
#undef q2
#undef s1
#undef s2

unif01_Gen *ulec_CreateCombTausLCG21 (unsigned int k1, unsigned int q1,
   unsigned int s1, unsigned int Y1, unsigned int k2, unsigned int q2,
   unsigned int s2, unsigned int Y2, long m, long a, long c, long Y3)
{
   unif01_Gen *gen1, *gen2;
   double x;

   gen1 = utaus_CreateCombTaus2 (k1, k2, q1, q2, s1, s2, Y1, Y2);
   x = a * (double) m;
   if (((x + c) >= num_TwoExp[53]) || ((-x) >= num_TwoExp[53]))
      gen2 = ulcg_CreateLCG (m, a, c, Y3 % m);
   else
      gen2 = ulcg_CreateLCGFloat (m, a, c, Y3 % m);
   return unif01_CreateCombAdd2 (gen1, gen2, "ulec_CreateCombTausLCG21:");
}


void ulec_DeleteCombTausLCG21 (unif01_Gen *gen)
{
   unif01_Comb2_Param *param = gen->param;

   ulcg_DeleteGen (param->gen2);
   utaus_DeleteGen (param->gen1);
   unif01_DeleteCombGen (gen);
}


/***************************************************************************/

#undef M1
#undef M2
#undef NORM
#undef MASK12
#undef MASK13
#undef MASK21
#define M1      2147483647
#define M2      2147462579
#define NORM    4.656612873077393e-10
#define MASK12  511
#define MASK13  16777215
#define MASK21  65535

static double MRG31k3p_U01 (void *junk, void *vsta)
{
   MRG31k3p_state *state = vsta;
   unsigned long t1, t2;

   /* First component */
   t1 = (((state->x11 & MASK12) << 22) + (state->x11 >> 9))
      + (((state->x12 & MASK13) << 7) + (state->x12 >> 24));
   if (t1 >= M1)
      t1 -= M1;
   t1 += state->x12;
   if (t1 >= M1)
      t1 -= M1;
   state->x12 = state->x11;
   state->x11 = state->x10;
   state->x10 = t1;

   /* Second component */
   t1 = ((state->x20 & MASK21) << 15) + 21069 * (state->x20 >> 16);
   if (t1 >= M2)
      t1 -= M2;
   t2 = ((state->x22 & MASK21) << 15) + 21069 * (state->x22 >> 16);
   if (t2 >= M2)
      t2 -= M2;
   t2 += state->x22;
   if (t2 >= M2)
      t2 -= M2;
   t2 += t1;
   if (t2 >= M2)
      t2 -= M2;
   state->x22 = state->x21;
   state->x21 = state->x20;
   state->x20 = t2;

   /* Combination */
   if (state->x10 <= state->x20)
      return ((state->x10 - state->x20 + M1) * NORM);
   else
      return ((state->x10 - state->x20) * NORM);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG31k3p_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG31k3p_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG31k3p (void *vsta)
{
   MRG31k3p_state *state = vsta;
   printf (
 " x10 = %10lu,  x11 = %10lu,  x12 = %10lu,\n x20 = %10lu,  x21 = %10lu,  x22 = %10lu\n\n",
      state->x10, state->x11, state->x12, state->x20, state->x21, state->x22);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * ulec_CreateMRG31k3p (long x10, long x11, 
  long x12, long x20, long x21, long x22)
{
   unif01_Gen *gen;
   MRG31k3p_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert (!((0 == x10) && (0 == x11) && (0 == x12)),
      "ulec_CreateMRG31k3p:   the first 3 seeds are all 0");
   util_Assert (!((0 == x20) && (0 == x21) && (0 == x22)),
      "ulec_CreateMRG31k3p:   the first 3 seeds are all 0");

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MRG31k3p_state));

   strncpy (name, "ulec_CreateMRG31k3p:   (x10, x11, x12, x20, x21, x22) = ",
     (size_t) LEN);
   addstr_Long (name, "(", x10);
   addstr_Long (name, ", ", x11);
   addstr_Long (name, ", ", x12);
   addstr_Long (name, ", ", x20);
   addstr_Long (name, ", ", x21);
   addstr_Long (name, ", ", x22);
   strncat (name, ")", (size_t) 1);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->x10 = x10 % M1;
   state->x11 = x11 % M1;
   state->x12 = x12 % M1;
   state->x20 = x20 % M2;
   state->x21 = x21 % M2;
   state->x22 = x22 % M2;

   gen->GetBits = &MRG31k3p_Bits;
   gen->GetU01 = &MRG31k3p_U01;
   gen->Write = &WrMRG31k3p;
   gen->param = NULL;
   gen->state = state;
   return gen;
}


/***************************************************************************/

void ulec_DeleteGen (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);

}
