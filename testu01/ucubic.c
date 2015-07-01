/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ucubic.c
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

#include "ucubic.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#define  LEN  150                 /* Max length of strings */



/*================================ Types ================================*/

typedef struct {
   unsigned long m, a, b, c, d;
   double Norm;
} Cubic_param;

typedef struct {
   unsigned long X;
} Cubic_state;


/*---------------*/
#ifdef USE_LONGLONG

typedef struct {
   ulonglong m, a, b, c, d;
   double Norm;
} CubicL_param;

typedef struct {
   ulonglong X;
} CubicL_state;

/*---------------*/
#else

typedef struct {
   long m, a, b, c, d;
   double Norm;
} CubicL_param;

typedef struct {
   long X;
} CubicL_state;

#endif
/*---------------*/


/*-------------------------------------------------------------------------*/

typedef struct {
   double a, b, c, d, m;
   double Norm;
} CubicFloat_param;

typedef struct {
   double X;
} CubicFloat_state;


/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long A, M;
   double Norm;
} Cubic1_param;

typedef struct {
   unsigned long X;
} Cubic1_state;


typedef struct {
   long M, A;
   double Norm;
} Cubic1L_param;

typedef struct {
   long X;
} Cubic1L_state;


/*-------------------------------------------------------------------------*/

typedef struct {
   double M, A;
   double Norm;
} Cubic1Float_param;

typedef struct {
   double X;
} Cubic1Float_state;


/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long M1, M2, A1, A2;
   double Norm1, Norm2;
} CombCubic2_param;

typedef struct {
   unsigned long S1, S2;
} CombCubic2_state;

typedef struct {
   long M1, M2, A1, A2;
   double Norm1, Norm2;
} CombCubic2L_param;

typedef struct {
   long S1, S2;
} CombCubic2L_state;


/*-------------------------------------------------------------------------*/

typedef struct {
   long M, A, C;
   double Norm;
} CubicOut_param;

typedef struct {
   long X;
} CubicOut_state;


/*-------------------------------------------------------------------------*/

typedef struct {
   double M, A, C;
   double Norm;
} CubicOutFloat_param;

typedef struct {
   double X;
} CubicOutFloat_state;




/**************************************************************************/

static double Cubic_U01 (void *vpar, void *vsta)
/*
 * Implementation used when all intermediary results hold in an int of 32
 * bits. This is the case when  m < 2^{16} since a, b, c, d are all < m.
 */
{
   Cubic_param *param = vpar;
   Cubic_state *state = vsta;

   state->X = (param->d + state->X * ((param->c + state->X * ((param->b +
                  state->X * param->a) % param->m)) % param->m)) % param->m;
   return (state->X * param->Norm);
}


/*---------------*/
#ifdef USE_LONGLONG

static double CubicL_U01 (void *vpar, void *vsta)
/*
 * Implementation of Cubic used in the general case.
 */
{
   CubicL_param *param = vpar;
   CubicL_state *state = vsta;

   state->X = (param->d + state->X * ((param->c + state->X * ((param->b +
                  state->X * param->a) % param->m)) % param->m)) % param->m;
   return state->X * param->Norm;
}

/*---------------*/
#else

static double CubicL_U01 (void *vpar, void *vsta)
/*
 * Implementation of Cubic used in the general case. Very slow.
 */
{
   CubicL_param *param = vpar;
   CubicL_state *state = vsta;
   long k;

   k = num_MultModL (param->a, state->X, param->b, param->m);
   k = num_MultModL (k, state->X, param->c, param->m);
   state->X = num_MultModL (k, state->X, param->d, param->m);
   return state->X * param->Norm;
}

#endif
/*---------------*/


static unsigned long Cubic_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Cubic_U01 (vpar, vsta));
}

static unsigned long CubicL_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicL_U01 (vpar, vsta));
}

static void WrCubic (void *vsta)
{
   Cubic_state *state = vsta;
   printf (" X = %lu\n", state->X);
}

static void WrCubicL (void *vsta)
{
   CubicL_state *state = vsta;
   printf (" X = %1lu\n", (unsigned long) state->X);
}


static void SetCubic_0 (unif01_Gen * gen, long m, long a, long b,
   long c, long d, long s)
{
   Cubic_param *param;
   Cubic_state *state;

   param = util_Malloc (sizeof (Cubic_param));
   state = util_Malloc (sizeof (Cubic_state));

   param->Norm = 1.0 / m;
   param->m = m;
   param->a = a;
   param->b = b;
   param->c = c;
   param->d = d;
   state->X = s;

   gen->GetBits = &Cubic_Bits;
   gen->GetU01  = &Cubic_U01;
   gen->Write   = &WrCubic;
   gen->param   = param;
   gen->state   = state;
}

static void SetCubicL_0 (unif01_Gen * gen, long m, long a, long b,
   long c, long d, long s)
{
   CubicL_param *param;
   CubicL_state *state;

   param = util_Malloc (sizeof (CubicL_param));
   state = util_Malloc (sizeof (CubicL_state));

   param->Norm = 1.0 / m;
   param->m = m;
   param->a = a;
   param->b = b;
   param->c = c;
   param->d = d;
   state->X = s;

   gen->GetBits = &CubicL_Bits;
   gen->GetU01  = &CubicL_U01;
   gen->Write   = &WrCubicL;
   gen->param   = param;
   gen->state   = state;
}

unif01_Gen *ucubic_CreateCubic (long m, long a, long b, long c, long d,
   long s)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   util_Assert (m > 0, "ucubic_CreateCubic:   m <= 0");
   util_Assert ((a > 0) && (a < m),
      "ucubic_CreateCubic:   a must be in (0, m)");
   util_Assert ((b >= 0) && (b < m),
      "ucubic_CreateCubic:   b must be in [0, m)");
   util_Assert ((c >= 0) && (c < m),
      "ucubic_CreateCubic:   c must be in [0, m)");
   util_Assert ((d >= 0) && (d < m),
      "ucubic_CreateCubic:   d must be in [0, m)");
   util_Assert ((s >= 0) && (s < m),
      "ucubic_CreateCubic:   s must be in [0, m)");

   gen = util_Malloc (sizeof (unif01_Gen));

   strncpy (name, "ucubic_CreateCubic:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   b = ", b);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   d = ", d);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (m < num_TwoExp[16])
      SetCubic_0 (gen, m, a, b, c, d, s);
   else
      SetCubicL_0 (gen, m, a, b, c, d, s);

   return gen;
}


/**************************************************************************/

static void WrCubicFloat (void *vsta)
{
   CubicFloat_state *state = vsta;
   printf (" X = %1ld\n", (long) state->X);
}

static double CubicFloatA_U01 (void *vpar, void *vsta)
/*
 * Implementation used when
 *    [ a*(m-1)^3 + b*(m-1)^2 + c*(m-1) + d ] / m < 2^31
 * else k may not hold in a long int.
 */
{
   CubicFloat_param *param = vpar;
   CubicFloat_state *state = vsta;
   long k;

   state->X = param->d + state->X * (param->c + state->X * (param->b +
              state->X * param->a));
   k = state->X * param->Norm;
   state->X -= k * param->m;
   return (state->X * param->Norm);
}

static double CubicFloatB_U01 (void *vpar, void *vsta)
/*
 * Implementation used when    (m-1)*m < 2^53
 * else, R*R + R may not hold in a double.
 */
{
   CubicFloat_param *param = vpar;
   CubicFloat_state *state = vsta;
   long k;
   double y;

   y = param->a * state->X + param->b;
   k = y * param->Norm;
   y -= k * param->m;
   y = y * state->X + param->c;
   k = y * param->Norm;
   y -= k * param->m;
   state->X = y * state->X + param->d;
   k = state->X * param->Norm;
   state->X -= k * param->m;
   return state->X * param->Norm;
}

static double CubicFloatC_U01 (void *vpar, void *vsta)
/*
 * Implementation used in the general case.
 */
{
   CubicFloat_param *param = vpar;
   CubicFloat_state *state = vsta;
   double y;

   y = num_MultModD (param->a, state->X, param->b, param->m);
   y = num_MultModD (y, state->X, param->c, param->m);
   state->X = num_MultModD (y, state->X, param->d, param->m);
   return (state->X * param->Norm);
}

static unsigned long CubicFloatA_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicFloatA_U01 (vpar, vsta));
}

static unsigned long CubicFloatB_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicFloatB_U01 (vpar, vsta));
}

static unsigned long CubicFloatC_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicFloatC_U01 (vpar, vsta));
}


unif01_Gen *ucubic_CreateCubicFloat (long m, long a, long b, long c, long d,
   long s)
{
   unif01_Gen *gen;
   CubicFloat_param *param;
   CubicFloat_state *state;
   size_t leng;
   char name[LEN + 1];
   double y, v;

   util_Assert (m > 0, "ucubic_CreateCubicFloat:   m <= 0");
   util_Assert ((a > 0) && (a < m),
      "ucubic_CreateCubicFloat:   a must be in (0, m)");
   util_Assert ((b >= 0) && (b < m),
      "ucubic_CreateCubicFloat:   b must be in [0, m)");
   util_Assert ((c >= 0) && (c < m),
      "ucubic_CreateCubicFloat:   c must be in [0, m)");
   util_Assert ((d >= 0) && (d < m),
      "ucubic_CreateCubicFloat:   d must be in [0, m)");
   util_Assert ((s >= 0) && (s < m),
      "ucubic_CreateCubicFloat:   s must be in [0, m)");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CubicFloat_param));
   state = util_Malloc (sizeof (CubicFloat_state));

   strncpy (name, "ucubic_CreateCubicFloat:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   b = ", b);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   d = ", d);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Norm = 1.0 / m;
   param->m = m;
   param->a = a;
   param->b = b;
   param->c = c;
   param->d = d;
   state->X = s;

   gen->Write   = &WrCubicFloat;
   gen->param   = param;
   gen->state   = state;

   y = m - 1;
   v = d + y * (c + y * (b + y * a));

   if (v / m < num_TwoExp[31]) {
      gen->GetU01  = CubicFloatA_U01;
      gen->GetBits = CubicFloatA_Bits;

   } else if (y * y < num_TwoExp[53]) {
      gen->GetU01  = CubicFloatB_U01;
      gen->GetBits = CubicFloatB_Bits;

   } else {
      gen->GetU01  = CubicFloatC_U01;
      gen->GetBits = CubicFloatC_Bits;
   }

   return gen;
}


/**************************************************************************/

static void WrCubic1 (void *vsta)
{
   Cubic1_state *state = vsta;
   printf (" X = %1lu\n", state->X);
}

static double Cubic1_U01 (void *vpar, void *vsta)
/*
 * Implementation used when all intermediary results hold in an int of 32
 * bits. This is the case when  M < 2^{16}.
 */
{
   Cubic1_param *param = vpar;
   Cubic1_state *state = vsta;

   state->X = (state->X * ((state->X * state->X) % param->M)) % param->M;
   state->X = (param->A * state->X + 1) % param->M;
   return (state->X * param->Norm);
}

static double Cubic1L_U01 (void *vpar, void *vsta)
/*
 * Implementation of Cubic1 used in the general case. Very slow.
 */
{
   Cubic1L_param *param = vpar;
   Cubic1L_state *state = vsta;
   long S;

   S = num_MultModL (state->X, state->X, 0L, param->M);
   S = num_MultModL (state->X, S, 0L, param->M);
   state->X = num_MultModL (param->A, S, 1L, param->M);
   return (state->X * param->Norm);
}

static unsigned long Cubic1_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Cubic1_U01 (vpar, vsta));
}

static unsigned long Cubic1L_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Cubic1L_U01 (vpar, vsta));
}


static void SetCubic1_0 (unif01_Gen * gen, long m, long a, long s)
{
   Cubic1_param *param;
   Cubic1_state *state;

   param = util_Malloc (sizeof (Cubic1_param));
   state = util_Malloc (sizeof (Cubic1_state));

   param->Norm = 1.0 / m;
   param->M = m;
   param->A = a;
   state->X = s;

   gen->GetBits = &Cubic1_Bits;
   gen->GetU01  = &Cubic1_U01;
   gen->Write   = &WrCubic1;
   gen->param   = param;
   gen->state   = state;
}

static void SetCubic1L_0 (unif01_Gen * gen, long m, long a, long s)
{
   Cubic1L_param *param;
   Cubic1L_state *state;

   param = util_Malloc (sizeof (Cubic1L_param));
   state = util_Malloc (sizeof (Cubic1L_state));

   param->Norm = 1.0 / m;
   param->M = m;
   param->A = a;
   state->X = s;

   gen->GetBits = &Cubic1L_Bits;
   gen->GetU01  = &Cubic1L_U01;
   gen->Write   = &WrCubic1;
   gen->param   = param;
   gen->state   = state;
}


unif01_Gen *ucubic_CreateCubic1 (long m, long a, long s)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   util_Assert ((m > 0), "ucubic_CreateCubic1:   m <= 0");
   util_Assert ((a > 0) && (a < m),
      "ucubic_CreateCubic1:   a must be in (0, m)");
   util_Assert ((s >= 0) && (s < m),
      "ucubic_CreateCubic1:   s must be in [0, m)");

   gen = util_Malloc (sizeof (unif01_Gen));

   strncpy (name, "ucubic_CreateCubic1:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (m < num_TwoExp[16])
      SetCubic1_0 (gen, m, a, s);
   else
      SetCubic1L_0 (gen, m, a, s);

   return gen;
}


/**************************************************************************/

static void WrCubic1Float (void *vsta)
{
   Cubic1Float_state *state = vsta;
   printf (" X = %1ld\n", (long) state->X);
}

static double Cubic1FloatA_U01 (void *vpar, void *vsta)
/*
 * Used when condition 1 + a*(m-1)^3/m < 2^31; 
 * else k may not hold in a long.
 */
{
   Cubic1Float_param *param = vpar;
   Cubic1Float_state *state = vsta;
   long k;

   state->X = param->A * state->X * state->X * state->X + 1.0;
   k = state->X * param->Norm;
   state->X -= k * param->M;
   return (state->X * param->Norm);
}

static double Cubic1FloatB_U01 (void *vpar, void *vsta)
/*
 * Used when condition (m-1)^2 < 2^53;
 * else R*R may not hold in a double.
 */
{
   Cubic1Float_param *param = vpar;
   Cubic1Float_state *state = vsta;
   long k;
   double y;

   y = param->A * state->X;
   k = y * param->Norm;
   y -= k * param->M;
   y *= state->X;
   k = y * param->Norm;
   y -= k * param->M;
   state->X = y * state->X + 1.0;
   k = state->X * param->Norm;
   state->X -= k * param->M;
   return (state->X * param->Norm);
}

static double Cubic1FloatC_U01 (void *vpar, void *vsta)
/*
 * Inmplementation used in the general case
 */
{
   Cubic1Float_param *param = vpar;
   Cubic1Float_state *state = vsta;
   double y;

   y = num_MultModD (state->X, state->X, 0.0, param->M);
   y = num_MultModD (state->X, y, 0.0, param->M);
   state->X = num_MultModD (param->A, y, 1.0, param->M);
   return (state->X * param->Norm);
}

static unsigned long Cubic1FloatA_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Cubic1FloatA_U01 (vpar, vsta));
}

static unsigned long Cubic1FloatB_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Cubic1FloatB_U01 (vpar, vsta));
}

static unsigned long Cubic1FloatC_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Cubic1FloatC_U01 (vpar, vsta));
}

unif01_Gen *ucubic_CreateCubic1Float (long m, long a, long s)
{
   unif01_Gen *gen;
   Cubic1Float_param *param;
   Cubic1Float_state *state;
   size_t leng;
   char name[LEN + 1];
   double y;

   util_Assert ((m > 0), "ucubic_CreateCubic1Float:   m <= 0");
   util_Assert ((a > 0) && (a < m),
      "ucubic_CreateCubic1Float:   a must be in (0, m)");
   util_Assert ((s >= 0) && (s < m),
      "ucubic_CreateCubic1Float:   s must be in [0, m)");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Cubic1Float_param));
   state = util_Malloc (sizeof (Cubic1Float_state));

   strncpy (name, "ucubic_CreateCubic1Float:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->Norm = 1.0 / m;
   param->M = m;
   param->A = a;
   state->X = s;

   y = m - 1;
   if (1.0 + a * y * y * y / m < num_TwoExp[31]) {
      gen->GetU01  = Cubic1FloatA_U01;
      gen->GetBits = Cubic1FloatA_Bits;

   } else if (y * y < num_TwoExp[53]) {
      gen->GetU01  = Cubic1FloatB_U01;
      gen->GetBits = Cubic1FloatB_Bits;

   } else {
      gen->GetU01  = Cubic1FloatC_U01;
      gen->GetBits = Cubic1FloatC_Bits;
   }

   gen->Write = &WrCubic1Float;
   gen->param = param;
   gen->state = state;

   return gen;
}


/**************************************************************************/

static void WrCombCubic2 (void *vsta)
{
   CombCubic2_state *state = vsta;
   printf (" X1 = %1lu,   X2 = %1lu\n", state->S1, state->S2);
}

static double CombCubic2_U01 (void *vpar, void *vsta)
/*
 * Implementation used when all intermediary results hold in an int of 32
 * bits. This is the case when  M < 2^{16}. 
 */
{
   CombCubic2_param *param = vpar;
   CombCubic2_state *state = vsta;
   unsigned long Y;
   double U;

   Y = (state->S1 * ((state->S1 * state->S1) % param->M1)) % param->M1;
   state->S1 = (param->A1 * Y + 1) % param->M1;
   Y = (state->S2 * ((state->S2 * state->S2) % param->M2)) % param->M2;
   state->S2 = (param->A2 * Y + 1) % param->M2;
   U = state->S1 * param->Norm1 + state->S2 * param->Norm2;
   if (U >= 1.0)
      return (U - 1.0);
   else
      return U;
}

static double CombCubic2L_U01 (void *vpar, void *vsta)
/*
 * Implementation of CombCubic2 used in the general case. Very slow.
 */
{
   CombCubic2L_param *param = vpar;
   CombCubic2L_state *state = vsta;
   long Y;
   double U;

   Y = num_MultModL (state->S1, state->S1, 0L, param->M1);
   Y = num_MultModL (state->S1, Y, 0L, param->M1);
   state->S1 = num_MultModL (param->A1, Y, 1L, param->M1);
   Y = num_MultModL (state->S2, state->S2, 0L, param->M2);
   Y = num_MultModL (state->S2, Y, 0L, param->M2);
   state->S2 = num_MultModL (param->A2, Y, 1L, param->M2);
   U = state->S1 * param->Norm1 + state->S2 * param->Norm2;
   if (U >= 1.0)
      return (U - 1.0);
   else
      return U;
}

static unsigned long CombCubic2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombCubic2_U01 (vpar, vsta));
}

static unsigned long CombCubic2L_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombCubic2L_U01 (vpar, vsta));
}


static void SetCombCubic2_0 (unif01_Gen * gen, long m1, long m2, long a1,
   long a2, long s1, long s2)
{
   CombCubic2_param *param;
   CombCubic2_state *state;

   param = util_Malloc (sizeof (CombCubic2_param));
   state = util_Malloc (sizeof (CombCubic2_state));

   param->Norm1 = 1.0 / m1;
   param->Norm2 = 1.0 / m2;
   param->M1 = m1;
   param->A1 = a1;
   param->M2 = m2;
   param->A2 = a2;
   state->S1 = s1;
   state->S2 = s2;

   gen->GetBits = &CombCubic2_Bits;
   gen->GetU01  = &CombCubic2_U01;
   gen->Write   = &WrCombCubic2;
   gen->param   = param;
   gen->state   = state;
}

static void SetCombCubic2L_0 (unif01_Gen * gen, long m1, long m2, long a1,
   long a2, long s1, long s2)
{
   CombCubic2L_param *param;
   CombCubic2L_state *state;

   param = util_Malloc (sizeof (CombCubic2L_param));
   state = util_Malloc (sizeof (CombCubic2L_state));

   param->Norm1 = 1.0 / m1;
   param->Norm2 = 1.0 / m2;
   param->M1 = m1;
   param->A1 = a1;
   param->M2 = m2;
   param->A2 = a2;
   state->S1 = s1;
   state->S2 = s2;

   gen->GetBits = &CombCubic2L_Bits;
   gen->GetU01  = &CombCubic2L_U01;
   gen->Write   = &WrCombCubic2;
   gen->param   = param;
   gen->state   = state;
}


unif01_Gen *ucubic_CreateCombCubic2 (long m1, long m2, long a1, long a2,
   long s1, long s2)
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   if ((a1 <= 0) || (a1 >= m1) || (s1 < 0) || (s1 >= m1) || (m1 <= 0) ||
      (a2 <= 0) || (a2 >= m2) || (s2 < 0) || (s2 >= m2) || (m2 <= 0)) {
      printf ("m1 = %1ld,  m2 = %1ld,  a1 = %1ld,  a2 = %1ld,"
         " s1 = %1ld,  s2 = %1ld\n", m1, m2, a1, a2, s1, s2);
      util_Error ("ucubic_CreateCombCubic2:   Invalid parameter");
   }

   gen = util_Malloc (sizeof (unif01_Gen));

   strncpy (name, "ucubic_CreateCombCubic2:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   a1 = ", a1);
   addstr_Long (name, ",   a2 = ", a2);
   addstr_Long (name, ",   s1 = ", s1);
   addstr_Long (name, ",   s2 = ", s2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if ((m1 < num_TwoExp[16]) && (m2 < num_TwoExp[16]))
      SetCombCubic2_0 (gen, m1, m2, a1, a2, s1, s2);
   else
      SetCombCubic2L_0 (gen, m1, m2, a1, a2, s1, s2);
   return gen;
}


/**************************************************************************/

static void WrCubicOut (void *vsta)
{
   CubicOut_state *state = vsta;
   printf (" X = %1ld\n", state->X);
}

static double CubicOut_U01 (void *vpar, void *vsta)
{
   CubicOut_param *param = vpar;
   CubicOut_state *state = vsta;
   long Y;

   state->X = num_MultModL (param->A, state->X, param->C, param->M);
   Y = num_MultModL (state->X, state->X, 0L, param->M);
   Y = num_MultModL (state->X, Y, 0L, param->M);
   return (Y * param->Norm);
}

static unsigned long CubicOut_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicOut_U01 (vpar, vsta));
}


unif01_Gen *ucubic_CreateCubicOut (long m, long a, long c, long s)
{
   unif01_Gen *gen;
   CubicOut_param *param;
   CubicOut_state *state;
   size_t leng;
   char name[LEN + 1];

   util_Assert ((m > 0), "ucubic_CreateCubicOut:   m <= 0");
   util_Assert ((a > 0) && (a < m),
      "ucubic_CreateCubicOut:   a must be in (0, m)");
   util_Assert ((c >= 0) && (c < m),
      "ucubic_CreateCubicOut:   c must be in [0, m)");
   util_Assert ((s >= 0) && (s < m),
      "ucubic_CreateCubicOut:   s must be in [0, m)");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CubicOut_param));
   state = util_Malloc (sizeof (CubicOut_state));

   strncpy (name, "ucubic_CreateCubicOut:", (size_t) LEN);
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
   state->X = s;

   gen->GetU01  = CubicOut_U01;
   gen->GetBits = CubicOut_Bits;
   gen->Write   = &WrCubicOut;
   gen->param   = param;
   gen->state   = state;
   return gen;
}


/**************************************************************************/

static void WrCubicOutFloat (void *vsta)
{
   CubicOutFloat_state *state = vsta;
   printf (" X = %1ld\n", (long) state->X);
}

static double CubicOutFloatA_U01 (void *vpar, void *vsta)
/*
 * Used when condition (m-1)^3/m < 2^31;
 * else k may not hold in a long.
 */
{
   CubicOutFloat_param *param = vpar;
   CubicOutFloat_state *state = vsta;
   long k;
   double y;

   state->X = param->A * state->X + param->C;
   k = state->X * param->Norm;
   state->X -= k * param->M;
   y = state->X * state->X * state->X;
   k = y * param->Norm;
   y -= k * param->M;
   return (y * param->Norm);
}

static double CubicOutFloatB_U01 (void *vpar, void *vsta)
/*
 * Floating-point implementation for the case when m <= 94906266 ~ 2^26.5;
 * follows the condition (m-1)^2 < 2^53; else, R*R may not hold in a double.
 */
{
   CubicOutFloat_param *param = vpar;
   CubicOutFloat_state *state = vsta;
   long k;
   double y;

   state->X = param->A * state->X + param->C;
   k = state->X * param->Norm;
   state->X -= k * param->M;
   y = state->X * state->X;
   k = y * param->Norm;
   y -= k * param->M;
   y = state->X * y;
   k = y * param->Norm;
   y -= k * param->M;
   return (y * param->Norm);
}

static double CubicOutFloatC_U01 (void *vpar, void *vsta)
{
   CubicOutFloat_param *param = vpar;
   CubicOutFloat_state *state = vsta;
   double y;

   state->X = num_MultModD (param->A, state->X, param->C, param->M);
   y = num_MultModD (state->X, state->X, 0.0, param->M);
   y = num_MultModD (state->X, y, 0.0, param->M);
   return (y * param->Norm);
}

static unsigned long CubicOutFloatA_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicOutFloatA_U01 (vpar, vsta));
}

static unsigned long CubicOutFloatB_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicOutFloatB_U01 (vpar, vsta));
}

static unsigned long CubicOutFloatC_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CubicOutFloatC_U01 (vpar, vsta));
}

unif01_Gen *ucubic_CreateCubicOutFloat (long m, long a, long c, long s)
{
   unif01_Gen *gen;
   CubicOutFloat_param *param;
   CubicOutFloat_state *state;
   size_t leng;
   char name[LEN + 1];
   double y;

   util_Assert ((m > 0), "ucubic_CreateCubicOutFloat:   m <= 0");
   util_Assert ((a > 0) && (a < m),
      "ucubic_CreateCubicOutFloat:   a must be in (0, m)");
   util_Assert ((c >= 0) && (c < m),
      "ucubic_CreateCubicOutFloat:   c must be in [0, m)");
   util_Assert ((s >= 0) && (s < m),
      "ucubic_CreateCubicOutFloat:   s must be in [0, m)");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CubicOutFloat_param));
   state = util_Malloc (sizeof (CubicOutFloat_state));

   strncpy (name, "ucubic_CreateCubicOutFloat:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   a = ", a);
   addstr_Long (name, ",   c = ", c);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->M = m;
   param->A = a;
   state->X = s;
   param->C = c;
   param->Norm = 1.0 / m;

   y = m - 1;
   if (y * y * y / m < num_TwoExp[31]) {
      gen->GetU01  = CubicOutFloatA_U01;
      gen->GetBits = CubicOutFloatA_Bits;

   } else if (y * y < num_TwoExp[53]) {
      gen->GetU01  = CubicOutFloatB_U01;
      gen->GetBits = CubicOutFloatB_Bits;

   } else {
      gen->GetU01  = CubicOutFloatC_U01;
      gen->GetBits = CubicOutFloatC_Bits;
   } 

   gen->Write = &WrCubicOutFloat;
   gen->param = param;
   gen->state = state;

   return gen;
}


/**************************************************************************/

void ucubic_DeleteGen (unif01_Gen * gen)
{
   unif01_DeleteGen (gen);
}
