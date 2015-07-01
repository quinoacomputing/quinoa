/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ucarry.c
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

#include "ucarry.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/*=============================== Constants ===============================*/

#define  LEN 300

#define  IDEF    314159265         /* Constants for Ranlux */
#define  ICONS   2147483563
#define  LUXS    10
#define  LUXR    24
#define  TwoM24  0.596046447753906250E-7     /* 1 / 2^24 */
#define  TwoM12  0.2441406250E-3             /* 1 / 2^12 */


#define MASK16  0xffff            /* Mask of 16 bits */
#define MASK32  0xffffffffU       /* Mask of 32 bits */




/*=============================== Types =================================*/

typedef struct {
   unsigned long M;                /* Modulus */
   double Norm;                    /* 1 / M */
   lebool rFlag;                  /* = TRUE if r >= s, FALSE if r < s */
   int skip;                       /* For luxury: skip = Lux - S numbers */
} AWC_param;

typedef struct {
   unsigned long *X;               /* State */
   unsigned int C;                 /* Carry */
   unsigned int R, S;              /* Lags */
   unsigned int n;                 /* how many since last skip */
   unsigned int Order;             /* Order */
} AWC_state;

/*-----------------------------------------------------------------------*/

typedef AWC_param SWB_param;
typedef AWC_state SWB_state;

/*-----------------------------------------------------------------------*/

typedef struct {
   int Next[LUXR + 1];
   int skip;
} Ranlux_param;

typedef struct {
   double X[LUXR + 1];             /* State */
   double C;                       /* Carry */
   unsigned int R, S;              /* Lags */
   unsigned int RR;                /* Skip index */
} Ranlux_state;

/*-----------------------------------------------------------------------*/
#ifdef USE_LONGLONG

typedef struct {
   unsigned long *A;
   unsigned int W;
   unsigned int shift;
   ulonglong mask;
} MWC_param;

typedef struct {
   unsigned long *X;               /* State */
   ulonglong C;                    /* Carry */
   unsigned int R;
   unsigned int Order;             /* Order */
} MWC_state;

#endif
/*-----------------------------------------------------------------------*/

typedef struct {
   unsigned long *A;
   unsigned int W;
   unsigned int shift;
   unsigned long mask;
   unsigned long aa;               /* a0^(-1) mod 2^w */
   double Norm;                    /* 1 / 2^w */
} MWCFloat_param;

typedef struct {
   unsigned long *X;               /* State */
   unsigned long C;                /* Carry */
   unsigned int R;
   unsigned int Order;             /* Order */
} MWCFloat_state;

/*-----------------------------------------------------------------------*/

typedef struct {
   unsigned int x1, x2, x3, x4,
                 x5, x6, x7, x8;   /* State */
   unsigned int C;                /* Carry */
} MWCfixCouture_state;

/*-----------------------------------------------------------------------*/

typedef struct {
   unsigned int *A;
   unsigned int Asize;
   unsigned int W;
   unsigned int shift;
   unsigned int mask;
} SWC_param;

typedef struct {
   unsigned int *X;               /* State */
   unsigned int C;                /* Carry */
   unsigned int R, S;
   unsigned int Order;             /* Order */
} SWC_state;

/*-----------------------------------------------------------------------*/

typedef struct {
   unsigned int a, b;
} MWC1616_param;

typedef struct {
   unsigned int x, y;
} MWC1616_state;






/*============================= Functions ===============================*/

/* CreateGen common to generators AWC et SWB */

static unif01_Gen *Create_AWC_SWB (unsigned int r, unsigned int s,
   unsigned long c, unsigned long m, unsigned long S[], const char nom[])
{
   unif01_Gen *gen;
   AWC_param *param;
   AWC_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int j, rs;

   if ((c != 0 && c != 1) || (s < 1) || (r < 1) || (r == s)) {
      strcpy (name, nom);
      strcat (name, ":   invalid parameter");
      util_Error (name);
   }
   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (AWC_param));
   state = util_Malloc (sizeof (AWC_state));
   if (r > s) {
      param->rFlag = TRUE;
      rs = r;
   } else {
      rs = s;
      param->rFlag = FALSE;
   }
/*   param->skip = Lux - rs; */
   state->X = util_Calloc ((size_t) rs + 1, sizeof (unsigned long));

   strncpy (name, nom, (size_t) LEN);
   strcat (name, ":   ");
   addstr_Uint (name, "r = ", r);
   addstr_Uint (name, ",   s = ", s);
   addstr_Ulong (name, ",   c = ", c);
   addstr_Ulong (name, ",   m = ", m);
/*   if (Lux > 0)
      addstr_Int (name, ",   Lux = ", Lux); */
   addstr_ArrayUlong (name, ",   S = ", (int) rs, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (m == 0) {
      for (j = 0; j < rs; j++)
         state->X[j] = S[j] & MASK32;
      util_Warning (1, "AWC or SWB:   m = 0;  I will assume m = 2^32");
   } else
      for (j = 0; j < rs; j++)
         state->X[j] = S[j] % m;
   
   if (param->rFlag)
      state->S = r - s;
   else
      state->S = s - r;
   state->R = 0;
   state->C = c;
   state->Order = rs;
   state->n = 0;
   param->M = m;
   if (m == 0)
      param->Norm = unif01_INV32;
   else
      param->Norm = 1.0 / m;
   gen->param = param;
   gen->state = state;
   return gen;
}

/***************************************************************************/

static double AWC_U01 (void *vpar, void *vsta)
{
   AWC_param *param = vpar;
   AWC_state *state = vsta;
   unsigned long x, temp;

   x = state->X[state->R] + state->C;
   temp = param->M - state->X[state->S];
   if (x < temp) {
      x += state->X[state->S];
      state->C = 0;
   } else {
      x -= temp;
      state->C = 1;
   }
   state->X[state->R] = x;
   state->R = (state->R + 1) % state->Order;
   state->S = (state->S + 1) % state->Order;
   return (x * param->Norm);
}

static unsigned long AWC_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * AWC_U01 (vpar, vsta));
}

static void WrAWC (void *vsta)
{
   AWC_state *state = vsta;
   unsigned int j;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < state->Order; j++) {
         printf (" %12lu", state->X[j]);
         if (j < state->Order - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("   }\n\n");
      printf (" c = %1u\n\n", state->C);
   } else
      unif01_WrLongStateDef ();
}

unif01_Gen *ucarry_CreateAWC (unsigned int r, unsigned int s,
   unsigned long c, unsigned long m, unsigned long S[])
{
   unif01_Gen *gen;

   gen = Create_AWC_SWB (r, s, c, m, S, "ucarry_CreateAWC");
   gen->GetBits = &AWC_Bits;
   gen->GetU01  = &AWC_U01;
   gen->Write   = &WrAWC;
   return gen;
}

void ucarry_DeleteAWC (unif01_Gen * gen)
{
   AWC_state *state;

   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->X);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/***************************************************************************/

static double SWB_U01 (void *vpar, void *vsta)
{
   SWB_param *param = vpar;
   SWB_state *state = vsta;
   unsigned long x, y;

   if (param->rFlag) {
      x = state->X[state->R];
      y = state->X[state->S] + state->C;
   } else {
      x = state->X[state->S];
      y = state->X[state->R] + state->C;
   }
   if (x < y) {
      x = (param->M - y) + x;
      state->C = 1;
   } else {
      x -= y;
      state->C = 0;
   }
   state->X[state->R] = x;

   state->S++;
   if (state->S == state->Order)
      state->S = 0;
   state->R++;
   if (state->R == state->Order)
      state->R = 0;
#if 0
   state->n++;
   if (state->n == state->Order) {
      unsigned long z;
      int i;
      for (i = 0; i < param->skip; i++) {
         if (param->rFlag) {
            z = state->X[state->R];
            y = state->X[state->S] + state->C;
         } else {
            z = state->X[state->S];
            y = state->X[state->R] + state->C;
         }
         if (z < y) {
            z = (param->M - y) + z;
            state->C = 1;
         } else {
            z -= y;
            state->C = 0;
         }
         state->X[state->R] = z;

         state->S++;
         if (state->S == state->Order)
            state->S = 0;
         state->R++;
         if (state->R == state->Order)
            state->R = 0;
      }   
      state->n = 0;
   }
#endif

   return (x * param->Norm);
}

static unsigned long SWB_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * SWB_U01 (vpar, vsta));
}

static void WrSWB (void *vsta)
{
   WrAWC (vsta);
}

unif01_Gen *ucarry_CreateSWB (unsigned int r, unsigned int s,
   unsigned long c, unsigned long m, unsigned long S[])
{
   unif01_Gen *gen;

   gen = Create_AWC_SWB (r, s, c, m, S, "ucarry_CreateSWB");
   gen->GetBits = &SWB_Bits;
   gen->GetU01  = &SWB_U01;
   gen->Write   = &WrSWB;
   return gen;
}

void ucarry_DeleteSWB (unif01_Gen * gen)
{
   ucarry_DeleteAWC (gen);
}


/**************************************************************************/

static double Ranlux_U01 (void *vpar, void *vsta)
{
   Ranlux_param *param = vpar;
   Ranlux_state *state = vsta;
   double Uni, wi;
   int i;

   Uni = state->X[state->S] - state->X[state->R] - state->C;
   state->C = 0.0;
   if (Uni < 0.0) {
      Uni += 1.0;
      state->C = TwoM24;
   }
   state->X[state->R] = Uni;
   state->R = param->Next[state->R];
   state->S = param->Next[state->S];

   wi = Uni;
   if (Uni < TwoM12) {
      wi += state->X[state->S] * TwoM24;
      if (wi == 0.0) {
         wi = TwoM24 * TwoM24;
      }
   }
   /* Skipping to luxury: generate LUXR values and skip (Lux - LUXR) values */
   state->RR++;
   if (state->RR == LUXR) {
      state->RR = 0;
      for (i = 1; i <= param->skip; i++) {
         Uni = state->X[state->S] - state->X[state->R] - state->C;
         state->C = 0.0;
         if (Uni < 0.0) {
            Uni += 1.0;
            state->C = TwoM24;
         }
         state->X[state->R] = Uni;
         state->R = param->Next[state->R];
         state->S = param->Next[state->S];
      }
   }
   return wi;
}

static unsigned long Ranlux_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * Ranlux_U01 (vpar, vsta));
}

static void WrRanlux (void *vsta)
{
   Ranlux_state *state = vsta;
   int j;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 1; j <= LUXR; j++) {
         printf (" %10.7f", state->X[j]);
         if (j < LUXR)
            printf (", ");
         if (((j - 1) % 5) == 4)
            printf ("\n ");
      };
      printf ("  }\n\n");
   } else
      unif01_WrLongStateDef ();
}

unif01_Gen *ucarry_CreateRanlux (unsigned int Lux, long s)
{
   unif01_Gen *gen;
   Ranlux_param *param;
   Ranlux_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;
   long z, x0, temp;

   if (Lux < 24)
      util_Error ("ucarry_CreateRanlux:   Lux < 24");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (Ranlux_param));
   state = util_Malloc (sizeof (Ranlux_state));

   strncpy (name, "ucarry_CreateRanlux:", (size_t) LEN);
   addstr_Ulong (name, "   Lux = ", Lux);
   addstr_Long (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->skip = Lux - LUXR;
   state->RR = 0;
   if (s == 0)
      s = IDEF;                     /* default initial value */
   z = s;
   for (i = 1; i <= LUXR; i++) {
      x0 = z / 53668;
      z = 40014 * (z - x0 * 53668) - x0 * 12211;
      if (z < 0)
         z += ICONS;
      temp = z % 16777216;
      state->X[i] = temp * TwoM24;
      param->Next[i] = i - 1;
   }
   param->Next[1] = LUXR;
   state->R = LUXR;
   state->S = LUXS;
   state->C = 0.0;
   if (state->X[LUXR] == 0.0)
      state->C = TwoM24;

   gen->GetBits = &Ranlux_Bits;
   gen->GetU01  = &Ranlux_U01;
   gen->Write   = &WrRanlux;
   gen->param   = param;
   gen->state   = state;
   return gen;
}

void ucarry_DeleteRanlux (unif01_Gen * gen)
{
   unif01_DeleteGen (gen);
}


/************************************************************************/
#ifdef USE_LONGLONG

static unsigned long MWC_Bits (void *vpar, void *vsta)
{
   MWC_param *param = vpar;
   MWC_state *state = vsta;
   unsigned int i, j;
   unsigned long SomC;

   for (j = 0; j < state->Order; j++) {
      if (param->A[j] != 0) {
         i = j + state->R;
         if (i >= state->Order)
            i -= state->Order;
         state->C += param->A[j] * (ulonglong) state->X[i];
      }
   }
   SomC = state->C & param->mask;
   state->X[state->R] = SomC;
   state->C >>= param->W;
   state->R++;
   if (state->R >= state->Order)
      state->R = 0;
   return (SomC << param->shift) & MASK32;
}

static double MWC_U01 (void *vpar, void *vsta)
{
   return (MWC_Bits (vpar, vsta) * unif01_INV32);
}


static void WrMWC (void *vsta)
{
   MWC_state *state = vsta;
   unsigned int j;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < state->Order - 1; j++) {
	 printf ("%10lu,  ", state->X[j]);
	 if (((j + 1) % 5) == 0)
	    printf ("\n ");
      }
      printf ("%10lu   }\n\n", state->X[state->Order - 1]);
      printf (" c = %1" PRIuLEAST64 "\n\n", state->C);
   } else
      unif01_WrLongStateDef ();
}

unif01_Gen *ucarry_CreateMWC (unsigned int r, unsigned long c,
   unsigned int w, unsigned long A[], unsigned long S[])
{
   unif01_Gen *gen;
   MWC_param *param;
   MWC_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int j;
   ulonglong Suma;

   if (w > 32)
      util_Error ("ucarry_CreateMWC:   w > 32");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (MWC_param));
   state = util_Malloc (sizeof (MWC_state));
   state->X = util_Calloc ((size_t) r, sizeof (unsigned long));
   param->A = util_Calloc ((size_t) r, sizeof (unsigned long));

   strncpy (name, "ucarry_CreateMWC:", (size_t) LEN);
   addstr_Uint (name, "   r = ", r);
   addstr_Ulong (name, ",   c = ", c);
   addstr_Uint (name, ",   w = ", w);
   addstr_ArrayUlong (name, ",   A = ", (int) r, A);
   addstr_ArrayUlong (name, ",   S = ", (int) r, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   Suma = 0;
   for (j = 0; j < r; j++) {
      if (A[j] >= num_TwoExp[w])
         util_Error ("ucarry_CreateMWC:   A[i] must be < 2^w");
      if (S[j] >= num_TwoExp[w])
         util_Error ("ucarry_CreateMWC:   S[i] must be < 2^w");
      Suma += A[j];
   }
   /* so that carry c doesn't lose digits of precision */
   Suma = Suma * ((ulonglong) num_TwoExp[w] - 1) + c;
   if (Suma >= num_TwoExp[64])
      util_Error ("ucarry_CreateMWC:   Sum over A[i] is too big");

   state->C = (ulonglong) c;
   state->R = 0;
   state->Order = r;
   param->W = w;
   param->shift = 32 - w;
   if (w < 32) {
      param->mask = (ulonglong) num_TwoExp[w] - 1;
   } else {
      param->mask = (ulonglong) 0xffffffffUL;
   }

   for (j = 0; j < r; j++) {
      param->A[j] = A[j];
      state->X[j] = S[j];
   }

   gen->param = param;
   gen->state = state;
   gen->GetBits = &MWC_Bits;
   gen->GetU01  = &MWC_U01;
   gen->Write   = &WrMWC;
   return gen;
}

void ucarry_DeleteMWC (unif01_Gen * gen)
{
   MWC_param *param;
   MWC_state *state;

   if (NULL == gen)
      return;
   state = gen->state;
   param = gen->param;
   util_Free (state->X);
   util_Free (param->A);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


#endif  /* USE_LONGLONG */
/**************************************************************************/

static unsigned long MWCFloat_Bits (void *vpar, void *vsta)
{
   MWCFloat_param *param = vpar;
   MWCFloat_state *state = vsta;

   unsigned int i, j;
   double dv;
   unsigned long SomC;

   dv = state->C;
   for (j = 0; j < state->Order; j++) {
      if (param->A[j] != 0) {
         i = j + state->R;
         if (i >= state->Order)
            i -= state->Order;
         state->C += param->A[j] * state->X[i];
         dv += (double) param->A[j] * (double) state->X[i];
      }
   }
   SomC = state->C & param->mask;
   state->X[state->R] = SomC;
   state->C = dv * param->Norm;
   state->R++;
   if (state->R >= state->Order)
      state->R = 0;
   return (SomC << param->shift) & MASK32;
}

static double MWCFloat_U01 (void *vpar, void *vsta)
{
   return (MWCFloat_Bits (vpar, vsta) * unif01_INV32);
}

static void WrMWCFloat (void *vsta)
{
   MWCFloat_state *state = vsta;
   unsigned int j;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < state->Order - 1; j++) {
	 printf ("%10lu,  ", state->X[j]);
	 if (((j + 1) % 5) == 0)
	    printf ("\n ");
      }
      printf ("%10lu   }\n\n", state->X[state->Order - 1]);
      printf (" c = %1lu\n\n", state->C);
   } else
      unif01_WrLongStateDef ();
}

unif01_Gen *ucarry_CreateMWCFloat (unsigned int r, unsigned long c,
   unsigned int w, unsigned long A[], unsigned long S[])
{
   unif01_Gen *gen;
   MWCFloat_param *param;
   MWCFloat_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int j;
   double Sum;

   util_Assert (w <= 32, "ucarry_CreateMWCFloat:   w > 32");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (MWCFloat_param));
   state = util_Malloc (sizeof (MWCFloat_state));
   state->X = util_Calloc ((size_t) r, sizeof (unsigned long));
   param->A = util_Calloc ((size_t) r, sizeof (unsigned long));

   strncpy (name, "ucarry_CreateMWCFloat:", (size_t) LEN);
   addstr_Uint (name, "   r = ", r);
   addstr_Ulong (name, ",   c = ", c);
   addstr_Uint (name, ",   w = ", w);
   addstr_ArrayUlong (name, ",   A = ", (int) r, A);
   addstr_ArrayUlong (name, ",   S = ", (int) r, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   Sum = 0.0;
   for (j = 0; j < r; j++) {
      if (A[j] >= num_TwoExp[w])
         util_Error ("ucarry_CreateMWCFloat:   A[i] must be < 2^w");
      if (S[j] >= num_TwoExp[w])
         util_Error ("ucarry_CreateMWCFloat:   S[i] must be < 2^w");
      Sum += A[j];
   }
   Sum = Sum * (num_TwoExp[w] - 1.0) + c;
   if (Sum >= num_TwoExp[53])
      util_Error
         ("ucarry_CreateMWCFloat:   c + (2^w - 1) * Sum A[i] >= 2^(53)");
   if (Sum >= num_TwoExp[32 + w])
      util_Error
         ("ucarry_CreateMWCFloat:   c + (2^w - 1) * Sum A[i] >= 2^(32 + w)");

   state->C = c;
   state->R = 0;
   state->Order = r;
   param->W = w;
   param->shift = 32 - w;
   if (w < 32) {
      param->Norm = 1.0 / num_TwoExp[w];
      param->mask = (unsigned long) num_TwoExp[w] - 1;
   } else {
      param->Norm = 1.0 / num_TwoExp[32];
      param->mask = 0xffffffffU;
   }
   for (j = 0; j < r; j++) {
      param->A[j] = A[j];
      state->X[j] = S[j];
   }

   gen->param = param;
   gen->state = state;
   gen->GetBits = &MWCFloat_Bits;
   gen->GetU01  = &MWCFloat_U01;
   gen->Write   = &WrMWCFloat;
   return gen;
}

void ucarry_DeleteMWCFloat (unif01_Gen * gen)
{
   MWCFloat_param *param;
   MWCFloat_state *state;

   if (NULL == gen)
      return;
   state = gen->state;
   param = gen->param;
   util_Free (state->X);
   util_Free (param->A);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}




/**************************************************************************/
#if 0
/* These are not used anywhere. Older versions?? (RS) */

#ifdef USE_LONGLONG

static double MWCfix8r4 (void)
{
   ulonglong S2 = UCA;
   S2 += 2736 * ((ulonglong) x1) + 7456 * ((ulonglong) x2)
      + 76438 * ((ulonglong) x5) + 84739 * ((ulonglong) x7);
   x1 = x2;
   x2 = x3;
   x3 = x4;
   x4 = x5;
   x5 = x6;
   x6 = x7;
   x7 = x8;
   x8 = S2 & e32m1;
   UCA = S2 >> 32;
   return (double) x8 *UnSur2e32;
}

static void WrMWCfix8 (void)
{
   printf
      ("   (x1, x2, x3, x4, x5, x6, x7, x8) = ( %lu, %lu, %lu,\n      %lu, %lu, %lu, %lu, %lu )\n   c  = %" PRIuLEAST64 "\n\n",
      x1, x2, x3, x4, x5, x6, x7, x8, UCA);
}

unif01_Gen *ucarry_CreateMWCfix8r4 (unsigned long c, unsigned long S[])
{
   strcpy (name, "ucarry_CreateMWCfix8r4: ");
   addstr_Ulong (name, " c = ", c);
   addstr_ArrayUlong (name, ",   S = ", 8, S);
   UCA = c;
   x1 = S[0];
   x2 = S[1];
   x3 = S[2];
   x4 = S[3];
   x5 = S[4];
   x6 = S[5];
   x7 = S[6];
   x8 = S[7];
   unif01_Gen = MWCfix8r4;
   unif01_WrGen = WrMWCfix8;
}

/**************************************************************************/

static double MWCfix8r8 (void)
{
   unsigned long S1;
   ulonglong S2 = UCA;
   S2 += 2736 * ((ulonglong) x1) + 7456 * ((ulonglong) x2)
      + 76438 * ((ulonglong) x5) + 84739 * ((ulonglong) x7)
      + 2736 * ((ulonglong) x3) + 7456 * ((ulonglong) x4)
      + 76438 * ((ulonglong) x6) + 84739 * ((ulonglong) x8);
   S1 = S2 & e32m1;
   UCA = S2 >> 32;
   x1 = x2;
   x2 = x3;
   x3 = x4;
   x4 = x5;
   x5 = x6;
   x6 = x7;
   x7 = x8;
   x8 = S1;
   return ((double) S1 * UnSur2e32);
}

unif01_Gen *ucarry_CreateMWCfix8r8 (unsigned long c, unsigned long S[])
{
   strcpy (name, "ucarry_CreateMWCfix8r8: ");
   addstr_Ulong (name, " c = ", c);
   addstr_ArrayUlong (name, ",   S = ", 8, S);
   UCA = c;
   x1 = S[0];
   x2 = S[1];
   x3 = S[2];
   x4 = S[3];
   x5 = S[4];
   x6 = S[5];
   x7 = S[6];
   x8 = S[7];
   unif01_Gen = MWCfix8r8;
   unif01_WrGen = WrMWCfix8;
}

#endif                            /* USE_LONGLONG */

#endif

/**************************************************************************/

static unsigned long MWCfixCouture_Bits (void *junk, void *vsta)
{
   MWCfixCouture_state *state = vsta;
   unsigned int s1, s3;
   unsigned int s2 = state->C;

   s2 += 14 * state->x1 + 18 * state->x2 + 144 * state->x3 + 1499 * state->x4
         + 2083 * state->x5 + 5273 * state->x6 + 10550 * state->x7 +
           45539 * state->x8;
   s1 = s2 & MASK16;
   s3 = s2 >> 16;
   s3 += 14 * state->x2 + 18 * state->x3 + 144 * state->x4 + 1499 * state->x5
      + 2083 * state->x6 + 5273 * state->x7 + 10550 * state->x8 + 45539 * s1;
   state->C = s3 >> 16;
   s3 &= MASK16;
   state->x1 = state->x3;
   state->x2 = state->x4;
   state->x3 = state->x5;
   state->x4 = state->x6;
   state->x5 = state->x7;
   state->x6 = state->x8;
   state->x7 = s1;
   state->x8 = s3;
   return s1 + (s3 << 16);
}

static double MWCfixCouture_U01 (void *vpar, void *vsta)
{
   return (MWCfixCouture_Bits (vpar, vsta) * unif01_INV32);
}

static void WrMWCfixCouture (void *vsta)
{
   MWCfixCouture_state *state = vsta;
   printf
      (" (x1, x2, x3, x4, x5, x6, x7, x8) = ( %u, %u, %u,\n      %u, %u, %u, %u, %u )\n c  = %u\n\n",
      state->x1, state->x2, state->x3, state->x4, state->x5, state->x6,
      state->x7, state->x8, state->C);
}

unif01_Gen *ucarry_CreateMWCfixCouture (unsigned int c, unsigned int S[])
{
   unif01_Gen *gen;
   MWCfixCouture_state *state;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (MWCfixCouture_state));

   strncpy (name, "ucarry_CreateMWCfixCouture:", (size_t) LEN);
   addstr_Uint (name, "   c = ", c);
   addstr_ArrayUint (name, ",   S = ", 8, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->C = c;
   state->x1 = S[0];
   state->x2 = S[1];
   state->x3 = S[2];
   state->x4 = S[3];
   state->x5 = S[4];
   state->x6 = S[5];
   state->x7 = S[6];
   state->x8 = S[7];

   gen->param = NULL;
   gen->state = state;
   gen->GetBits = &MWCfixCouture_Bits;
   gen->GetU01  = &MWCfixCouture_U01;
   gen->Write   = &WrMWCfixCouture;
   return gen;
}

void ucarry_DeleteMWCfixCouture (unif01_Gen *gen)
{
   if (NULL == gen) return;
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/****************************************************************************/

static unsigned long SWC_Bits (void *vpar, void *vsta)
{
/*
 * Generator Shift With Carry on 32 bits
 */
   SWC_param *param = vpar;
   SWC_state *state = vsta;
   unsigned int Som;
   unsigned int j;

   Som = state->C;
   state->C = 0;
   for (j = 0; j <= param->Asize / 2 - 1; j++) {
      state->S = (state->R + param->A[2 * j]) % state->Order;
      Som ^= state->X[state->S] << param->A[2 * j + 1];
      state->C ^= state->X[state->S] >> (32 - param->A[2 * j + 1]);
   }

   if (param->W != 32)
      state->C = (state->C << (32 - param->W)) | (Som >> param->W);

   Som &= param->mask;
   state->X[state->R] = Som;
   state->R = (state->R + 1) % state->Order;
   return Som << param->shift;
}

static double SWC_U01 (void *vpar, void *vsta)
{
   return (SWC_Bits (vpar, vsta) * unif01_INV32);
}

static unsigned long SWCshort_Bits (void *vpar, void *vsta)
{
/*
 * Generator Shift With Carry on 16 bits
 */
   SWC_param *param = vpar;
   SWC_state *state = vsta;
   unsigned int Som;
   unsigned int j;

   Som = state->C;
   for (j = 0; j <= param->Asize / 2 - 1; j++) {
      state->S = (state->R + param->A[2 * j]) % state->Order;
      Som ^= (state->X[state->S] << param->A[2 * j + 1]);
   }

   state->C = Som >> param->W;
   Som &= param->mask;
   state->X[state->R] = Som;
   state->R = (state->R + 1) % state->Order;
   return Som << param->shift;
}

static double SWCshort_U01 (void *vpar, void *vsta)
{
   return (SWCshort_Bits (vpar, vsta) * unif01_INV32);
}

static void WrSWC (void *vsta)
{
   SWC_state *state = vsta;
   unsigned int j;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < state->Order - 1; j++) {
	 printf ("%10u,  ", state->X[j]);
	 if (((j + 1) % 5) == 0)
	    printf ("\n ");
      }
      printf ("%10u   }\n\n", state->X[state->Order - 1]);
      printf (" c = %1u\n\n", state->C);
   } else
      unif01_WrLongStateDef ();
}


unif01_Gen *ucarry_CreateSWC (unsigned int r, unsigned int h,
   unsigned int c, unsigned int w, unsigned int A[], unsigned int S[])
{
   unif01_Gen *gen;
   SWC_param *param;
   SWC_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int j;

   if ((w > 32U) || (h > 31U * r)) {
      util_Error ("ucarry_CreateSWC:   invalid parameter");
   }
   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (SWC_param));
   state = util_Malloc (sizeof (SWC_state));
   state->X = util_Calloc ((size_t) r, sizeof (unsigned int));
   param->A = util_Calloc ((size_t) h, sizeof (unsigned int));

   strncpy (name, "ucarry_CreateSWC:", (size_t) LEN);
   addstr_Uint (name, "   r = ", r);
   addstr_Uint (name, ",   h = ", h);
   addstr_Uint (name, ",   c = ", c);
   addstr_Uint (name, ",   w = ", w);
   addstr_ArrayUint (name, ",   A = ", (int) h, A);
   addstr_ArrayUint (name, ",   S = ", (int) r, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->R = 0;
   state->S = r;
   state->Order = r;
   state->C = c;
   param->Asize = h;
   param->W = w;
   param->shift = 32 - w;
   if (w < 32)
      param->mask = num_TwoExp[w] - 1;
   else
      param->mask = 0xffffffffU;               /* 2^32 - 1 */

   for (j = 0; j < h; j++)
      param->A[j] = A[j] & 0xff;
   for (j = 0; j < r; j++)
      state->X[j] = S[j] & param->mask;

   if (w > 15) {
      gen->GetBits = &SWC_Bits;
      gen->GetU01 = &SWC_U01;
   } else {
      gen->GetBits = &SWCshort_Bits;
      gen->GetU01 = &SWCshort_U01;
   }

   gen->param = param;
   gen->state = state;
   gen->Write = &WrSWC;
   return gen;
}

void ucarry_DeleteSWC (unif01_Gen * gen)
{
   SWC_param *param;
   SWC_state *state;

   if (NULL == gen)
      return;
   state = gen->state;
   param = gen->param;
   util_Free (state->X);
   util_Free (param->A);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/***************************************************************************/

static unsigned long MWC1616_Bits (void *par, void *sta)
{
   MWC1616_state *state = sta;
   MWC1616_param *param = par;
   state->x = param->a * (state->x & MASK16) + (state->x >> 16);
   state->y = param->b * (state->y & MASK16) + (state->y >> 16);
   return (state->x << 16) + (state->y & MASK16);
}

/*-----------------------------------------------------------------------*/

static double MWC1616_U01 (void *par, void *sta)
{
   return MWC1616_Bits (par, sta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrMWC1616 (void *sta)
{
   MWC1616_state *state = sta;
   printf (" x = %u,   y = %u\n", state->x, state->y);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ucarry_CreateMWC1616 (unsigned int a, unsigned int b,
   unsigned int x, unsigned int y)
{
   unif01_Gen *gen;
   MWC1616_state *state;
   MWC1616_param *param;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   gen->state = state = util_Malloc (sizeof (MWC1616_state));
   gen->param = param = util_Malloc (sizeof (MWC1616_param));
   state->x = x & MASK32;
   state->y = y & MASK32;
   param->a = a & MASK16;
   param->b = b & MASK16;
   gen->Write = WrMWC1616;
   gen->GetU01 = MWC1616_U01;
   gen->GetBits = MWC1616_Bits;

   strcpy (name, "ucarry_CreateMWC1616:");
   addstr_Uint (name, "   a = ", a);
   addstr_Uint (name, ",   b = ", b);
   addstr_Uint (name, ",   x = ", x);
   addstr_Uint (name, ",   y = ", y);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);
   return gen;
}

/***************************************************************************/

void ucarry_DeleteGen (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}
