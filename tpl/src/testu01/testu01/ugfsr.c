/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ugfsr.c
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

#include "ugfsr.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>





/*=============================== Macros ===============================*/

#define LEN  300                     /* Max length of strings */

#define MATRIX_A 0x9908b0dfUL        /* constant vector a for MT19937 */
#define UPPER_MASK 0x80000000UL      /* most significant w - r bits   */
#define LOWER_MASK 0x7fffffff        /* least significant r bits      */

/* for tempering */
#define TEMPERING_MASK_B 0x9d2c5680UL
#define TEMPERING_MASK_C 0xefc60000UL
#define TEMPERING_SHIFT_U(y)  ((y) >> 11)
#define TEMPERING_SHIFT_S(y)  ((y) << 7)
#define TEMPERING_SHIFT_T(y)  ((y) << 15)
#define TEMPERING_SHIFT_L(y)  ((y) >> 18)






/*================================ Types ================================*/

typedef struct {
   unsigned long Shift, Mask, mag01[2];
   double Norm;
} GFSR_param;

typedef struct {
   unsigned long *X;
   unsigned int r, s, K; 
} GFSR_state;

/*------------------------------------*/

typedef struct {
   unsigned long Shift, Mask, mag01[2];
   unsigned long b, c, S2, T2;
} TGFSR2_param;

typedef GFSR_state TGFSR2_state;

/*------------------------------------*/

typedef struct {
   unsigned long mag01[2];
} MT19937_param;

typedef GFSR_state MT19937_state;

/*------------------------------------*/

typedef struct {
   unsigned long Shift;
} GFSR5_param;

typedef struct {
   unsigned long *X;
   unsigned int r1, r2, r3, s, K; 
} GFSR5_state;

/*------------------------------------*/






/***************************************************************************/

static void PRODUIT (int k, unsigned int F[], unsigned int G[],
   unsigned int H[])
{
   int i, j;
   unsigned int Work;

   for (i = 0; i < k; i++) {
      Work = 0;
      for (j = 0; j < k; j++) {
         if (G[j] == 1)
            Work += F[i + j];
      }
      H[i] = Work % 2;
   }
}


/***************************************************************************/

#define Fushia 69069              /* multiplier */

static void InitFushimi (int k, int r, int s, GFSR_state *state)
/*
 * Initialization copied from Fushimi90
 */
{
   unsigned int *x, *x0, *x1, *x2, *x3;
   unsigned int E0, E1;
   int wi, ix, B[32];
   int i, iMD2, k2, k3, j;

   E0 = 0;
   E1 = 1;
   k2 = k * 2;
   k3 = k * 3;
   state->K = k3;
   state->r = 3 * (k - r);
   state->s = 0;

   x = calloc ((size_t) 3 * (k + 1), sizeof (unsigned int));
   x0 = calloc ((size_t) 2 * (k + 1), sizeof (unsigned int));
   x1 = calloc ((size_t) 2 * (k + 1), sizeof (unsigned int));
   x2 = calloc ((size_t) 1 * (k + 1), sizeof (unsigned int));
   x3 = calloc ((size_t) 3 * (k + 1), sizeof (unsigned int));

   /* Initialization */
   for (j = 0; j < k2; j++) {
      x0[j] = E0;
      x3[j] = E0;
   }

   /* Table B contains 2^30, 2^29, ... , 1, 0. */
   B[31] = 0;
   B[30] = 1;
   for (j = 0; j <= 29; j++)
      B[29 - j] = B[30 - j] * 2;

   ix = s;
   for (j = 0; j < k; j++) {
      ix *= Fushia;
      if (ix > 0)
         x0[j] = E1;
   }
   for (j = k; j < k2; j++)
      x0[j] = x0[j - k] ^ x0[j - r];
   x3[1] = E1;
   for (i = 0; i <= k - 2; i++) {
      iMD2 = i % 2;
      for (j = k - 1; j >= 0; j--) {
         x3[2 * j + iMD2] = x3[j];
         x3[2 * j + 1 - iMD2] = E0;
      }
      for (j = k - 1; j >= 0; j--) {
         x3[j] ^= x3[j + k];
         x3[j + k - r] ^= x3[j + k];
         x3[j + k] = E0;
      }
   }
   PRODUIT (k, x0, x3, x1);
   for (j = k; j < k2; j++)
      x1[j] = x1[j - k] ^ x1[j - r];
   PRODUIT (k, x1, x3, x2);
   for (j = 0; j <= k; j++) {
      x[3 * j] = x0[j];
      x[3 * j + 1] = x1[j];
      x[3 * j + 2] = x2[j];
   }
   for (i = 0; i < k3; i++) {
      wi = 0;
      for (j = 0; j < 32; j++) {
         if (x[state->s] != E0)
            wi += B[j];
         x[state->s] ^= x[state->r];
         state->s++;
         if (state->s == state->K)
            state->s = 0;
         state->r++;
         if (state->r == state->K)
            state->r = 0;
      }
      state->X[i] = wi & unif01_MASK32;
   }
   free (x);
   free (x0);
   free (x1);
   free (x2);
   free (x3);
}


/***************************************************************************/

static unsigned long GFSR_Bits (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;

   if (++state->s == state->K)
      state->s = 0;                     /* s = (s + 1) % K */
   if (++state->r == state->K)
      state->r = 0;                     /* r = (r + 1) % K */
   state->X[state->s] ^= state->X[state->r];
   return state->X[state->s] << param->Shift;
}

/*-----------------------------------------------------------------------*/

static double GFSR_U01 (void *vpar, void *vsta)
{
   return GFSR_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrGFSR (void *vsta)
{
   GFSR_state *state = vsta;
   unsigned int j;
   unsigned int s = state->s;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < state->K; j++) {
         if (++s >= state->K)
            s = 0;
         printf (" %12lu", state->X[s]);
         if (j < state->K - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

static unif01_Gen *CreateGFSR0 (unsigned int k, unsigned int r,
   unsigned int l, unsigned long S[], const char nom[])
/*
 * CreateGen common to many generators GFSR and TGFSR
 */
{
   unif01_Gen *gen;
   GFSR_param *param;
   GFSR_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int j;
   unsigned long Mask;

   strcpy (name, nom);
   addstr_Uint (name, "   k = ", k);
   addstr_Uint (name, ",   r = ", r);
   addstr_Uint (name, ",   l = ", l);
   addstr_ArrayUlong (name, ",   S = ", (int) k, S);
   if ((l < 1) || (r >= k) || (l > 32))
      util_Error (name);

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (GFSR_param));
   state = util_Malloc (sizeof (GFSR_state));

   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (l == 32)
      Mask = unif01_MASK32;
   else
      Mask = num_TwoExp[l] - 1.0;
   state->X = util_Calloc ((size_t) k, sizeof (unsigned long));
   for (j = 0; j < k; j++)
      state->X[j] = S[j] & Mask;

   state->s = 0;
   state->r = k - r;
   state->K = k;
   param->Shift = 32 - l;
   param->Mask = Mask;

   gen->param = param;
   gen->state = state;
   gen->GetBits = &GFSR_Bits;
   gen->GetU01  = &GFSR_U01;
   gen->Write   = &WrGFSR;
   return gen;
}

/*=========================================================================*/

unif01_Gen *ugfsr_CreateGFSR3 (unsigned int k, unsigned int r,
   unsigned int l, unsigned long S[])
{
   return CreateGFSR0 (k, r, l, S, "ugfsr_CreateGFSR3:");
}


/*=========================================================================*/

#define RIPa      16807           /* multiplier */
#define RIPk        521           /* Parametres for */
#define RIPr         32           /* Ripley90 ....  */
#define RIPk2  (2*RIPk)
#define RIPc     127773


static double Ripley90_U01 (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long temp;

   state->s--;
   state->r--;
   temp = state->X[state->r];
   state->X[state->r] ^= state->X[state->s];
   if (state->s == 0)
      state->s = RIPk;
   if (state->r == 0)
      state->r = RIPk;
   return temp * param->Norm;
}

/*-----------------------------------------------------------------------*/

static unsigned long Ripley90_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (Ripley90_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrRipley90 (void *vsta)
{
   GFSR_state *state = vsta;
   int j;
   int r = state->r;

   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < RIPk; j++) {
         r--;
         printf (" %12lu", state->X[r]);
         if (r <= 0)
            r = RIPk;
         if (j < RIPk - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ugfsr_CreateRipley90 (long s)
{
   unif01_Gen *gen;
   GFSR_param *param;
   GFSR_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned long x0[RIPk2];
   unsigned long wi, wt;
   long vt, ix, i, j;                /* Work variables */

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (GFSR_param));
   state = util_Malloc (sizeof (GFSR_state));

   strcpy (name, "ugfsr_CreateRipley90:");
   addstr_Long (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->K = RIPk;
   state->r = RIPk - RIPr;
   state->s = RIPk;
   param->Norm = 1.0 / (num_TwoExp[31] - 1.0);
   state->X = util_Calloc ((size_t) RIPk, sizeof (unsigned long));

   /* Initialization */
   ix = s;
   for (j = 0; j < RIPk; j++) {
      x0[j] = 0;
      vt = ix / RIPc;
      ix = RIPa * (ix - RIPc * vt) - 2836 * vt;
      if (ix < 0)
         ix += 2147483647;
      if (ix > 1073741824)
         x0[j] = 1;
   }
   for (j = RIPk; j < RIPk2; j++)
      x0[j] = x0[j - RIPk] ^ x0[j - RIPr];
   for (i = 0; i < RIPk; i++) {
      wi = 0;
      for (j = 0; j <= 30; j++) {
         wt = x0[i + 16 * (j + 1)];
         wt <<= j;
#ifndef IS_ULONG32
         wt &= unif01_MASK32;
#endif
         wi += wt;
#ifndef IS_ULONG32
         wi &= unif01_MASK32;
#endif
      }
      state->X[i] = wi & unif01_MASK32;
   }
   gen->param   = param;
   gen->state   = state;
   gen->GetBits = &Ripley90_Bits;
   gen->GetU01  = &Ripley90_U01;
   gen->Write   = &WrRipley90;
   return gen;
}


/***************************************************************************/

unif01_Gen * ugfsr_CreateToot73 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;
   GFSR_state *state;
   size_t leng;
   char name[LEN + 1];
   long q, l, m, n;
   int h, t, j, i, K = 607, R = 273;
   unsigned long Q[700], Mask;

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (GFSR_param));
   state = util_Malloc (sizeof (GFSR_state));

   strcpy (name, "ugfsr_CreateToot73:");
   addstr_ArrayUlong (name, "   S = ", K, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->X = util_Calloc ((size_t) 700, sizeof (unsigned long));
   state->r = K - R;
   state->s = 0;
   state->K = K;
   Mask = num_TwoExp[23] - 1.0;
   param->Shift = 9;
   param->Mask = Mask;

   q = S[19];
   m = S[10];
   for (j = 1; j <= 19; j++)
      Q[j] = S[j];
   for (j = 19; j <= K; j++) {
      l = (long) Q[j - 18];
      n = (long) Q[j - 8];
      Q[j] = (((unsigned long) q) << 1) + (((unsigned long) l) >> 31);
      Q[j] ^= (((unsigned long) m) << 15) + (((unsigned long) n) >> 17);
      Q[j] &= unif01_MASK32;
      q = l;
      m = n;
      Q[j - 18] = l & Mask;
   }

   for (j = 590; j <= K; j++)
      Q[j] &= Mask;
   i = 1;
   t = 0;
   do {
      for (j = i; j <= K; j += 16) {
         state->X[t] = Q[j];
         t++;
      }
      i++;
      /* Recompute table Q so that the equation */
      /* Q(n) =( ( x^15*Q(n- 9)+x^- 17*Q(n- 8) ) XOR */
      /* ( x*Q(n-19)+x^- 31*Q(n-18) ) */
      /* remains valid */
      for (h = 1; h <= R; h++)
         Q[h] ^= Q[h + K - R];
      for (h = R + 1; h <= K; h++)
         Q[h] ^= Q[h - R];
   } while (t <= K);

   gen->param = param;
   gen->state = state;
   gen->GetBits = &GFSR_Bits;
   gen->GetU01  = &GFSR_U01;
   gen->Write   = &WrGFSR;
   return gen;
}


/**************************************************************************/

static unsigned long Fushimi90_Bits (void *junk, void *vsta)
{
   GFSR_state *state = vsta;
   unsigned long V;

   V = state->X[state->s];
   state->X[state->s] ^= state->X[state->r];
   if (++state->s == state->K)
      state->s = 0;
   if (++state->r == state->K)
      state->r = 0;
   return V << 1;
}

/*-----------------------------------------------------------------------*/

static double Fushimi90_U01 (void *vpar, void *vsta)
{
   return Fushimi90_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateFushimi90 (int s)
{
   unif01_Gen *gen;
   GFSR_state *state;
   GFSR_param *param;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (GFSR_param));
   state = util_Malloc (sizeof (GFSR_state));

   state->K = 521;
   state->r = 521 - 32;
   state->s = 0;

   strcpy (name, "ugfsr_CreateFushimi90:");
   addstr_Int (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->X = util_Calloc ((size_t) 3 * 522, sizeof (unsigned long));
   InitFushimi (521, 32, s, state);
   gen->param = param;
   gen->state = state;
   gen->GetBits = &Fushimi90_Bits;
   gen->GetU01  = &Fushimi90_U01;
   gen->Write   = &WrGFSR;
   return gen;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateFushimi (int k, int r, int s)
{
   unif01_Gen *gen;
   GFSR_state *state;
   GFSR_param *param;
   size_t leng;
   char name[LEN + 1];

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (GFSR_param));
   state = util_Malloc (sizeof (GFSR_state));

   state->K = k;
   state->r = k - r;
   state->s = 0;

   strcpy (name, "ugfsr_CreateFushimi:");
   addstr_Int (name, "   k = ", k);
   addstr_Int (name, ",   r = ", r);
   addstr_Int (name, ",   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->X = util_Calloc ((size_t) k * 3, sizeof (unsigned long));
   InitFushimi (k, r, s, state);
   gen->param = param;
   gen->state = state;
   gen->GetBits = &Fushimi90_Bits;
   gen->GetU01  = &Fushimi90_U01;
   gen->Write   = &WrGFSR;
   return gen;
}


/**************************************************************************/

unif01_Gen *ugfsr_CreateKirk81 (long s)
{
   unif01_Gen *gen;
   GFSR_state *state;
   GFSR_param *param;
   size_t leng;
   char name[LEN + 1];
   unsigned long mask, msb;
   unsigned int j, k;
   long vt;

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (GFSR_param));
   state = util_Malloc (sizeof (GFSR_state));

   strcpy (name, "ugfsr_CreateKirk81:");
   addstr_Long (name, "   s = ", s);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   state->K = 250;
   state->s = 0;
   state->r = 147;
   state->X = util_Calloc ((size_t) state->K, sizeof (unsigned long));

   msb = 2147483648UL;                    /* TwoExp31 */
   mask = 4294967295UL;                   /* TwoExp32 - 1 */
   for (j = 0; j < state->K; j++) {
      vt = s / 127773;
      s = 16807 * (s - 127773 * vt) - 2836 * vt;
      if (s < 0)
         s += 2147483647;
      state->X[j] = 2 * ((unsigned long) s);
      if (s > 1000000000)
         state->X[j] += 1;
      /* When s > 1000000000, then set last bit to 1 */
   }
   for (j = 1; j <= 31; j++) {
      k = 7 * j + 3;
      state->X[k] = (state->X[k] & mask) | msb;
      mask >>= 1;
      msb >>= 1;
   }
   param->Shift = 0;
   gen->param = param;
   gen->state = state;
   gen->GetBits = &GFSR_Bits;
   gen->GetU01  = &GFSR_U01;
   gen->Write   = &WrGFSR;
   return gen;
}


/**************************************************************************/

static unsigned long TGFSR_Bits (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long v;

   state->X[state->s] = state->X[state->r] ^ (state->X[state->s] >> 1)
      ^ param->mag01[state->X[state->s] % 2];
   v = state->X[state->s] & param->Mask;
   if (++state->s == state->K)
      state->s = 0;                     /* s = (s + 1) % K */
   if (++state->r == state->K)
      state->r = 0;                     /* r = (r + 1) % K */
   return v << param->Shift;
}

/*-----------------------------------------------------------------------*/

static double TGFSR_U01 (void *vpar, void *vsta)
{
   return TGFSR_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTGFSR (unsigned int k, unsigned int r,
   unsigned int l, unsigned long Av, unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;
   char name[LEN + 1] = "";
   size_t leng;

   gen = CreateGFSR0 (k, r, l, S, "ugfsr_CreateTGFSR:");
   addstr_Ulong (name, ",   Av = ", Av);
   leng = strlen (gen->name) + strlen (name);
   gen->name = util_Realloc (gen->name, (leng + 1) * sizeof (char));
   strncat (gen->name, name, leng);

   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = Av;
   gen->GetBits = &TGFSR_Bits;
   gen->GetU01  = &TGFSR_U01;
   return gen;
}


/**************************************************************************/

static unsigned long T800_Bits (void *vpar, void *vsta)
{
   GFSR_state *state = vsta;
   GFSR_param *param = vpar;
   unsigned long v;

   state->X[state->s] = state->X[state->r] ^ (state->X[state->s] >> 1)
      ^ param->mag01[state->X[state->s] % 2];
   v = state->X[state->s];
#ifndef IS_ULONG32
   v &= 0xffffffffUL;        /* you may delete this line if word size = 32 */
#endif
   if (++state->s == state->K)
      state->s = 0;                     /* s = (s + 1) % K */
   if (++state->r == state->K)
      state->r = 0;                     /* r = (r + 1) % K */
   return v;
}

/*-----------------------------------------------------------------------*/

static double T800_U01 (void *vpar, void *vsta)
{
   return T800_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateT800 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;
   char name[LEN + 1] = "";
   size_t leng;

   gen = CreateGFSR0 (25, 18, 32, S, "ugfsr_CreateT800:");
   addstr_Ulong (name, ",   Av = ", 2394935336UL);
   leng = strlen (gen->name) + strlen (name);
   gen->name = util_Realloc (gen->name, (leng + 1) * sizeof (char));
   strncat (gen->name, name, leng);

   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = 0x8ebfd028UL;       /* this is magic vector 'a' */
   gen->GetBits = &T800_Bits;
   gen->GetU01  = &T800_U01;
   return gen;
}


/***************************************************************************/

static unsigned long TGFSR2_Bits (void *vpar, void *vsta)
{
   TGFSR2_param *param = vpar;
   TGFSR2_state *state = vsta;
   unsigned long v;

   v = state->X[state->s];
   state->X[state->s] =
      state->X[state->r] ^ (v >> 1) ^ param->mag01[v % 2];
   if (++state->s == state->K)
      state->s = 0;                      /* s = (s + 1) % K */
   if (++state->r == state->K)
      state->r = 0;                      /* r = (r + 1) % K */
   v ^= (v << param->S2) & param->b;     /* s and b, magic vectors */
   v ^= (v << param->T2) & param->c;     /* t and c, magic vectors */
   v &= param->Mask;                     /* l least significant bits */
   return v << param->Shift;
}

/*-----------------------------------------------------------------------*/

static double TGFSR2_U01 (void *vpar, void *vsta)
{
   return TGFSR2_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTGFSR2 (unsigned int k, unsigned int r,
   unsigned int l, unsigned int s, unsigned int t, unsigned long Av,
   unsigned long Bv, unsigned long Cv, unsigned long S[])
{
   unif01_Gen *gen;
   TGFSR2_param *param;
   char name[LEN + 1];
   size_t leng;

   gen = CreateGFSR0 (k, r, l, S, "");
   util_Free (gen->name);
   strcpy (name, "ugfsr_CreateTGFSR2:");
   addstr_Uint (name, "   k = ", k);
   addstr_Uint (name, ",   r = ", r);
   addstr_Uint (name, ",   l = ", l);
   addstr_Ulong (name, ",   Av = ", Av);
   addstr_Ulong (name, ",   Bv = ", Bv);
   addstr_Ulong (name, ",   Cv = ", Cv);
   addstr_Uint (name, ",   s = ", s);
   addstr_Uint (name, ",   t = ", t);
   addstr_ArrayUlong (name, ",   S", (int) k, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   util_Free (gen->param);
   gen->param = param = util_Malloc (sizeof (TGFSR2_param));
   param->S2 = s;
   param->T2 = t;
   param->mag01[0] = 0x0;
   param->mag01[1] = Av;                 /* this is magic vector 'a' */
   param->b = Bv;
   param->c = Cv;
   param->Shift = 32 - l;
   param->Mask = num_TwoExp[l] - 1.0;
   if (l == 32)
      param->Mask = unif01_MASK32;
   gen->GetBits = &TGFSR2_Bits;
   gen->GetU01  = &TGFSR2_U01;
   gen->Write   = &WrGFSR;
   return gen;
}


/***************************************************************************/

static double TT400_U01 (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long v;

   v = state->X[state->s];
   state->X[state->s] = state->X[state->r] ^ (v >> 1) ^ param->mag01[v % 2];
   v ^= (v << 2) & 0x6A68;        /* s and b, magic vectors */
   v ^= (v << 7) & 0x7500;        /* t and c, magic vectors */
   v &= 0xffff;
   if (++state->s == state->K)
      state->s = 0;
   if (++state->r == state->K)
      state->r = 0;
   return (double) v / 0xffff;
}

/*-----------------------------------------------------------------------*/

static unsigned long TT400_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (TT400_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTT400 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;
   unsigned int N, M;

   N = 25;
   M = 25 - 11;
   gen = CreateGFSR0 (N, M, 16, S, "ugfsr_CreateTT400:");

   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = 0xA875;             /* this is magic vector 'a' */
   gen->GetBits = &TT400_Bits;
   gen->GetU01  = &TT400_U01;
   return gen;
}


/***************************************************************************/

static double TT403_U01 (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long v;

   v = state->X[state->s];
   state->X[state->s] = state->X[state->r] ^ (v >> 1) ^ param->mag01[v % 2];
   v ^= (v << 8) & 0x102D1200;          /* s and b, magic vectors */
   v ^= (v << 14) & 0x66E50000;         /* t and c, magic vectors */
   v &= 0x7fffffff;
   if (++state->s == state->K)
      state->s = 0;
   if (++state->r == state->K)
      state->r = 0;
   return (double) v / 0x7fffffff;
}

/*-----------------------------------------------------------------------*/

static unsigned long TT403_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (TT403_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTT403 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;
   unsigned int N, M;

   N = 13;
   M = 13 - 2;
   gen = CreateGFSR0 (N, M, 31, S, "ugfsr_CreateTT403:");

   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = 0x6B5ECCF6;         /* this is magic vector 'a' */
   gen->GetBits = &TT403_Bits;
   gen->GetU01  = &TT403_U01;
   return gen;
}


/***************************************************************************/

static double TT775_U01 (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long v;

   v = state->X[state->s];
   state->X[state->s] = state->X[state->r] ^ (v >> 1) ^ param->mag01[v % 2];
   v ^= (v << 6) & 0x1ABD5900;            /* s and b, magic vectors */
   v ^= (v << 14) & 0x776A0000;           /* t and c, magic vectors */
   v &= 0x7fffffff;
   if (++state->s == state->K)
      state->s = 0;
   if (++state->r == state->K)
      state->r = 0;
   return (double) v / 0x7fffffff;
}

/*-----------------------------------------------------------------------*/

static unsigned long TT775_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (TT775_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTT775 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;
   unsigned int N, M;

   N = 25;
   M = 25 - 8;
   gen = CreateGFSR0 (N, M, 31, S, "ugfsr_CreateTT775:");

   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = 0x6C6CB38C;         /* this is magic vector 'a' */
   gen->GetBits = &TT775_Bits;
   gen->GetU01  = &TT775_U01;
   return gen;
}


/**************************************************************************
 *
 * Our version of TT800
 *
 **************************************************************************/

static unsigned long TT800_Bits (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long v;

   v = state->X[state->s];
   state->X[state->s] = state->X[state->r] ^ (v >> 1) ^ param->mag01[v % 2];
   v ^= (v << 7) & 0x2b5b2500;          /* s and b, magic vectors */
   v ^= (v << 15) & 0xdb8b0000UL;       /* t and c, magic vectors */
#ifndef IS_ULONG32
   v &= 0xffffffffUL;         /* you may delete this line if word size = 32 */
#endif
   if (++state->s == state->K)
      state->s = 0;
   if (++state->r == state->K)
      state->r = 0;
   return v;
}

/*-----------------------------------------------------------------------*/

static double TT800_U01 (void *vpar, void *vsta)
{
   return TT800_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTT800 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;
   unsigned int N, M;

   N = 25;
   M = 25 - 7;
   gen = CreateGFSR0 (N, M, 32, S, "ugfsr_CreateTT800:");

   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = 0x8ebfd028UL;       /* this is magic vector 'a' */
   gen->GetBits = &TT800_Bits;
   gen->GetU01  = &TT800_U01;
   return gen;
}


/***************************************************************************/

/* A C-program for TT800 : July 8th 1994 Version                        */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp                     */
/* genrand() generate one pseudorandom number with double precision     */
/* which is uniformly distributed on [0,1]-interval                     */
/* for each call.  One may choose any initial 25 seeds                  */
/* except all zeros.                                                    */

#define NN 25
#define MM 7

static double TT800M94_U01 (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long y;
   unsigned int j;

   if (state->s == NN) {
      /* generate NN words at one time */
      for (j = 0; j < NN - MM; j++) {
         state->X[j] = state->X[j + MM] ^ (state->X[j] >> 1) ^
            param->mag01[state->X[j] % 2];
      }
      for (; j < NN; j++) {
         state->X[j] = state->X[j + (MM - NN)] ^ (state->X[j] >> 1) ^
            param->mag01[state->X[j] % 2];
      }
      state->s = 0;
   }
   y = state->X[state->s++];
   y ^= (y << 7) & 0x2b5b2500;       /* s and b, magic vectors */
   y ^= (y << 15) & 0xdb8b0000UL;    /* t and c, magic vectors */
#ifndef IS_ULONG32
   y &= 0xffffffffUL;      /* you may delete this line if word size = 32 */
#endif

   return (double) y / 0xffffffffUL;
}

/*-----------------------------------------------------------------------*/

static unsigned long TT800M94_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (TT800M94_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/


/* A C-program for TT800 : July 8th 1996 Version                        */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp                     */
/* genrand() generate one pseudorandom number with double precision     */
/* which is uniformly distributed on [0,1]-interval                     */
/* for each call.  One may choose any initial 25 seeds                  */
/* except all zeros.                                                    */

static double TT800M96_U01 (void *vpar, void *vsta)
{
   GFSR_param *param = vpar;
   GFSR_state *state = vsta;
   unsigned long y;
   unsigned int j;

   if (state->s == NN) {
      /* generate NN words at one time */
      for (j = 0; j < NN - MM; j++) {
         state->X[j] = state->X[j + MM] ^ (state->X[j] >> 1) ^
            param->mag01[state->X[j] % 2];
      }
      for (; j < NN; j++) {
         state->X[j] = state->X[j + (MM - NN)] ^ (state->X[j] >> 1) ^
            param->mag01[state->X[j] % 2];
      }
      state->s = 0;
   }
   y = state->X[state->s++];
   y ^= (y << 7) & 0x2b5b2500;       /* s and b, magic vectors */
   y ^= (y << 15) & 0xdb8b0000UL;    /* t and c, magic vectors */
#ifndef IS_ULONG32
   y &= 0xffffffffUL;                /* you may delete this line if word
                                        size = 32 */
#endif

   /* the following line was added by Makoto Matsumoto in the 1996 version
      to improve lower bit's correlation. This line differs from the code
      published in 1994.  */
   y ^= (y >> 16);                /* added to the 1994 version */

   return (double) y / 0xffffffffUL;
}

/*-----------------------------------------------------------------------*/

static unsigned long TT800M96_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (TT800M96_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTT800M94 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;

   gen = CreateGFSR0 (NN, MM, 32, S, "ugfsr_CreateTT800M94:");
   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = 0x8ebfd028UL;       /* this is magic vector 'a' */
   gen->GetBits = &TT800M94_Bits;
   gen->GetU01  = &TT800M94_U01;
   return gen;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateTT800M96 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_param *param;

   gen = CreateGFSR0 (NN, MM, 32, S, "ugfsr_CreateTT800M96:");
   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = 0x8ebfd028UL;       /* this is magic vector 'a' */
   gen->GetBits = &TT800M96_Bits;
   gen->GetU01  = &TT800M96_U01;
   return gen;
}


/**************************************************************************
 *
 * The code for the two versions of MT19937 is copyrighted by Makoto
 * Matsumoto and Takuji Nishimura. I have changed the name of some variables
 * to make it consistent with TestU01. (R. Simard)
 */

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

/* Period parameters */
#undef NN
#undef MM
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0dfUL     /* constant vector a */
#define UPPER_MASK 0x80000000UL   /* most significant w-r bits */
#undef LOWER_MASK
#define LOWER_MASK 0x7fffffffUL   /* least significant r bits */


/*-----------------------------------------------------------------------*/

/* initializes mt[NN] with a seed */
static void init_genrand (void *vsta, unsigned long s)
{
   MT19937_state *state = vsta;
   int j;

   state->X[0] = s & 0xffffffffUL;
   for (j = 1; j < NN; j++) {
      state->X[j] =
         (1812433253UL * (state->X[j - 1] ^ (state->X[j - 1] >> 30)) + j);
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array state->X[].                  */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      state->X[j] &= 0xffffffffUL;
      /* for >32 bit machines */
   }
   state->s = NN;
}


/*-----------------------------------------------------------------------*/

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
static void init_by_array (void *vsta, unsigned long init_key[],
   int key_length)
{
   MT19937_state *state = vsta;
   int i, j, k;
   init_genrand (vsta, 19650218UL);
   i = 1;
   j = 0;
   k = (NN > key_length ? NN : key_length);
   for (; k; k--) {
      state->X[i] = (state->X[i] ^ ((state->X[i - 1] ^ (state->X[i - 1] >>
                     30)) * 1664525UL)) + init_key[j] + j;  /* non linear */
      state->X[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      j++;
      if (i >= NN) {
         state->X[0] = state->X[NN - 1];
         i = 1;
      }
      if (j >= key_length)
         j = 0;
   }
   for (k = NN - 1; k; k--) {
      state->X[i] = (state->X[i] ^ ((state->X[i - 1] ^ (state->X[i - 1] >>
                     30)) * 1566083941UL)) - i;          /* non linear */
      state->X[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if (i >= NN) {
         state->X[0] = state->X[NN - 1];
         i = 1;
      }
   }

   state->X[0] = 0x80000000UL;    /* MSB is 1; assuring non-zero initial
                                     array */
}


/*-----------------------------------------------------------------------*/

/* generates a random number on [0,0xffffffff]-interval
unsigned long genrand_int32(void) */
static unsigned long MT19937_02_Bits (void *vpar, void *vsta)
{
   MT19937_param *param = vpar;
   MT19937_state *state = vsta;
   unsigned long y;

   if (state->s >= NN) {          /* generate NN words at one time */
      int kk;

      if (state->s == NN + 1)     /* if init_genrand() has not been called, */
         init_genrand (state, 5489UL);    /* a default initial seed is used */

      for (kk = 0; kk < NN - MM; kk++) {
         y = (state->X[kk] & UPPER_MASK) | (state->X[kk + 1] & LOWER_MASK);
         state->X[kk] =
            state->X[kk + MM] ^ (y >> 1) ^ param->mag01[y & 0x1UL];
      }
      for (; kk < NN - 1; kk++) {
         y = (state->X[kk] & UPPER_MASK) | (state->X[kk + 1] & LOWER_MASK);
         state->X[kk] = state->X[kk + (MM - NN)] ^ (y >> 1) ^
            param->mag01[y & 0x1UL];
      }
      y = (state->X[NN - 1] & UPPER_MASK) | (state->X[0] & LOWER_MASK);
      state->X[NN - 1] =
         state->X[MM - 1] ^ (y >> 1) ^ param->mag01[y & 0x1UL];

      state->s = 0;
   }

   y = state->X[state->s++];

   /* Tempering */
   y ^= (y >> 11);
   y ^= (y << 7) & 0x9d2c5680UL;
   y ^= (y << 15) & 0xefc60000UL;
   y ^= (y >> 18);

   return y;
}


/*-----------------------------------------------------------------------*/

/* generates a random number on (0,1)-real-interval
double genrand_real3(void) */
static double MT19937_02_U01 (void *vpar, void *vsta)
{
   return ((double) MT19937_02_Bits (vpar, vsta) + 0.5) *
           (1.0 / 4294967296.0);                /* divided by 2^32 */
}

/* These real versions are due to Isaku Wada, 2002/01/09 added */


/*-----------------------------------------------------------------------*/

static void WrMT19937 (void *);

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateMT19937_02 (unsigned long seed,
   unsigned long Key[], int len)
{
   unif01_Gen *gen;
   MT19937_param *param;
   MT19937_state *state;
   unsigned long S[NN];
   char name[LEN + 1];
   size_t leng;

   gen = CreateGFSR0 (NN, MM, 32, S, "");
   param = gen->param;
   state = gen->state;
   param->mag01[0] = 0x0UL;
   param->mag01[1] = MATRIX_A;
   gen->GetBits = &MT19937_02_Bits;
   gen->GetU01 = &MT19937_02_U01;
   gen->Write = &WrMT19937;
   strcpy (name, "ugfsr_CreateMT19937_02:");
   if ((len <= 0) || (NULL == Key)) {
      init_genrand (state, seed);
      addstr_Ulong (name, "   seed = ", seed);
   } else {
      state->s = NN + 1;
      init_by_array (state, Key, len);
      addstr_ArrayUlong (name, "   Key = ", len, Key);
   }
   leng = strlen (name);
   gen->name = util_Realloc (gen->name, (leng + 1) * sizeof (char));
   strncpy (gen->name, name, leng);
   gen->name[leng] = '\0';
   return gen;
}



/**************************************************************************/
/* A C-program for MT19937: Real number version -     */
/* generates one pseudorandom number with double precision    */
/* which is uniformly distributed on [0,1]-interval for each call.   */
/* CreateMT19937(seed) set initial values to the working area of 624 words.  */
/* (seed is any integer except for 0).       */

#undef NN
#undef MM
#define NN 624
#define MM 397

static double MT19937_98_U01 (void *vpar, void *vsta)
{
   MT19937_param *param = vpar;
   MT19937_state *state = vsta;
   unsigned long y;
   unsigned int j;

   if (state->s == NN) {                /* generate NN words at one time */
      for (j = 0; j < NN - MM; j++) {
         y = (state->X[j] & UPPER_MASK) | (state->X[j + 1] & LOWER_MASK);
         state->X[j] = state->X[j + MM] ^ (y >> 1) ^ param->mag01[y & 0x1];
      }
      for (; j < NN - 1; j++) {
         y = (state->X[j] & UPPER_MASK) | (state->X[j + 1] & LOWER_MASK);
         state->X[j] = state->X[j + (MM - NN)] ^ (y >> 1) ^
              param->mag01[y & 0x1];
      }
      y = (state->X[NN - 1] & UPPER_MASK) | (state->X[0] & LOWER_MASK);
      state->X[NN - 1] = state->X[MM - 1] ^ (y >> 1) ^ param->mag01[y & 0x1];

      state->s = 0;
   }
   y = state->X[state->s++];
   y ^= TEMPERING_SHIFT_U (y);
   y ^= TEMPERING_SHIFT_S (y) & TEMPERING_MASK_B;
   y ^= TEMPERING_SHIFT_T (y) & TEMPERING_MASK_C;
#ifndef IS_ULONG32
   y &= 0xffffffffUL;        /* you may delete this line if word size = 32 */
#endif
   y ^= TEMPERING_SHIFT_L (y);
   return (double) y / 0xffffffffUL;
}

/*-----------------------------------------------------------------------*/

static unsigned long MT19937_98_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (MT19937_98_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrMT19937 (void *vsta)
{
   GFSR_state *state = vsta;
   unsigned int j;
   if (unif01_WrLongStateFlag) {
      printf ("S = {\n ");
      for (j = 0; j < state->K; j++) {
         printf (" %12lu", state->X[j]);
         if (j < state->K - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateMT19937_98 (unsigned long seed)
/* 
 * Setting initial seeds to X[NN] using
 * the generator Line 25 of Table 1 in          
 * [KNUTH 1981, The Art of Computer Programming 
 * Vol. 2 (2nd Ed.), pp102]                    
 */
{
   unif01_Gen *gen;
   MT19937_param *param;
   unsigned int j;
   unsigned long S[NN];
   char name[LEN + 1];
   size_t leng;

   S[0] = seed & 0xffffffffUL;
   for (j = 1; j < NN; j++)
      S[j] = (69069 * S[j - 1]) & 0xffffffffUL;
   gen = CreateGFSR0 (NN, MM, 32, S, "");
   param = gen->param;
   param->mag01[0] = 0x0;
   param->mag01[1] = MATRIX_A;
   gen->GetBits = &MT19937_98_Bits;
   gen->GetU01  = &MT19937_98_U01;
   gen->Write   = &WrMT19937;
   strcpy (name, "ugfsr_CreateMT19937_98:");
   addstr_Ulong (name, "   seed = ", seed);
   leng = strlen (name);
   gen->name = util_Realloc (gen->name, (leng + 1) * sizeof (char));
   strncpy (gen->name, name, leng);
   gen->name[leng] = '\0';
   return gen;
}


/**************************************************************************/

static unsigned long GFSR5_Bits (void *vpar, void *vsta)
{
   GFSR5_param *param = vpar;
   GFSR5_state *state = vsta;

   if (++state->s == state->K)
      state->s = 0;
   if (++state->r1 == state->K)
      state->r1 = 0;
   if (++state->r2 == state->K)
      state->r2 = 0;
   if (++state->r3 == state->K)
      state->r3 = 0;
   state->X[state->s] ^= state->X[state->r1] ^ state->X[state->r2] ^
      state->X[state->r3];
   return state->X[state->s] << param->Shift;
}

/*-----------------------------------------------------------------------*/

static double GFSR5_U01 (void *vpar, void *vsta)
{
   return GFSR5_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateGFSR5 (unsigned int k, unsigned int r1,
   unsigned int r2, unsigned int r3, unsigned int l, unsigned long S[])
{
   unif01_Gen *gen;
   GFSR5_param *param;
   GFSR5_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int j;
   unsigned long Mask;

   util_Assert ((l >= 1) && (l <= 32),
      "ugfsr_CreateGFSR5:   l must be in [1..32]");
   util_Assert ((r3 > 0) && (r3 < r2),
      "ugfsr_CreateGFSR5:   we must have  0 < r3 < r2");
   util_Assert ((r1 > r2) && (r1 < k),
      "ugfsr_CreateGFSR5:   we must have  r2 < r1 < k");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (GFSR5_param));
   state = util_Malloc (sizeof (GFSR5_state));

   strcpy (name, "ugfsr_CreateGFSR5:");
   addstr_Uint (name, "   k = ", k);
   addstr_Uint (name, ",   r1 = ", r1);
   addstr_Uint (name, ",   r2 = ", r2);
   addstr_Uint (name, ",   r3 = ", r3);
   addstr_Uint (name, ",   l = ", l);
   addstr_ArrayUlong (name, ",   S = ", (int) k, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (l == 32)
      Mask = unif01_MASK32;
   else
      Mask = num_TwoExp[l] - 1.0;
   state->X = util_Calloc ((size_t) k, sizeof (unsigned long));
   for (j = 0; j < k; j++)
      state->X[j] = S[j] & Mask;

   state->r1 = k - r1;
   state->r2 = k - r2;
   state->r3 = k - r3;
   state->s = 0;
   state->K = k;
   param->Shift = 32 - l;

   gen->param = param;
   gen->state = state;
   gen->GetBits = &GFSR5_Bits;
   gen->GetU01  = &GFSR5_U01;
   gen->Write   = &WrGFSR;
   return gen;
}


/**************************************************************************/

#define A 471
#define B 1586
#define C 6988
#define D 9689
#define M 16383

#define RandomInteger (++state->s, state->X[state->s & M] = \
  state->X[(state->s-A) & M] ^ state->X[(state->s-B) & M] ^ \
  state->X[(state->s-C) & M] ^ state->X[(state->s-D) & M])

static void WrZiff98 (void *vsta)
{
   GFSR_state *state = vsta;
   unsigned int j;
   int s = state->s;

   s = (s - D) % (M + 1);
   if (unif01_WrLongStateFlag) {
      printf (" S = {\n ");
      for (j = 0; j < state->K; j++) {
         if (++s > M)
            s = 0;
         printf (" %12lu", state->X[s]);
         if (j < state->K - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

static unsigned long Ziff98_Bits (void *junk, void *vsta)
{
   GFSR_state *state = vsta;
   return RandomInteger;
}

/*-----------------------------------------------------------------------*/

static double Ziff98_U01 (void *junk, void *vsta)
{
   GFSR_state *state = vsta;
   return  RandomInteger * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

unif01_Gen *ugfsr_CreateZiff98 (unsigned long S[])
{
   unif01_Gen *gen;
   GFSR_state *state;
   size_t leng;
   char name[LEN + 1];
   unsigned int j, l = 32, Mask;

   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (GFSR_state));

   Mask = unif01_MASK32;
   state->X = util_Calloc ((size_t) M + 1, sizeof (unsigned long));
   for (j = 0; j < D; j++)
      state->X[j] = S[j] & Mask;
   state->s = D;
   state->K = D;

   strcpy (name, "ugfsr_CreateZiff98:");
   addstr_Uint (name, "   k = ", D);
   addstr_Uint (name, ",   r1 = ", C);
   addstr_Uint (name, ",   r2 = ", B);
   addstr_Uint (name, ",   r3 = ", A);
   addstr_Uint (name, ",   l = ", l);
   addstr_ArrayUlong (name, ",   S = ", D, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->param   = NULL;
   gen->state   = state;
   gen->GetBits = &Ziff98_Bits;
   gen->GetU01  = &Ziff98_U01;
   gen->Write   = &WrZiff98;
   return gen;
}

/**************************************************************************/

void ugfsr_DeleteGFSR5 (unif01_Gen *gen)
{
   GFSR5_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->X);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/**************************************************************************/

void ugfsr_DeleteGen (unif01_Gen *gen)
{
   GFSR_state *state;
   if (NULL == gen)
      return;
   state = gen->state;
   util_Free (state->X);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}
