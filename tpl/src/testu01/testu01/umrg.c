/*************************************************************************\
 *
 * Package:        TestU01
 * File:           umrg.c
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

#include "umrg.h"
#include "unif01.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


#ifdef USE_GMP
#include <gmp.h>
#endif



#define  LEN  300                 /* Max length of strings */
#define  MASK32 0xffffffffUL      /* 2^32 - 1 */





/*================================= Types =================================*/

typedef struct {
   int kind;
   long a1, q1, r1, a2, q2, r2;
   long M;
   double Norm;
} MRG2_param;

typedef struct {
   long x1, x2;
} MRG2_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   int kind;
   long a1, q1, r1, a3, q3, r3;
   long M;
   double Norm;
} MRG3_param;

typedef struct {
   long x1, x2, x3;
} MRG3_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   int kind;
   long a1, q1, r1, a5, q5, r5;
   long M;
   double Norm;
} MRG5_param;

typedef struct {
   long x1, x2, x3, x4, x5;
} MRG5_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   int kind;
   long a1, q1, r1, a7, q7, r7;
   long M;
   double Norm;
} MRG7_param;

typedef struct {
   long x1, x2, x3, x4, x5, x6, x7;
} MRG7_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   int kind;
   long *A, *Q, *R;
   long M;
   double Norm;
} MRG_param;

typedef struct {
   long *S;
   int k;
} MRG_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   double *A;
   double M;
   double Norm;
} MRGFloat_param;

typedef struct {
   double *S;
   int k;
} MRGFloat_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   long a12, a13, a21, a23;               /* Multipliers */
   long q12, q13, q21, q23;               /* qij = mi DIV aj */
   long r12, r13, r21, r23;               /* rij = mi MOD aj */
   long M1, M2;                           /* Modules */
   double Norm;
} CombMRG3_param;

typedef struct {
   long x10, x11, x12, x20, x21, x22;     /* State */
} CombMRG3_state;

/*-------------------------------------------------------------------------*/
#ifdef USE_GMP

typedef struct {
   mpz_t M, *A, W, T;
   mpf_t F, Norm;
   lebool *AnonZero;
} BigMRG_param;

typedef struct {
   mpz_t *S;
   int k;
} BigMRG_state;

typedef struct {
   mpz_t M1, M2, *A1, *A2, W, T1, T2;
   mpf_t F1, F2, Norm1, Norm2;
   lebool *A1nonZero, *A2nonZero;
} BigC2MRG_param;

typedef struct {
   mpz_t *S1, *S2;
   int k;
} BigC2MRG_state;

#endif
/*-------------------------------------------------------------------------*/

typedef struct {
   lebool Flag;         /* TRUE if k > r, FALSE if k < r */
   int Skip;
} LagFibFloat_param;

typedef struct {
   double *X;
   int r, s, RR;
   int Lag;
} LagFibFloat_state;

/*-------------------------------------------------------------------------*/

typedef struct {
   unsigned long Mask;   /* 2^t - 1 */
   int b;                /* shift |t - 32| bits, right or left */
   lebool LeftShift;    /* TRUE for left shift, FALSE for right shift */
   lebool Flag;         /* TRUE if k > r, FALSE if k < r */
   int Skip;
} LagFib_param;

typedef struct {
   unsigned long *X;
   int r, s, RR;
   int Lag;
} LagFib_state;

/*-------------------------------------------------------------------------*/






/*============================ functions ==================================*/


static void AddArrayString (char *name, const char *add, int high, char *A[])
/*
 * Add the string array A of size = high to name
 */
{
   int j;
   strcat (name, add);
   strcat (name, "(");
   strcat (name, A[0]);
   for (j = 1; (j < high) && (j < 7); j++) {
      strcat (name, ", ");
      if (A[j])
         strcat (name, A[j]);
      else
         strcat (name, "0");
   }
   if (high > 7)
      strcat (name, ", ... )");
   else
      strcat (name, ")");
}


/**************************************************************************/

static double MRG2_U01 (void *vpar, void *vsta)
/* 
 * Implementation used for k = 2,  a1 > 0,  a2 > 0
 */
{
   MRG2_param *param = vpar;
   MRG2_state *state = vsta;
   long h;
   long p;

   h = state->x2 / param->q2;
   p = param->a2 * (state->x2 - h * param->q2) - h * param->r2;
   if (p < 0)
      p += param->M;
   state->x2 = state->x1;

   h = state->x1 / param->q1;
   state->x1 = param->a1 * (state->x1 - h * param->q1) - h * param->r1;
   if (state->x1 <= 0)
      state->x1 += p;
   else
      state->x1 = (state->x1 - param->M) + p;
   if (state->x1 < 0)
      state->x1 += param->M;
   return state->x1 * param->Norm;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG2_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG2 (void *vsta)
{
   MRG2_state *state = vsta;
   printf (" S[1] = %10ld,   S[2] = %10ld\n\n", state->x1, state->x2);
}


/**************************************************************************/

static double MRG3_U01 (void *vpar, void *vsta)
/* 
 * Implementation used for k = 3,  a1 > 0,  a3 > 0 and a2 = 0
 */
{
   MRG3_param *param = vpar;
   MRG3_state *state = vsta;
   long h;
   long p;

   h = state->x3 / param->q3;
   p = param->a3 * (state->x3 - h * param->q3) - h * param->r3;
   if (p < 0)
      p += param->M;
   state->x3 = state->x2;
   state->x2 = state->x1;
   h = state->x1 / param->q1;
   state->x1 = param->a1 * (state->x1 - h * param->q1) - h * param->r1;
   if (state->x1 <= 0)
      state->x1 += p;
   else
      state->x1 = (state->x1 - param->M) + p;
   if (state->x1 < 0)
      state->x1 += param->M;
   return state->x1 * param->Norm;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG3_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG3 (void *vsta)
{
   MRG3_state *state = vsta;
   printf (" S[1] = %10ld,   S[2] = %10ld,   S[3] = %10ld\n\n",
      state->x1, state->x2, state->x3);
}


/**************************************************************************/

static double MRG5_U01 (void *vpar, void *vsta)
/* 
 * Implementation used for k = 5, a1 > 0, a5 > 0 and a2 = a3 = a4 = 0
 */
{
   MRG5_param *param = vpar;
   MRG5_state *state = vsta;
   long h;
   long p;

   h = state->x5 / param->q5;
   p = param->a5 * (state->x5 - h * param->q5) - h * param->r5;
   if (p < 0)
      p += param->M;
   state->x5 = state->x4;
   state->x4 = state->x3;
   state->x3 = state->x2;
   state->x2 = state->x1;
   h = state->x1 / param->q1;
   state->x1 = param->a1 * (state->x1 - h * param->q1) - h * param->r1;
   if (state->x1 <= 0)
      state->x1 += p;
   else
      state->x1 = (state->x1 - param->M) + p;
   if (state->x1 < 0)
      state->x1 += param->M;
   return state->x1 * param->Norm;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG5_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG5_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG5 (void *vsta)
{
   MRG5_state *state = vsta;
   printf (" S[1] = %10ld,   S[2] = %10ld,   S[3] = %10ld,\n",
      state->x1, state->x2, state->x3);
   printf (" S[4] = %10ld,   S[5] = %10ld\n\n", state->x4, state->x5);
}


/**************************************************************************/

static double MRG7_U01 (void *vpar, void *vsta)
/* 
 * Implementation used for k = 7, a1 > 0, a7 > 0, a2 = a3 = a4 = a5 = a6 = 0
 */
{
   MRG7_param *param = vpar;
   MRG7_state *state = vsta;
   long h;
   long p;

   h = state->x7 / param->q7;
   p = param->a7 * (state->x7 - h * param->q7) - h * param->r7;
   if (p < 0)
      p += param->M;
   state->x7 = state->x6;
   state->x6 = state->x5;
   state->x5 = state->x4;
   state->x4 = state->x3;
   state->x3 = state->x2;
   state->x2 = state->x1;
   h = state->x1 / param->q1;
   state->x1 = param->a1 * (state->x1 - h * param->q1) - h * param->r1;
   if (state->x1 <= 0)
      state->x1 += p;
   else
      state->x1 = (state->x1 - param->M) + p;
   if (state->x1 < 0)
      state->x1 += param->M;
   return state->x1 * param->Norm;
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG7_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG7_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG7 (void *vsta)
{
   MRG7_state *state = vsta;
   printf (" S[1] = %10ld,   S[2] = %10ld,   S[3] = %10ld,\n",
      state->x1, state->x2, state->x3);
   printf (" S[4] = %10ld,   S[5] = %10ld,   S[6] = %10ld,\n",
      state->x4, state->x5, state->x6);
   printf (" S[7] = %10ld\n\n", state->x7);
}


/**************************************************************************/

static double MRG_U01 (void *vpar, void *vsta)
/*
 * Generator MRG of order k when R[i] < Q[i] for all i. The values returned
 * by the generator MRG are in the set { 0, 1/m, 2/m, ..., (m-1)/m } 
 */
{
   MRG_param *param = vpar;
   MRG_state *state = vsta;
   long i, p, h, t;

   p = 0;
   for (i = state->k; i > 0; i--) {
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
   return (p * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRG_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRG (void *vsta)
{
   MRG_state *state = vsta;
   int i;
   if (unif01_WrLongStateFlag || (state->k < 8)) {
      printf (" S = {\n ");
      for (i = 1; i <= state->k; i++) {
	 printf ("   %12ld", state->S[i]);
         if (i < state->k)
            printf (",");
	 if ((i % 4) == 0)
	    printf ("\n ");
      }
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-------------------------------------------------------------------------*/

enum {MRG_ALL = 987654321};

static unif01_Gen * CreateMRG_all (long m, int k, long AA[], long SS[])
{
   unif01_Gen *gen;
   MRG_param *param;
   MRG_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;
   long *A, *S, *Q, *R;

   if ((k < 2) || (m <= 1))
      util_Error ("umrg_CreateMRG:   k < 2  or  m < 2");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (MRG_param));
   state = util_Malloc (sizeof (MRG_state));

   strncpy (name, "umrg_CreateMRG:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   k = ", k);
   addstr_ArrayLong (name, ",   A = ", k, AA);
   addstr_ArrayLong (name, ",   S = ", k, SS);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   A = util_Calloc ((size_t) k + 1, sizeof (long));
   R = util_Calloc ((size_t) k + 1, sizeof (long));
   Q = util_Calloc ((size_t) k + 1, sizeof (long));
   S = util_Calloc ((size_t) k + 1, sizeof (long));

   for (i = 1; i <= k; i++) {
      A[i] = AA[i - 1];
      S[i] = SS[i - 1];
      if (A[i] != 0) {
         R[i] = m % labs (A[i]);
         Q[i] = m / labs (A[i]);
      }
   }

   param->kind = MRG_ALL;
   param->M = m;
   param->Norm = 1.0 / m;
   param->A = A;
   param->R = R;
   param->Q = Q;
   state->k = k;
   state->S = S;

   gen->param = param;
   gen->state = state;
   gen->GetBits = &MRG_Bits;
   gen->GetU01 = &MRG_U01;
   gen->Write = &WrMRG;

   return gen;
}

/*-------------------------------------------------------------------------*/

static void DeleteMRG_all (unif01_Gen * gen)
{
   MRG_param *param;
   MRG_state *state;

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

unif01_Gen * umrg_CreateMRG (long m, int k, long A[], long S[])
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   int i, n;
   long r, q;
   long R[8] = { 0 };
   long Q[8] = { 0 };

   if ((k < 2) || (m <= 1))
      util_Error ("umrg_CreateMRG:   k < 2  or  m < 2");

   n = 0;
   for (i = 0; i < k; i++) {
      util_Assert (labs (A[i]) < m, "umrg_CreateMRG:   |A[i]| >= m");
      util_Assert (S[i] < m, "umrg_CreateMRG:   S[i] >= m");
      util_Assert (S[i] >= 0, "umrg_CreateMRG:   S[i] < 0");
      if (A[i] != 0) {
         r = m % labs (A[i]);
         q = m / labs (A[i]);
         util_Assert (r <= q, "umrg_CreateMRG:   r[i] > q[i]");
         if (i < 7) {
            R[i] = r;
            Q[i] = q;
         }
      }
      if (S[i] != 0)
         n++;
   }
   util_Assert (n > 0, "umrg_CreateMRG:   all S[i] are 0");


   if ((k == 2) && (A[0] > 0) && (A[1] > 0)) {
      MRG2_param *param;
      MRG2_state *state;

      gen = util_Malloc (sizeof (unif01_Gen));
      param = util_Malloc (sizeof (MRG2_param));
      state = util_Malloc (sizeof (MRG2_state));
      param->kind = 2;
      param->M = m;
      param->Norm = 1.0 / m;
      param->a1 = A[0];
      param->r1 = R[0];
      param->q1 = Q[0];
      param->a2 = A[1];
      param->r2 = R[1];
      param->q2 = Q[1];
      state->x1 = S[0];
      state->x2 = S[1];

      gen->param = param;
      gen->state = state;
      gen->GetBits = &MRG2_Bits;
      gen->GetU01 = &MRG2_U01;
      gen->Write = &WrMRG2;

   } else if ((k == 3) && (A[0] > 0) && (A[1] == 0) && (A[2] > 0)) {
      MRG3_param *param;
      MRG3_state *state;

      gen = util_Malloc (sizeof (unif01_Gen));
      param = util_Malloc (sizeof (MRG3_param));
      state = util_Malloc (sizeof (MRG3_state));
      param->kind = 3;
      param->M = m;
      param->Norm = 1.0 / m;
      param->a1 = A[0];
      param->r1 = R[0];
      param->q1 = Q[0];
      param->a3 = A[2];
      param->r3 = R[2];
      param->q3 = Q[2];
      state->x1 = S[0];
      state->x2 = S[1];
      state->x3 = S[2];

      gen->param = param;
      gen->state = state;
      gen->GetBits = &MRG3_Bits;
      gen->GetU01 = &MRG3_U01;
      gen->Write = &WrMRG3;

   } else if ((k == 5) && (A[0] > 0) && (A[1] == 0) && (A[2] == 0) &&
	      (A[3] == 0) && (A[4] > 0)) {
      MRG5_param *param;
      MRG5_state *state;

      gen = util_Malloc (sizeof (unif01_Gen));
      param = util_Malloc (sizeof (MRG5_param));
      state = util_Malloc (sizeof (MRG5_state));
      param->kind = 5;
      param->M = m;
      param->Norm = 1.0 / m;
      param->a1 = A[0];
      param->r1 = R[0];
      param->q1 = Q[0];
      param->a5 = A[4];
      param->r5 = R[4];
      param->q5 = Q[4];
      state->x1 = S[0];
      state->x2 = S[1];
      state->x3 = S[2];
      state->x4 = S[3];
      state->x5 = S[4];

      gen->param = param;
      gen->state = state;
      gen->GetBits = &MRG5_Bits;
      gen->GetU01 = &MRG5_U01;
      gen->Write = &WrMRG5;

   } else if ((k == 7) && (A[0] > 0) && (A[1] == 0) && (A[2] == 0) &&
	      (A[3] == 0) && (A[4] == 0) && (A[5] == 0) && (A[6] > 0)) {
      MRG7_param *param;
      MRG7_state *state;

      gen = util_Malloc (sizeof (unif01_Gen));
      param = util_Malloc (sizeof (MRG7_param));
      state = util_Malloc (sizeof (MRG7_state));
      param->kind = 7;
      param->M = m;
      param->Norm = 1.0 / m;
      param->a1 = A[0];
      param->r1 = R[0];
      param->q1 = Q[0];
      param->a7 = A[6];
      param->r7 = R[6];
      param->q7 = Q[6];
      state->x1 = S[0];
      state->x2 = S[1];
      state->x3 = S[2];
      state->x4 = S[3];
      state->x5 = S[4];
      state->x6 = S[5];
      state->x7 = S[6];

      gen->param = param;
      gen->state = state;
      gen->GetBits = &MRG7_Bits;
      gen->GetU01 = &MRG7_U01;
      gen->Write = &WrMRG7;

   } else {
      return CreateMRG_all (m, k, A, S);
   }

   strncpy (name, "umrg_CreateMRG:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   k = ", k);
   addstr_ArrayLong (name, ",   A = ", k, A);
   addstr_ArrayLong (name, ",   S = ", k, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   return gen;
}

/*-------------------------------------------------------------------------*/

void umrg_DeleteMRG (unif01_Gen * gen)
{
   MRG_param *param;
   if (NULL == gen)
      return;
   param = (MRG_param *) gen->param;
   if (param->kind == MRG_ALL)
      DeleteMRG_all (gen);
   else
      unif01_DeleteGen (gen);
}


/**************************************************************************/

static double MRGFloat_U01 (void *vpar, void *vsta)
/*
 * Generator MRG of order k. Implementation uses floating-point.
 */
{
   MRGFloat_param *param = vpar;
   MRGFloat_state *state = vsta;
   double p;
   long j;
   int i;

   p = 0.0;
   for (i = state->k; i > 0; i--) {
      if (param->A[i] != 0.0)
         p += param->A[i] * state->S[i];
      if (i > 1)
         state->S[i] = state->S[i - 1];
   }
   j = p / param->M;
   if (p >= 0.0)
      p -= j * param->M;          /* p = sum ...  mod m */
   else {
      p += (1 - j) * param->M;
      /* When truncating, if p % M = 0, we will get p = M */
      if (p >= param->M)
         p -= param->M;
   }
   state->S[1] = p;               /* Must be in {0,..., m-1} */
   return (p * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long MRGFloat_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * MRGFloat_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrMRGFloat (void *vsta)
{
   MRGFloat_state *state = vsta;
   int i;
   if (unif01_WrLongStateFlag || (state->k < 8)) {
      printf (" S = {\n ");
      for (i = 1; i <= state->k; i++) {
	 printf ("   %12.0f", state->S[i]);
         if (i < state->k)
            printf (",");
	 if ((i % 4) == 0)
	    printf ("\n ");
      }
      printf ("    }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-------------------------------------------------------------------------*/

unif01_Gen * umrg_CreateMRGFloat (long m, int k, long AA[], long SS[])
{
   unif01_Gen *gen;
   MRGFloat_param *param;
   MRGFloat_state *state;
   size_t leng;
   char name[LEN + 1];
   int i, n;
   double *A, *S;
   double pr1, pr2;

   if ((k < 2) || (m < 2))
      util_Error ("umrg_CreateMRGFloat:    k < 2  or  m < 2");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (MRGFloat_param));
   state = util_Malloc (sizeof (MRGFloat_state));

   strncpy (name, "umrg_CreateMRGFloat:", (size_t) LEN);
   addstr_Long (name, "   m = ", m);
   addstr_Long (name, ",   k = ", k);
   addstr_ArrayLong (name, ",   A = ", k, AA);
   addstr_ArrayLong (name, ",   S = ", k, SS);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   A = util_Calloc ((size_t) k + 1, sizeof (double));
   S = util_Calloc ((size_t) k + 1, sizeof (double));

   n = 0;
   pr1 = 0.0;
   pr2 = 0.0;
   for (i = 1; i <= k; i++) {
      if ((AA[i - 1] >= m) || (-AA[i - 1] >= m))
         util_Error ("umrg_CreateMRGFloat:   |A[i]| >= m");
      if ((SS[i - 1] >= m) || (-SS[i - 1] >= m))
         util_Error ("umrg_CreateMRGFloat:   |S[i]| >= m");
      A[i] = AA[i - 1];
      S[i] = SS[i - 1];
      if (SS[i - 1] < 0)
         S[i] += m;
      if (AA[i - 1] < 0)
         pr2 -= A[i];
      else
         pr1 += A[i];
      if (SS[i - 1] != 0)
         n++;
   }
   if (n == 0)
      util_Error (" umrg_CreateMRGFloat:   all S[i] are 0");
   if ((pr1 * m >= num_TwoExp[53]) || (pr2 * m >= num_TwoExp[53]))
      util_Error ("umrg_CreateMRGFloat:   Condition on a_i not valid");

   param->M = m;
   param->Norm = 1.0 / m;
   param->A = A;
   state->k = k;
   state->S = S;

   gen->param = param;
   gen->state = state;
   gen->GetBits = &MRGFloat_Bits;
   gen->GetU01 = &MRGFloat_U01;
   gen->Write = &WrMRGFloat;

   return gen;
}

/*-------------------------------------------------------------------------*/

void umrg_DeleteMRGFloat (unif01_Gen * gen)
{
   MRGFloat_param *param;
   MRGFloat_state *state;

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


/**************************************************************************/

static double CombMRG3_U01 (void *vpar, void *vsta)
/*
 * Implementation used for k = 3. It assumes that
 *    a11 = 0, a12 > 0, a13 < 0 and a21 > 0, a22 = 0, a23 < 0.
 */
{
   CombMRG3_param *param = vpar;
   CombMRG3_state *state = vsta;
   long h, p12, p13, p21, p23, Z;

   /* Component 1 */
   h = state->x10 / param->q13;
   p13 = param->a13 * (state->x10 - h * param->q13) - h * param->r13;
   if (p13 < 0)
      p13 += param->M1;
   util_Assert (p13 >= 0,
      "umrg_CreateC2MRG:   invalid parameters for a_{1,3}");
   h = state->x11 / param->q12;
   p12 = param->a12 * (state->x11 - h * param->q12) - h * param->r12;
   if (p12 < 0)
      p12 += param->M1;
   util_Assert (p12 >= 0,
      "umrg_CreateC2MRG:   invalid parameters for a_{1,2}");
   state->x10 = state->x11;
   state->x11 = state->x12;
   state->x12 = p12 - p13;
   if (state->x12 < 0)
      state->x12 += param->M1;

   /* Component 2 */
   h = state->x20 / param->q23;
   p23 = param->a23 * (state->x20 - h * param->q23) - h * param->r23;
   if (p23 < 0)
      p23 += param->M2;
   util_Assert (p23 >= 0,
      "umrg_CreateC2MRG:   invalid parameters for a_{2,3}");
   h = state->x22 / param->q21;
   p21 = param->a21 * (state->x22 - h * param->q21) - h * param->r21;
   if (p21 < 0)
      p21 += param->M2;
   util_Assert (p21 >= 0,
      "umrg_CreateC2MRG:   invalid parameters for a_{2,1}");
   state->x20 = state->x21;
   state->x21 = state->x22;
   state->x22 = p21 - p23;
   if (state->x22 < 0)
      state->x22 += param->M2;

   /* Combinaison */
   if (state->x12 < state->x22)
      Z = state->x12 - state->x22 + param->M1;
   else
      Z = state->x12 - state->x22;
   return (Z * param->Norm);
}

/*-------------------------------------------------------------------------*/

static unsigned long CombMRG3_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * CombMRG3_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrCombMRG3 (void *vsta)
{
   CombMRG3_state *state = vsta;
   printf ("   S1[0] = %1ld   S1[1] = %1ld   S1[2] = %1ld\n",
      state->x10, state->x11, state->x12);
   printf ("   S2[0] = %1ld   S2[1] = %1ld   S2[2] = %1ld\n\n",
      state->x20, state->x21, state->x22);
}

/*-------------------------------------------------------------------------*/

unif01_Gen * umrg_CreateC2MRG (long m1, long m2, int k, long A1[], long A2[],
   long S1[], long S2[])
{
   unif01_Gen *gen;
   CombMRG3_param *param;
   CombMRG3_state *state;
   size_t leng;
   char name[LEN + 1];
   int i;
   long cs1[4], cs2[4], cq1[4] = {0}, cq2[4] = {0}, cr1[4] = {0},
        cr2[4] = {0}, ca1[4], ca2[4];

   if (k != 3)
      util_Error ("umrg_CreateC2MRG:   k != 3");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (CombMRG3_param));
   state = util_Malloc (sizeof (CombMRG3_state));

   strncpy (name, "umrg_CreateC2MRG:", (size_t) LEN);
   addstr_Long (name, "   m1 = ", m1);
   addstr_Long (name, ",   m2 = ", m2);
   addstr_Long (name, ",   k = ", k);
   addstr_ArrayLong (name, ",   A1 = ", k, A1);
   addstr_ArrayLong (name, ",   S1 = ", k, S1);
   addstr_ArrayLong (name, ",   A2 = ", k, A2);
   addstr_ArrayLong (name, ",   S2 = ", k, S2);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   for (i = 1; i <= k; i++) {
      ca1[i] = A1[i - 1];
      ca2[i] = A2[i - 1];
      cs1[i] = S1[i - 1];
      cs2[i] = S2[i - 1];
      if (ca1[i] != 0) {
         cr1[i] = m1 % labs (ca1[i]);
         cq1[i] = m1 / labs (ca1[i]);
         util_Assert ((double) cr1[i] * labs(ca1[i]) < m1,
            "umrg_CreateC2MRG:   |A1[i]| * (m1 mod |A1[i]|) >= m1");
      }
      if (ca2[i] != 0) {
         cr2[i] = m2 % labs (ca2[i]);
         cq2[i] = m2 / labs (ca2[i]);
         util_Assert ((double) cr2[i] * labs(ca2[i]) < m2,
            "umrg_CreateC2MRG:pp   |A2[i]| * (m2 mod |A2[i]|) >= m2");
      }
   }
   param->M1 = m1;
   param->M2 = m2;
   param->Norm = 1.0 / m1;
   if (k == 3) {
      param->a12 = ca1[2];
      param->a13 = ca1[3];
      param->a21 = ca2[1];
      param->a23 = ca2[3];
      param->q12 = cq1[2];
      param->q13 = cq1[3];
      param->q21 = cq2[1];
      param->q23 = cq2[3];
      param->r12 = cr1[2];
      param->r13 = cr1[3];
      param->r21 = cr2[1];
      param->r23 = cr2[3];
      state->x10 = cs1[1];
      state->x11 = cs1[2];
      state->x12 = cs1[3];
      state->x20 = cs2[1];
      state->x21 = cs2[2];
      state->x22 = cs2[3];

      gen->GetBits = &CombMRG3_Bits;
      gen->GetU01 = &CombMRG3_U01;
      gen->Write = &WrCombMRG3;
   }
   gen->param = param;
   gen->state = state;
   return gen;
}

/*-------------------------------------------------------------------------*/

void umrg_DeleteC2MRG (unif01_Gen *gen)
{
   unif01_DeleteGen (gen);
}


/***************************************************************************/

#ifdef USE_GMP

static void WrBigMRG (void *vsta)
{
   BigMRG_state *state = vsta;
   int i;

   printf (" S = {\n   ");
   for (i = 1; i <= state->k; i++) {
      mpz_out_str (NULL, 10, state->S[i]);
      if (i < state->k)
         printf (",\n   ");
   }
   printf ("\n   }\n");
}

/*-------------------------------------------------------------------------*/

static double BigMRG_U01 (void *vpar, void *vsta)
{
   BigMRG_param *param = vpar;
   BigMRG_state *state = vsta;
   int i;

   mpz_set_ui (param->T, 0);     /* T = 0 */

   for (i = state->k; i >= 1; i--) {
      if (param->AnonZero[i]) {
         mpz_mul (param->W, param->A[i], state->S[i]);  /* W = A[i]*S[i] */
         mpz_mod (param->W, param->W, param->M);        /* W = W % M */
         mpz_add (param->T, param->W, param->T);        /* T = W + T */
      }
      if (i > 1) {
	 mpz_set (state->S[i], state->S[i - 1]);      /* S[i] = S[i - 1] */
      } else {
         mpz_mod (state->S[i], param->T, param->M);   /* S[i] = T % M */
      }
   }
   mpf_set_z (param->F, state->S[1]);                 /* F = S[1] */
   mpf_mul (param->F, param->F, param->Norm);         /* F = F * Norm */
   return  mpf_get_d (param->F);                      /* U  = F */
}

/*-------------------------------------------------------------------------*/

static unsigned long BigMRG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * BigMRG_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

unif01_Gen *umrg_CreateBigMRG (char *m, int k, char *A[], char *S[])
{
   unif01_Gen *gen;
   BigMRG_param *param;
   BigMRG_state *state;
   size_t len1, len2;
   char name[LEN + 1];
   int i, flag;

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (BigMRG_param));
   state = util_Malloc (sizeof (BigMRG_state));

   /* These flags are set to 0 if the corresponding element is 0 */
   param->AnonZero = util_Calloc ((size_t) k + 1, sizeof (lebool));

   strncpy (name, "umrg_CreateBigMRG:", (size_t) LEN);
   strcat (name, "   m = ");
   strcat (name, m);
   addstr_Long (name, ",   k = ", k);
   len1 = strlen (name);
   len2 = len1 + 4 * strlen (",    A =  ");

   for (i = 0; i < k; i++) {
      if (A[i]) {
         len2 += strlen (A[i]);
         param->AnonZero[i + 1] = TRUE;
      } else
         /* If NULL pointer, set corresponding element to 0 */
         param->AnonZero[i + 1] = FALSE;

      if (S[i])
         len2 += strlen (S[i]);
   }

   gen->name = util_Calloc (len2 + LEN, sizeof (char));
   strncpy (gen->name, name, len2);
   AddArrayString (gen->name, ",   A = ", k, A);
   AddArrayString (gen->name, ",   S = ", k, S);
   len1 = strlen (gen->name);
   gen->name = util_Realloc (gen->name, (len1 + 1) * sizeof (char));

   state->k = k;
   mpz_init (param->M);
   mpz_init (param->T);
   mpz_init (param->W);
   mpz_set_str (param->M, m, 0);

   param->A = util_Calloc ((size_t) k + 1, sizeof (mpz_t));
   state->S = util_Calloc ((size_t) k + 1, sizeof (mpz_t));

   for (i = 1; i <= k; i++) {
      if (param->AnonZero[i]) {
         mpz_init (param->A[i]);
         mpz_set_str (param->A[i], A[i - 1], 0);
         mpz_abs (param->T, param->A[i]);
         flag = mpz_cmp (param->T, param->M);
         util_Assert (flag < 0, "umrg_CreateBigMRG:   one A >= m");
         /* Check if element is 0 */
         if (mpz_sgn (param->A[i]) == 0) {
            param->AnonZero[i] = FALSE;
            mpz_clear (param->A[i]);
         }
      }
      mpz_init (state->S[i]);
      if (S[i - 1]) {
         mpz_set_str (state->S[i], S[i - 1], 0);
         mpz_abs (param->T, state->S[i]);
         flag = mpz_cmp (param->T, param->M);
         util_Assert (flag < 0, "umrg_CreateBigMRG:   one S >= m");
      }
   }

   mpf_set_default_prec (DBL_MANT_DIG);     /* Set precision to 53 bits */
   mpf_init (param->F);
   mpf_init (param->Norm);
   mpf_set_z (param->Norm, param->M);
   mpf_ui_div (param->Norm, 1, param->Norm);    /* Norm = 1 / m */

   gen->param = param;
   gen->state = state;
   gen->GetBits = &BigMRG_Bits;
   gen->GetU01 = &BigMRG_U01;
   gen->Write = &WrBigMRG;

   return gen;
}

/*-------------------------------------------------------------------------*/

void umrg_DeleteBigMRG (unif01_Gen * gen)
{
   BigMRG_param *param;
   BigMRG_state *state;
   int i;

   if (NULL == gen)
      return;
   param = gen->param;
   state = gen->state;
   mpf_clear (param->F);
   mpf_clear (param->Norm);

   for (i = 1; i <= state->k; i++) {
      if (param->AnonZero[i])
         mpz_clear (param->A[i]);
      mpz_clear (state->S[i]);
   }

   mpz_clear (param->M);
   mpz_clear (param->T);
   mpz_clear (param->W);

   util_Free (param->AnonZero);
   util_Free (state->S);
   util_Free (param->A);

   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/***************************************************************************/

static void WrBigC2MRG (void *vsta)
{
   BigC2MRG_state *state = vsta;
   int i;

   printf (" S1 = {\n   ");
   for (i = 1; i <= state->k; i++) {
      mpz_out_str (NULL, 10, state->S1[i]);
      if (i < state->k)
         printf (",\n   ");
   }
   printf ("\n   }\n\n S2 = {\n   ");
   for (i = 1; i <= state->k; i++) {
      mpz_out_str (NULL, 10, state->S2[i]);
      if (i < state->k)
         printf (",\n   ");
      else
         printf ("\n   }\n");
   }
}

/*-------------------------------------------------------------------------*/

static double BigC2MRG_U01 (void *vpar, void *vsta)
{
   BigC2MRG_param *param = vpar;
   BigC2MRG_state *state = vsta;
   int i;
   double U, U2;

   mpz_set_ui (param->T1, 0);     /* T1 = 0 */
   mpz_set_ui (param->T2, 0);     /* T2 = 0 */

   for (i = state->k; i >= 1; i--) {
      if (param->A1nonZero[i]) {
         mpz_mul (param->W, param->A1[i], state->S1[i]); /* W = A1[i]*S1[i] */
         mpz_mod (param->W, param->W, param->M1);        /* W = W % M1 */
         mpz_add (param->T1, param->W, param->T1);       /* T1 = W + T1 */
      }
      if (param->A2nonZero[i]) {
         mpz_mul (param->W, param->A2[i], state->S2[i]); /* W = A2[i]*S2[i] */
         mpz_mod (param->W, param->W, param->M2);        /* W = W % M2 */
         mpz_add (param->T2, param->W, param->T2);       /* T2 = W + T2 */
      }
      if (i > 1) {
	 mpz_set (state->S1[i], state->S1[i - 1]);     /* S1[i] = S1[i - 1] */
         mpz_set (state->S2[i], state->S2[i - 1]);     /* S2[i] = S2[i - 1] */
      } else {
         mpz_mod (state->S1[i], param->T1, param->M1); /* S1[i] = T1 % M1 */
         mpz_mod (state->S2[i], param->T2, param->M2); /* S2[i] = T2 % M2 */
      }
   }
   mpf_set_z (param->F1, state->S1[1]);                /* F1 = S1[1] */
   mpf_set_z (param->F2, state->S2[1]);                /* F2 = S2[1] */
   mpf_mul (param->F1, param->F1, param->Norm1);       /* F1 = F1 * Norm1 */
   mpf_mul (param->F2, param->F2, param->Norm2);       /* F2 = F2 * Norm2 */
   U = mpf_get_d (param->F1);                          /* U  = F1 */
   U2 = mpf_get_d (param->F2);                         /* U2 = F2 */
   
   U = U - U2;
   if (U < 0.0)
      return U + 1.0;
   if (U < 1.0)
      return U;
   return U - 1.0;
}

/*-------------------------------------------------------------------------*/

static unsigned long BigC2MRG_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * BigC2MRG_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

unif01_Gen *umrg_CreateBigC2MRG (char *m1, char *m2, int k,
   char *A1[], char *A2[], char *S1[], char *S2[])
{
   unif01_Gen *gen;
   BigC2MRG_param *param;
   BigC2MRG_state *state;
   size_t len1, len2;
   char name[LEN + 1];
   int i, flag;

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (BigC2MRG_param));
   state = util_Malloc (sizeof (BigC2MRG_state));

   /* These flags are set to 0 if the corresponding element is 0 */
   param->A1nonZero = util_Calloc ((size_t) k + 1, sizeof (lebool));
   param->A2nonZero = util_Calloc ((size_t) k + 1, sizeof (lebool));

   strncpy (name, "umrg_CreateBigC2MRG:", (size_t) LEN);
   strcat (name, "   m1 = ");
   strcat (name, m1);
   strcat (name, ",   m2 = ");
   strcat (name, m2);
   addstr_Long (name, ",   k = ", k);
   len1 = strlen (name);
   len2 = len1 + 4 * strlen (",    A1 =  ");

   for (i = 0; i < k; i++) {
      if (A1[i]) {
         len2 += strlen (A1[i]);
         param->A1nonZero[i + 1] = TRUE;
      } else
         /* If NULL pointer, set corresponding element to 0 */
         param->A1nonZero[i + 1] = FALSE;

      if (A2[i]) {
         len2 += strlen (A2[i]);
         param->A2nonZero[i + 1] = TRUE;
      } else
         param->A2nonZero[i + 1] = FALSE;

      if (S1[i])
         len2 += strlen (S1[i]);

      if (S2[i])
         len2 += strlen (S1[i]);
   }

   gen->name = util_Calloc (len2 + LEN, sizeof (char));
   strncpy (gen->name, name, len2);
   AddArrayString (gen->name, ",   A1 = ", k, A1);
   AddArrayString (gen->name, ",   A2 = ", k, A2);
   AddArrayString (gen->name, ",   S1 = ", k, S1);
   AddArrayString (gen->name, ",   S2 = ", k, S2);
   len1 = strlen (gen->name);
   gen->name = util_Realloc (gen->name, (len1 + 1) * sizeof (char));

   state->k = k;
   mpz_init (param->M1);
   mpz_init (param->M2);
   mpz_init (param->T1);
   mpz_init (param->T2);
   mpz_init (param->W);
   mpz_set_str (param->M1, m1, 0);
   mpz_set_str (param->M2, m2, 0);

   param->A1 = util_Calloc ((size_t) k + 1, sizeof (mpz_t));
   param->A2 = util_Calloc ((size_t) k + 1, sizeof (mpz_t));
   state->S1 = util_Calloc ((size_t) k + 1, sizeof (mpz_t));
   state->S2 = util_Calloc ((size_t) k + 1, sizeof (mpz_t));

   for (i = 1; i <= k; i++) {
      if (param->A1nonZero[i]) {
         mpz_init (param->A1[i]);
         mpz_set_str (param->A1[i], A1[i - 1], 0);
         mpz_abs (param->T1, param->A1[i]);
         flag = mpz_cmp (param->T1, param->M1);
         util_Assert (flag < 0, "umrg_CreateBigC2MRG:   one A1 >= m1");
         /* Check if element is 0 */
         if (mpz_sgn (param->A1[i]) == 0) {
            param->A1nonZero[i] = FALSE;
            mpz_clear (param->A1[i]);
         }
      }
      if (param->A2nonZero[i]) {
         mpz_init (param->A2[i]);
         mpz_set_str (param->A2[i], A2[i - 1], 0);
         mpz_abs (param->T1, param->A2[i]);
         flag = mpz_cmp (param->T1, param->M2);
         util_Assert (flag < 0, "umrg_CreateBigC2MRG:   one A2 >= m2");
         if (mpz_sgn (param->A2[i]) == 0) {
            param->A2nonZero[i] = FALSE;
            mpz_clear (param->A2[i]);
         }
      }
      mpz_init (state->S1[i]);
      if (S1[i - 1]) {
         mpz_set_str (state->S1[i], S1[i - 1], 0);
         mpz_abs (param->T1, state->S1[i]);
         flag = mpz_cmp (param->T1, param->M1);
         util_Assert (flag < 0, "umrg_CreateBigC2MRG:   one S1 >= m1");
      }

      mpz_init (state->S2[i]);
      if (S2[i - 1]) {
         mpz_set_str (state->S2[i], S2[i - 1], 0);
         mpz_abs (param->T1, state->S2[i]);
         flag = mpz_cmp (param->T1, param->M2);
         util_Assert (flag < 0, "umrg_CreateBigC2MRG:   one S2 >= m2");
      }
   }

   mpf_set_default_prec (DBL_MANT_DIG);     /* Set precision to 53 bits */
   mpf_init (param->F1);
   mpf_init (param->F2);
   mpf_init (param->Norm1);
   mpf_init (param->Norm2);
   mpf_set_z (param->Norm1, param->M1);
   mpf_set_z (param->Norm2, param->M2);
   mpf_ui_div (param->Norm1, 1, param->Norm1);    /* Norm1 = 1 / m1 */
   mpf_ui_div (param->Norm2, 1, param->Norm2);    /* Norm2 = 1 / m2 */

   gen->param = param;
   gen->state = state;
   gen->GetBits = &BigC2MRG_Bits;
   gen->GetU01 = &BigC2MRG_U01;
   gen->Write = &WrBigC2MRG;

   return gen;
}

/*-------------------------------------------------------------------------*/

void umrg_DeleteBigC2MRG (unif01_Gen * gen)
{
   BigC2MRG_param *param;
   BigC2MRG_state *state;
   int i;

   if (NULL == gen)
      return;
   param = gen->param;
   state = gen->state;
   mpf_clear (param->F1);
   mpf_clear (param->F2);
   mpf_clear (param->Norm1);
   mpf_clear (param->Norm2);

   for (i = 1; i <= state->k; i++) {
      if (param->A1nonZero[i])
         mpz_clear (param->A1[i]);
      if (param->A2nonZero[i])
         mpz_clear (param->A2[i]);
      mpz_clear (state->S1[i]);
      mpz_clear (state->S2[i]);
   }

   mpz_clear (param->M1);
   mpz_clear (param->M2);
   mpz_clear (param->T1);
   mpz_clear (param->T2);
   mpz_clear (param->W);

   util_Free (param->A1nonZero);
   util_Free (param->A2nonZero);

   util_Free (state->S1);
   util_Free (state->S2);
   util_Free (param->A1);
   util_Free (param->A2);

   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}

#endif

/***************************************************************************/

static double LagFibAddFloat_U01 (void *junk, void *vsta)
{
   LagFibFloat_state *state = vsta;
   double temp;

   temp = state->X[state->r] + state->X[state->s];
   if (temp >= 1.0)
      temp -= 1.0;
   state->X[state->r] = temp;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static double LagFibSousFloat_U01 (void *vpar, void *vsta)
{
   LagFibFloat_param *param = vpar;
   LagFibFloat_state *state = vsta;
   double temp;

   if (param->Flag)
      temp = state->X[state->r] - state->X[state->s];
   else
      temp = state->X[state->s] - state->X[state->r];
   if (temp < 0.0)
      temp += 1.0;
   state->X[state->r] = temp;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static double LagFibAddFloatLux_U01 (void *vpar, void *vsta)
{
   LagFibFloat_param *param = vpar;
   LagFibFloat_state *state = vsta;
   double temp;
   int j;

   if (--state->RR == 0) {
      state->RR = state->Lag;
      for (j = 0; j < param->Skip; j++) {
         temp = state->X[state->r] + state->X[state->s];
         state->X[state->r] = (temp >= 1.0) ? temp - 1.0 : temp;
         if (--state->r == 0)
            state->r = state->Lag;
         if (--state->s == 0)
            state->s = state->Lag;
      }
   }
   temp = state->X[state->r] + state->X[state->s];
   if (temp >= 1.0)
      temp -= 1.0;
   state->X[state->r] = temp;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static double LagFibSousFloatLux_U01 (void *vpar, void *vsta)
{
   LagFibFloat_param *param = vpar;
   LagFibFloat_state *state = vsta;
   double temp;
   int j;

   if (--state->RR == 0) {
      state->RR = state->Lag;
      for (j = 0; j < param->Skip; j++) {
         if (param->Flag)
            temp = state->X[state->r] - state->X[state->s];
         else
            temp = state->X[state->s] - state->X[state->r];
            state->X[state->r] = (temp < 0.0) ? temp + 1.0 : temp;
         if (--state->r == 0)
            state->r = state->Lag;
         if (--state->s == 0)
            state->s = state->Lag;
      }
   }
   if (param->Flag)
      temp = state->X[state->r] - state->X[state->s];
   else
      temp = state->X[state->s] - state->X[state->r];
   if (temp < 0.0)
      temp += 1.0;
   state->X[state->r] = temp;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibAddFloat_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LagFibAddFloat_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibSousFloat_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LagFibSousFloat_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibAddFloatLux_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LagFibAddFloatLux_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibSousFloatLux_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (unif01_NORM32 * LagFibSousFloatLux_U01 (vpar, vsta));
}

/*-------------------------------------------------------------------------*/

static void WrLagFibFloat (void *vsta)
{
   LagFibFloat_state *state = vsta;
   int j;
   unsigned long Z;

   if (unif01_WrLongStateFlag) {
      printf ("S = {\n");
      for (j = 0; j < state->Lag; j++) {
         Z = unif01_NORM32 * state->X[state->r];
         printf (" %12lu", Z);
         if (--state->r == 0)
            state->r = state->Lag;
         if (j < state->Lag - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n");
      }
      printf ("   }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-------------------------------------------------------------------------*/

unif01_Gen * umrg_CreateLagFibFloat (int k, int r, char Op, int Lux,
   unsigned long S[])
{
   unif01_Gen *gen;
   LagFibFloat_param *param;
   LagFibFloat_state *state;
   size_t leng;
   char name[LEN + 1];
   char chaine[2];
   int j;
   double *X;

/* util_Assert (k > r, "umrg_CreateLagFibFloat:   k <= r"); */
   util_Assert ((Op == '-') || (Op == '+'),
      "umrg_CreateLagFibFloat:  only + and - are implemented");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LagFibFloat_param));
   state = util_Malloc (sizeof (LagFibFloat_state));

   strncpy (name, "umrg_CreateLagFibFloat:", (size_t) LEN);
   addstr_Long (name, "   k = ", k);
   addstr_Long (name, ",   r = ", r);
   strcat (name, ",   Op = ");
   sprintf (chaine, "%c", Op);
   strcat (name, chaine);
   addstr_Long (name, ",   Lux = ", Lux);
   if (k >= r)
      addstr_ArrayUlong (name, ",   S = ", k, S);
   else
      addstr_ArrayUlong (name, ",   S = ", r, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   if (k >= r) {
      state->Lag = k;
      state->r = k;
      state->s = r;
      param->Flag = TRUE;
   } else {
      state->Lag = r;
      state->r = r;
      state->s = k;
      param->Flag = FALSE;
   }
   param->Skip = Lux - state->Lag;

   if (param->Skip > 0) {
      X = util_Calloc ((size_t) Lux + 1, sizeof (double));
      state->RR = state->Lag;
      switch (Op) {
      case '+':
         gen->GetBits = &LagFibAddFloatLux_Bits;
         gen->GetU01 = &LagFibAddFloatLux_U01;
         break;
      case '-':
         gen->GetBits = &LagFibSousFloatLux_Bits;
         gen->GetU01 = &LagFibSousFloatLux_U01;
         break;
      }
   } else {
      X = util_Calloc ((size_t) state->Lag + 1, sizeof (double));
      switch (Op) {
      case '+':
         gen->GetBits = &LagFibAddFloat_Bits;
         gen->GetU01 = &LagFibAddFloat_U01;
         break;
      case '-':
         gen->GetBits = &LagFibSousFloat_Bits;
         gen->GetU01 = &LagFibSousFloat_U01;
         break;
      }
   }

   for (j = 0; j < state->Lag; j++)
      X[state->Lag - j] = (S[j] & MASK32) / unif01_NORM32;

   state->X = X;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrLagFibFloat;

   return gen;
}

/*-------------------------------------------------------------------------*/

void umrg_DeleteLagFibFloat (unif01_Gen * gen)
{
   LagFibFloat_state *state;

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

static unsigned long LagFibAdd_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   state->X[state->r] = (state->X[state->r] + state->X[state->s]) & param->Mask;
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibAddLux_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   if (--state->RR == 0) {
      int j;
      state->RR = state->Lag;
      for (j = 0; j < param->Skip; j++) {
         state->X[state->r] = (state->X[state->r] + state->X[state->s]) & param->Mask;
         if (--state->r == 0)
            state->r = state->Lag;
         if (--state->s == 0)
            state->s = state->Lag;
      }
   }
   state->X[state->r] = (state->X[state->r] + state->X[state->s]) & param->Mask;
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static double LagFibAdd_U01 (void *vpar, void *vsta)
{
   return LagFibAdd_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static double LagFibAddLux_U01 (void *vpar, void *vsta)
{
   return LagFibAddLux_Bits (vpar, vsta) * unif01_INV32;
}


/**************************************************************************/

static unsigned long LagFibSub_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   if (param->Flag)
      state->X[state->r] = (state->X[state->r] - state->X[state->s]) & param->Mask;
   else
      state->X[state->r] = (state->X[state->s] - state->X[state->r]) & param->Mask;
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibSubLux_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   if (--state->RR == 0) {
      int j;
      state->RR = state->Lag;
      for (j = 0; j < param->Skip; j++) {
         if (param->Flag)
            state->X[state->r] = (state->X[state->r] - state->X[state->s]) &
                                  param->Mask;
         else
            state->X[state->r] = (state->X[state->s] - state->X[state->r]) &
                                  param->Mask;
         if (--state->r == 0)
            state->r = state->Lag;
         if (--state->s == 0)
            state->s = state->Lag;
      }
   }
   if (param->Flag)
      state->X[state->r] = (state->X[state->r] - state->X[state->s]) & param->Mask;
   else
      state->X[state->r] = (state->X[state->s] - state->X[state->r]) & param->Mask;
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static double LagFibSub_U01 (void *vpar, void *vsta)
{
   return LagFibSub_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static double LagFibSubLux_U01 (void *vpar, void *vsta)
{
   return LagFibSubLux_Bits (vpar, vsta) * unif01_INV32;
}


/**************************************************************************/

static unsigned long LagFibXor_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   state->X[state->r] = state->X[state->s] ^ state->X[state->r];
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibXorLux_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   if (--state->RR == 0) {
      int j;
      state->RR = state->Lag;
      for (j = 0; j < param->Skip; j++) {
         state->X[state->r] = (state->X[state->r] ^ state->X[state->s]);
         if (--state->r == 0)
            state->r = state->Lag;
         if (--state->s == 0)
            state->s = state->Lag;
      }
   }
   state->X[state->r] = (state->X[state->r] ^ state->X[state->s]);
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static double LagFibXor_U01 (void *vpar, void *vsta)
{
   return LagFibXor_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static double LagFibXorLux_U01 (void *vpar, void *vsta)
{
   return LagFibXorLux_Bits (vpar, vsta) * unif01_INV32;
}


/**************************************************************************/

static unsigned long LagFibMult_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   state->X[state->r] = (state->X[state->r] * state->X[state->s]) & param->Mask;
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static unsigned long LagFibMultLux_Bits (void *vpar, void *vsta)
{
   LagFib_param *param = vpar;
   LagFib_state *state = vsta;
   unsigned long temp;

   if (--state->RR == 0) {
      int j;
      state->RR = state->Lag;
      for (j = 0; j < param->Skip; j++) {
         state->X[state->r] = (state->X[state->r] * state->X[state->s]) & param->Mask;
         if (--state->r == 0)
            state->r = state->Lag;
         if (--state->s == 0)
            state->s = state->Lag;
      }
   }
   state->X[state->r] = (state->X[state->r] * state->X[state->s]) & param->Mask;
   if (param->LeftShift)
      temp = state->X[state->r] << param->b;
   else
      temp = state->X[state->r] >> param->b;
   if (--state->r == 0)
      state->r = state->Lag;
   if (--state->s == 0)
      state->s = state->Lag;
   return temp;
}

/*-------------------------------------------------------------------------*/

static double LagFibMult_U01 (void *vpar, void *vsta)
{
   return LagFibMult_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static double LagFibMultLux_U01 (void *vpar, void *vsta)
{
   return LagFibMultLux_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrLagFib (void *vsta)
{
   LagFib_state *state = vsta;
   int j;

   if (unif01_WrLongStateFlag) {
      printf ("S = {\n");
      for (j = 0; j < state->Lag; j++) {
         printf (" %12lu", state->X[state->r]);
         if (--state->r == 0)
            state->r = state->Lag;
         if (j < state->Lag - 1)
            printf (",");
         if ((j % 5) == 4)
            printf ("\n");
      }
      printf ("   }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-------------------------------------------------------------------------*/

unif01_Gen * umrg_CreateLagFib (int t, int k, int r, char Op, int Lux,
   unsigned long S[])
{
   unif01_Gen *gen;
   LagFib_param *param;
   LagFib_state *state;
   size_t leng;
   char name[LEN + 1];
   char chaine[2];
   int j;
   unsigned long *X;

   util_Assert (t > 0, "umrg_CreateLagFib:   t <= 0");
#ifdef IS_ULONG32
   util_Assert (t <= 32, "umrg_CreateLagFib:   t > 32");
#else
   util_Assert (t <= 64, "umrg_CreateLagFib:   t > 64");
#endif
   /* util_Assert (k > r, "umrg_CreateLagFib:   k <= r"); */
   util_Assert ((Op == '*') || (Op == '+') || (Op == '-') || (Op == 'x'),
      "umrg_CreateLagFib:   Op must be one of { +, -, *, x }");

   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (LagFib_param));
   state = util_Malloc (sizeof (LagFib_state));

   strcpy (name, "umrg_CreateLagFib:");
   addstr_Long (name, "   t = ", t);
   addstr_Long (name, ",   k = ", k);
   addstr_Long (name, ",   r = ", r);
   strcat (name, ",   Op = ");
   sprintf (chaine, "%c", Op);
   strcat (name, chaine);
   addstr_Long (name, ",   Lux = ", Lux);

   if (k >= r) {
      state->Lag = k;
      state->r = k;
      state->s = r;
      param->Flag = TRUE;
   } else {
      state->Lag = r;
      state->r = r;
      state->s = k;
      param->Flag = FALSE;
   }
   param->Skip = Lux - state->Lag;

   addstr_ArrayUlong (name, ",   S = ", state->Lag, S);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   param->b = t - 32;
   if (param->b <= 0) {
      param->b = -param->b;
      param->LeftShift = TRUE;
   } else {
      param->LeftShift = FALSE;
   }
   param->Mask = (1UL << t) - 1;
#ifndef IS_ULONG32
   if (64 == t)
      param->Mask = 0xffffffffffffffffUL;
#endif
   if (32 == t)
      param->Mask = 0xffffffffUL;

   if (param->Skip > 0) {
      X = util_Calloc ((size_t) Lux + 1, sizeof (unsigned long));
      state->RR = state->Lag;

      switch (Op) {
      case '*':
         gen->GetBits = &LagFibMultLux_Bits;
         gen->GetU01 = &LagFibMultLux_U01;
         break;
      case '-':
         gen->GetBits = &LagFibSubLux_Bits;
         gen->GetU01 = &LagFibSubLux_U01;
         break;
      case '+':
         gen->GetBits = &LagFibAddLux_Bits;
         gen->GetU01 = &LagFibAddLux_U01;
         break;
      case 'x':
         gen->GetBits = &LagFibXorLux_Bits;
         gen->GetU01 = &LagFibXorLux_U01;
         break;
      }

   } else {
      X = util_Calloc ((size_t) state->Lag + 1, sizeof (unsigned long));
      switch (Op) {
      case '*':
         gen->GetBits = &LagFibMult_Bits;
         gen->GetU01 = &LagFibMult_U01;
         break;
      case '-':
         gen->GetBits = &LagFibSub_Bits;
         gen->GetU01 = &LagFibSub_U01;
         break;
      case '+':
         gen->GetBits = &LagFibAdd_Bits;
         gen->GetU01 = &LagFibAdd_U01;
         break;
      case 'x':
         gen->GetBits = &LagFibXor_Bits;
         gen->GetU01 = &LagFibXor_U01;
         break;
      }
   }

   if (Op == '*') {
      int flag = 0;
      for (j = 0; j < state->Lag; j++) {
         X[state->Lag - j] = (S[j] & param->Mask) | 1;
         if (X[state->Lag - j] % 4 != 1)
            flag = 1;
      }
      /* Make sure that not all seeds are = 1 mod 4 */
      if (!flag) 
         X[1] = (X[1] + 2) & param->Mask;
   } else {
      for (j = 0; j < state->Lag; j++)
         X[state->Lag - j] = S[j] & param->Mask;
   }

   state->X = X;
   gen->param = param;
   gen->state = state;
   gen->Write = &WrLagFib;
   return gen;
}


/*-------------------------------------------------------------------------*/

void umrg_DeleteLagFib (unif01_Gen * gen)
{
   LagFib_state *state;

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
