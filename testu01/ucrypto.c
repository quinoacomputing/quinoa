/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ucrypto.c
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
#include "num.h"

#include "rijndael-alg-fst.h"
#include "tu01_sha1.h"
#include "ucrypto.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#define LEN  200                  /* Max length of strings */
#define MAX_SIZE 64               /* number of ASCII char's in a key */



/*-------------------------------------------------------------------------*/

typedef struct
{
   int Nr;                        /* Number of rounds */
   int Nk;                        /* Number of columns (key length in words) */
   int Nb;                        /* Number of bytes (block size in bytes) */
   int r;                         /* Number of initial bytes dropped */
   int s;                         /* Number of bytes kept */
} AES_param;

typedef struct
{
   uint32_t *K;                   /* Key */
   uint8_t *PT;                   /* Plain text */
   uint8_t *CT;                   /* Cypher text */
   int i;                         /* Byte number */
   int Nk;
   ucrypto_Mode mode;
} AES_state;


typedef struct
{
   SHA1_CTX context;
   int i;                         /* byte used for output */
   int r;                         /* Number of initial bytes dropped */
   int s;                         /* Number of bytes kept */
   ucrypto_Mode mode;
   unsigned char V[64];
   unsigned char output[20];
} SHA1_state;



/*============================ functions ==================================*/

static void increment8 (uint8_t A[], int n)
{
/* A is an array of n char's; it is viewed here as an 8n-bit integer.
   This function increments the integer by 1 mod 2^(8n). */
   int i;
   for (i = n-1; i >= 0; --i) {
      ++A[i];
#if CHAR_BIT > 8
      A[i] &= 0xff;
#endif
      if (A[i] > 0)
         return;
   }
}


/*-------------------------------------------------------------------------*/

static void increment32 (uint32_t A[], int n)
{
/* A is an array of n unsigned int's; it is viewed here as an 32n-bit integer.
   This function increments this integer by 1 mod 2^(32n). */
   int i;
   for (i = n-1; i >= 0; --i) {
      ++A[i];
#if UINT_MAX > 4294967295U
      A[i] &= 0xffffffffU;
#endif
      if (A[i] > 0)
         return;
   }
}


/*-------------------------------------------------------------------------*/
#if 0

int main ()
{
#if 0
   unsigned char A[] = { 25, 255, 254 };
   unsigned char B[] = { 111, 124 };
#else
   unsigned int A[] = { 25, 255, 0xfffffff9 };
#endif
   int i, j;
   for (j = 0; j < 10; ++j) {
      increment32 (A, 3);
      for (i = 0; i < 3; i++) {
         printf (" %13u", (unsigned int) A[i]);
      }
      printf ("\n");
   }
}
#endif


/*=========================================================================*/

static unsigned long AES_Bits (void *vpar, void *vsta)
{
   AES_param *param = vpar;
   AES_state *state = vsta;
   int k;
   unsigned long Z = 0;

   /* 4 bytes = one 32-bit random number */
   for (k = 0; k < 4; ++k) {
      if (state->i >= param->s) {
         int j;
         rijndaelEncrypt (state->K, param->Nr, state->PT, state->CT);

         switch (state->mode) {
         case ucrypto_OFB:
            for (j = 0; j < 16; j++)
               state->PT[j] = state->CT[j];
            break;
         case ucrypto_CTR:
            increment8 (state->PT, 16);    /* pt = pt + 1 mod 2^128 */
            break;
         case ucrypto_KTR:
            increment32 (state->K, param->Nk);   /* K = K + 1 mod 2^(32Nk) */
            break;
         default:
            util_Error ("ucrypto_CreateAES:   no such mode");
         }
         state->i = param->r;
      }
      Z <<= 8;
      Z |= (unsigned long) state->CT[state->i++];  
   }

   return Z;
}


/*-------------------------------------------------------------------------*/

static double AES_U01 (void *vpar, void *vsta)
{
   return AES_Bits (vpar, vsta) * unif01_INV32;
}


/*-------------------------------------------------------------------------*/

static void WrAES (void *vsta)
{
   AES_state *state = vsta;
   unsigned char A[32];
   int i;
   printf (" Char's are written as small decimal integers\n");

   switch (state->mode) {
   case ucrypto_OFB:
   case ucrypto_CTR:
      printf ("  T = {\n  ");
      for (i = 0; i < 16; i++) {
         printf ("  %3d", (unsigned int) state->PT[i]);
         if ((i + 1) % 8 == 0)
            printf ("\n  ");
      }
      break;

   case ucrypto_KTR:
      num_Uint2Uchar (A, state->K, state->Nk);
      printf ("  Key = {\n  ");
      for (i = 0; i < 4*state->Nk; i++) {
         printf ("  %3d", (unsigned int) A[i]);
         if ((i + 1) % 8 == 0)
            printf ("\n  ");
      }
      break;

   default:
      util_Error ("ucrypto_CreateAES:   no such mode");
   }
   printf ("}\n");
}


/*-------------------------------------------------------------------------*/

static void getStringMode (ucrypto_Mode mode, char *str)
/* Transforms the mode in a string */
{
   switch (mode) {
   case ucrypto_OFB:
      strcpy (str, "OFB");
      break;
   case ucrypto_CTR:
      strcpy (str, "CTR");
      break;
   case ucrypto_KTR:
      strcpy (str, "KTR");
      break;
   default:
      util_Error ("ucrypto_Mode:   no such case");
   }
}


/*-------------------------------------------------------------------------*/

unif01_Gen * ucrypto_CreateAES (unsigned char *Key, int KeyLen,
   unsigned char *Seed, ucrypto_Mode mode, int r, int s)
{
   unif01_Gen *gen;
   AES_param *param;
   AES_state *state;
   size_t len1;
   char name[LEN + 1] = {0};
   char str[16] = {0};
   int i;
   unsigned int D[64];

   util_Assert ((KeyLen == 128) || (KeyLen == 192) || (KeyLen == 256),
      "ucrypto_CreateAES:   klen must be in { 128, 192, 256 }");
   util_Assert (r < 16, "ucrypto_CreateAES:   r > 15");
   util_Assert (0 < s, "ucrypto_CreateAES:   s <= 0");
   util_Assert (s <= 16, "ucrypto_CreateAES:   s > 16");
   util_Assert (r + s <= 16, "ucrypto_CreateAES:   r + s > 16");
   gen = util_Malloc (sizeof (unif01_Gen));
   param = util_Malloc (sizeof (AES_param));
   state = util_Malloc (sizeof (AES_state));

   if (r < 0) r = 0;
   switch (KeyLen) {
   case 128:
      param->Nb = 16;
      param->Nk = 4;
      param->Nr = 10;
      break;
   case 192:
      param->Nb = 16;
      param->Nk = 6;
      param->Nr = 12;
      break;
   case 256:
      param->Nb = 16;
      param->Nk = 8;
      param->Nr = 14;
      break;
   default:
      util_Error ("ucrypto_CreateAES, klen:   no such case");
   }

   strncpy (name, "ucrypto_CreateAES:   mode = ", (size_t) LEN);
   getStringMode (mode, str);
   strncat (name, str, (size_t) LEN);
   addstr_Int (name, ",   r = ", r);
   addstr_Int (name, ",   s = ", s);
   addstr_Long (name, ",   klen = ", (long) KeyLen);
   for (i = 0; i < KeyLen / 8; i++)
      D[i] = Key[i];
   addstr_ArrayUint (name, ",   Key = ", KeyLen / 8, D);
   for (i = 0; i < param->Nb; i++)
      D[i] = Seed[i];
   addstr_ArrayUint (name, ",   Seed = ", param->Nb, D);
   len1 = strlen (name);
   gen->name = util_Calloc (len1 + 1, sizeof (char));
   strncpy (gen->name, name, len1);

   state->PT = util_Calloc ((size_t) param->Nb, sizeof (uint8_t));
   state->CT = util_Calloc ((size_t) param->Nb, sizeof (uint8_t));
   state->K = util_Calloc ((size_t) 4 * (param->Nr + 1), sizeof (uint32_t));
   rijndaelKeySetupEnc (state->K, Key, KeyLen);

   for (i = 0; i < param->Nb; i++)
      state->PT[i] = Seed[i];
   state->mode = mode;
   param->r = r;
   param->s = s + r;
   state->i = 16;
   state->Nk = param->Nk;
   gen->param = param;
   gen->state = state;
   gen->GetBits = &AES_Bits;
   gen->GetU01 = &AES_U01;
   gen->Write = &WrAES;
   return gen;
}


/*-------------------------------------------------------------------------*/

void ucrypto_DeleteAES (unif01_Gen * gen)
{
   AES_param *param;
   AES_state *state;

   if (NULL == gen)
      return;
   param = gen->param;
   state = gen->state;
   util_Free (state->CT);
   util_Free (state->PT);
   util_Free (state->K);
   gen->state = util_Free (gen->state);
   gen->param = util_Free (gen->param);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}


/*=========================================================================*/
#define SLEN 55

static unsigned long SHA1_Bits (void *junk, void *vsta)
{
   SHA1_state *state = vsta;
   int k;
   unsigned long Z = 0;

   /* 4 bytes = one 32-bit random number */
   for (k = 0; k < 4; ++k) {
      if (state->i >= state->s) {
         switch (state->mode) {
         case ucrypto_OFB:
            SHA1Init (&state->context);
            SHA1Update (&state->context, state->output, 20);
            SHA1Final (state->output, &state->context);
            break;
         case ucrypto_CTR:
            SHA1Init (&state->context);
            SHA1Update (&state->context, state->V, SLEN);
            SHA1Final (state->output, &state->context);
            increment8 (state->V, SLEN);      /* V = V + 1 mod 2^440 */
            break;
         default:
            util_Error ("ucrypto_CreateSHA1:   no such mode");
         }

         state->i = state->r;
      }

      Z <<= 8;
      Z |= (unsigned long) state->output[state->i++];  
   }

   return Z;
}

/*-------------------------------------------------------------------------*/

static double SHA1_U01 (void *vpar, void *vsta)
{
   return SHA1_Bits (vpar, vsta) * unif01_INV32;
}

/*-------------------------------------------------------------------------*/

static void WrSHA1 (void *vsta)
{
   SHA1_state *state = vsta;
   int i;
   printf (" Char's are written as small decimal integers\n");
   printf ("  T = {\n  ");

   switch (state->mode) {
   case ucrypto_OFB:
      for (i = 0; i < 20; i++) {
         printf ("  %3d", (unsigned int) state->output[i]);
         if ((i + 1) % 10 == 0)
            printf ("\n  ");
      }
      break;

   case ucrypto_CTR:
      for (i = 0; i < SLEN; i++) {
         printf ("  %3d", (unsigned int) state->V[i]);
         if ((i + 1) % 10 == 0)
            printf ("\n  ");
      }
      break;

   default:
      util_Error ("ucrypto_SHA1:   no such mode");
   }
   printf (" }\n");
}


/*-------------------------------------------------------------------------*/

unif01_Gen * ucrypto_CreateSHA1 (unsigned char *Seed, int len, ucrypto_Mode mode,
    int r, int s)
{
   unif01_Gen *gen;
   SHA1_state *state;
   int i;
   unsigned int D[SLEN];
   size_t len1;
   char name[LEN + 1] = {0};
   char str[16] = {0};

   util_Assert (4 == sizeof (uint32_t),
      "ucrypto_CreateSHA1:   uint32 must be exactly 32 bits wide");
   util_Assert (r < 20, "ucrypto_CreateSHA1:   r > 19");
   util_Assert (0 < s, "ucrypto_CreateSHA1:   s <= 0");
   util_Assert (s <= 20, "ucrypto_CreateSHA1:   s > 20");
   util_Assert (r + s <= 20, "ucrypto_CreateSHA1:   r + s > 20");
   gen = util_Malloc (sizeof (unif01_Gen));
   state = util_Malloc (sizeof (SHA1_state));
   memset (state, 0, sizeof (SHA1_state));
   if (r < 0) r = 0;

   strncpy (name, "ucrypto_CreateSHA1:   mode = ", (size_t) LEN);
   getStringMode (mode, str);
   strncat (name, str, (size_t) LEN);
   addstr_Int (name, ",   r = ", r);
   addstr_Int (name, ",   s = ", s);
   addstr_Int (name, ",   len = ", len);
   if (len > SLEN) len = SLEN;
   for (i = 0; i < len; i++)
      D[i] = Seed[i];
   addstr_ArrayUint (name, ",   Seed = ", len, D);
   len1 = strlen (name);
   gen->name = util_Calloc (len1 + 1, sizeof (char));
   strncpy (gen->name, name, len1);

   switch (mode) {
   case ucrypto_OFB:
      SHA1Init (&state->context);
      SHA1Update (&state->context, Seed, len);
      SHA1Final (state->output, &state->context);
      break;
   case ucrypto_CTR:
      for (i = 0; i < len; i++)
         state->V[i] = Seed[i];
      break;
   default:
      util_Error ("ucrypto_CreateSHA1:   no such mode");
   }

   state->mode = mode;
   state->r = r;
   state->s = s + r;
   state->i = 20;
   gen->param = NULL;
   gen->state = state;
   gen->GetBits = &SHA1_Bits;
   gen->GetU01 = &SHA1_U01;
   gen->Write = &WrSHA1;
   return gen;
}


/*-------------------------------------------------------------------------*/

void ucrypto_DeleteSHA1 (unif01_Gen * gen)
{
   SHA1_state *state;

   if (NULL == gen)
      return;
   state = gen->state;
   SHA1Final (state->output, &state->context);
   gen->state = util_Free (gen->state);
   gen->name = util_Free (gen->name);
   util_Free (gen);
}

#undef SLEN

/*=========================================================================*/
/* Generator ISAAC */

#include "ucryptoIS.c"
