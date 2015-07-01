// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gq58x1_k 8
#define gq58x1_q 48
#define gq58x1_g       288230374541099008ULL
#define gq58x1_halfg   144115187270549504ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[32] __attribute__ ((aligned(16))),
     xP[32] __attribute__ ((aligned(16)));
} gq58x1_state;

typedef gq58x1_state gq58x1_sse_state;

lt gq58x1_sse_Consts[10] __attribute__ ((aligned(16))) =
   {13835057977972752384ULL,13835057977972752384ULL,1610612736ULL,1610612736ULL,
    288230371923853311ULL,288230371923853311ULL,288230374541099008ULL,288230374541099008ULL,
    805306368ULL,805306368ULL};

unsigned int gq58x1_sse_generate_(gq58x1_sse_state* state){
  unsigned output1,output2;
  asm volatile("movaps (%4),%%xmm0\n" \

      "movaps (%3),%%xmm1\n" \
      "movaps (%2),%%xmm4\n" \
      "movaps %%xmm4,(%3)\n" \
      "psllq  $3,%%xmm4\n" \
      "paddq  %%xmm0,%%xmm4\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "paddq  %%xmm3,%%xmm4\n" \

      "movaps %%xmm4,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm4\n" \
      "movaps %%xmm4,(%2)\n" \
      "paddq  64(%4),%%xmm4\n" \

      "movaps 16(%3),%%xmm1\n" \
      "movaps 16(%2),%%xmm5\n" \
      "movaps %%xmm5,16(%3)\n" \
      "psllq  $3,%%xmm5\n" \
      "paddq  %%xmm0,%%xmm5\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "paddq  %%xmm3,%%xmm5\n" \

      "movaps %%xmm5,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,16(%2)\n" \
      "paddq  64(%4),%%xmm5\n" \

      "movaps 32(%3),%%xmm1\n" \
      "movaps 32(%2),%%xmm6\n" \
      "movaps %%xmm6,32(%3)\n" \
      "psllq  $3,%%xmm6\n" \
      "paddq  %%xmm0,%%xmm6\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "paddq  %%xmm3,%%xmm6\n" \

      "movaps %%xmm6,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,32(%2)\n" \
      "paddq  64(%4),%%xmm6\n" \

      "movaps 48(%3),%%xmm1\n" \
      "movaps 48(%2),%%xmm7\n" \
      "movaps %%xmm7,48(%3)\n" \
      "psllq  $3,%%xmm7\n" \
      "paddq  %%xmm0,%%xmm7\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "paddq  %%xmm3,%%xmm7\n" \

      "movaps %%xmm7,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,48(%2)\n" \
      "paddq  64(%4),%%xmm7\n" \

      "psrlq  $57,%%xmm4\n" \
      "psrlq  $57,%%xmm5\n" \
      "psrlq  $57,%%xmm6\n" \
      "psrlq  $57,%%xmm7\n" \
      "packssdw  %%xmm5,%%xmm4\n" \
      "packssdw  %%xmm7,%%xmm6\n" \
      "packssdw  %%xmm6,%%xmm4\n" \

      "movaps 64(%3),%%xmm1\n" \
      "movaps 64(%2),%%xmm5\n" \
      "movaps %%xmm5,64(%3)\n" \
      "psllq  $3,%%xmm5\n" \
      "paddq  %%xmm0,%%xmm5\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "paddq  %%xmm3,%%xmm5\n" \

      "movaps %%xmm5,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,64(%2)\n" \
      "paddq  64(%4),%%xmm5\n" \

      "movaps 80(%3),%%xmm1\n" \
      "movaps 80(%2),%%xmm6\n" \
      "movaps %%xmm6,80(%3)\n" \
      "psllq  $3,%%xmm6\n" \
      "paddq  %%xmm0,%%xmm6\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "paddq  %%xmm3,%%xmm6\n" \

      "movaps %%xmm6,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,80(%2)\n" \
      "paddq  64(%4),%%xmm6\n" \

      "movaps 96(%3),%%xmm1\n" \
      "movaps 96(%2),%%xmm7\n" \
      "movaps %%xmm7,96(%3)\n" \
      "psllq  $3,%%xmm7\n" \
      "paddq  %%xmm0,%%xmm7\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "paddq  %%xmm3,%%xmm7\n" \

      "movaps %%xmm7,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,96(%2)\n" \
      "paddq  64(%4),%%xmm7\n" \

      "movaps 112(%3),%%xmm1\n" \
      "movaps 112(%2),%%xmm2\n" \
      "movaps %%xmm2,112(%3)\n" \
      "psllq  $3,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm2\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm3,%%xmm2\n" \

      "movaps %%xmm2,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm2\n" \
      "movaps %%xmm2,112(%2)\n" \
      "paddq  64(%4),%%xmm2\n" \

      "psrlq  $57,%%xmm5\n" \
      "psrlq  $57,%%xmm6\n" \
      "psrlq  $57,%%xmm7\n" \
      "psrlq  $57,%%xmm2\n" \
      "packssdw  %%xmm6,%%xmm5\n" \
      "packssdw  %%xmm2,%%xmm7\n" \
      "packssdw  %%xmm7,%%xmm5\n" \
      "packsswb  %%xmm5,%%xmm4\n" \
      "psllw  $7,%%xmm4\n" \
      "pmovmskb  %%xmm4,%0\n" \




      "movaps 128(%3),%%xmm1\n" \
      "movaps 128(%2),%%xmm4\n" \
      "movaps %%xmm4,128(%3)\n" \
      "psllq  $3,%%xmm4\n" \
      "paddq  %%xmm0,%%xmm4\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "paddq  %%xmm3,%%xmm4\n" \

      "movaps %%xmm4,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm4\n" \
      "movaps %%xmm4,128(%2)\n" \
      "paddq  64(%4),%%xmm4\n" \

      "movaps 144(%3),%%xmm1\n" \
      "movaps 144(%2),%%xmm5\n" \
      "movaps %%xmm5,144(%3)\n" \
      "psllq  $3,%%xmm5\n" \
      "paddq  %%xmm0,%%xmm5\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "paddq  %%xmm3,%%xmm5\n" \

      "movaps %%xmm5,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,144(%2)\n" \
      "paddq  64(%4),%%xmm5\n" \

      "movaps 160(%3),%%xmm1\n" \
      "movaps 160(%2),%%xmm6\n" \
      "movaps %%xmm6,160(%3)\n" \
      "psllq  $3,%%xmm6\n" \
      "paddq  %%xmm0,%%xmm6\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "paddq  %%xmm3,%%xmm6\n" \

      "movaps %%xmm6,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,160(%2)\n" \
      "paddq  64(%4),%%xmm6\n" \

      "movaps 176(%3),%%xmm1\n" \
      "movaps 176(%2),%%xmm7\n" \
      "movaps %%xmm7,176(%3)\n" \
      "psllq  $3,%%xmm7\n" \
      "paddq  %%xmm0,%%xmm7\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "paddq  %%xmm3,%%xmm7\n" \

      "movaps %%xmm7,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,176(%2)\n" \
      "paddq  64(%4),%%xmm7\n" \

      "psrlq  $57,%%xmm4\n" \
      "psrlq  $57,%%xmm5\n" \
      "psrlq  $57,%%xmm6\n" \
      "psrlq  $57,%%xmm7\n" \
      "packssdw  %%xmm5,%%xmm4\n" \
      "packssdw  %%xmm7,%%xmm6\n" \
      "packssdw  %%xmm6,%%xmm4\n" \

      "movaps 192(%3),%%xmm1\n" \
      "movaps 192(%2),%%xmm5\n" \
      "movaps %%xmm5,192(%3)\n" \
      "psllq  $3,%%xmm5\n" \
      "paddq  %%xmm0,%%xmm5\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "paddq  %%xmm3,%%xmm5\n" \

      "movaps %%xmm5,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,192(%2)\n" \
      "paddq  64(%4),%%xmm5\n" \

      "movaps 208(%3),%%xmm1\n" \
      "movaps 208(%2),%%xmm6\n" \
      "movaps %%xmm6,208(%3)\n" \
      "psllq  $3,%%xmm6\n" \
      "paddq  %%xmm0,%%xmm6\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "paddq  %%xmm3,%%xmm6\n" \

      "movaps %%xmm6,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,208(%2)\n" \
      "paddq  64(%4),%%xmm6\n" \

      "movaps 224(%3),%%xmm1\n" \
      "movaps 224(%2),%%xmm7\n" \
      "movaps %%xmm7,224(%3)\n" \
      "psllq  $3,%%xmm7\n" \
      "paddq  %%xmm0,%%xmm7\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "paddq  %%xmm3,%%xmm7\n" \

      "movaps %%xmm7,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,224(%2)\n" \
      "paddq  64(%4),%%xmm7\n" \

      "movaps 240(%3),%%xmm1\n" \
      "movaps 240(%2),%%xmm2\n" \
      "movaps %%xmm2,240(%3)\n" \
      "psllq  $3,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "psllq  $4,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm2\n" \
      "psllq  $1,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $58,%%xmm1\n" \
      "psllq  $29,%%xmm1\n" \
      "movaps %%xmm1,%%xmm3\n" \
      "psllq  $1,%%xmm3\n" \
      "paddq  %%xmm1,%%xmm3\n" \
      "psllq  $29,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm3,%%xmm2\n" \

      "movaps %%xmm2,%%xmm1\n" \
      "paddq  16(%4),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%4),%%xmm3\n" \
      "pand    48(%4),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm2\n" \
      "movaps %%xmm2,240(%2)\n" \
      "paddq  64(%4),%%xmm2\n" \

      "psrlq  $57,%%xmm5\n" \
      "psrlq  $57,%%xmm6\n" \
      "psrlq  $57,%%xmm7\n" \
      "psrlq  $57,%%xmm2\n" \
      "packssdw  %%xmm6,%%xmm5\n" \
      "packssdw  %%xmm2,%%xmm7\n" \
      "packssdw  %%xmm7,%%xmm5\n" \
      "packsswb  %%xmm5,%%xmm4\n" \
      "psllw  $7,%%xmm4\n" \
      "pmovmskb  %%xmm4,%1\n" \

      "shll $16,%1\n" \
      "addl %0,%1\n" \


      "":"=&r"(output1),"=&r"(output2):"r"(state->xN),"r"(state->xP),"r"(gq58x1_sse_Consts));
      return output2;
}

void gq58x1_get_sse_state_(gq58x1_state* state,gq58x1_sse_state* sse_state){
  int i; for(i=0;i<32;i++) {sse_state->xN[i]=state->xN[i]; sse_state->xP[i]=state->xP[i];}
}

lt gq58x1_mod_g(lt x){ // returns x (mod g)
  lt F,G; F = (x>>58); G = x-(F<<58)+(F<<29)+(F<<30);
  return ((G>=gq58x1_g) ? (G-gq58x1_g) : G);
}

lt gq58x1_MyMult(lt A,lt B){ // returns AB (mod gq58x1_g), where it is implied that A,B<gq58x1_g;
  lt A1,A0,B1,B0,curr,x,m;
  A1=A>>32; B1=B>>32; A0=A-(A1<<32)+(12*A1); B0=B-(B1<<32)+(12*B1);
  if(A0>>32) {A0-=4294967284ULL; A1++;}
  if(B0>>32) {B0-=4294967284ULL; B1++;}
  curr=A1*B0+B1*A0; m=curr>>26; x=curr-(m<<26);
  curr=((3*m+(x<<4))<<28)+(gq58x1_g-12*x)+(144*A1*B1)+(gq58x1_mod_g(A0*B0));
  return gq58x1_mod_g(curr);
}

lt gq58x1_CNext2(lt N,lt P,lt myk,lt myq){   // returns (myk*N-myq*P) (mod gq58x1_g)
  lt curr1,curr2;
  curr1=gq58x1_MyMult(myk,N); curr2=gq58x1_MyMult(myq,P);
  if(curr1>=curr2) return (curr1-curr2); else return (gq58x1_g+curr1-curr2);
}

lt gq58x1_CNext(lt N,lt P){ // returns (8N-48P) (mod gq58x1_g)
  return gq58x1_mod_g((N+6*(gq58x1_g-P))<<3);
}

lt gq58x1_GetNextN(lt x0,lt x1,unsigned int n){ //returns x_{2^n}
  lt myk=gq58x1_k,myq=gq58x1_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gq58x1_CNext2(x,x0,myk,myq);
    myk=gq58x1_CNext2(myk,2,myk,myq);
    myq=gq58x1_CNext2(myq,0,myq,0);
  }
  return x;
}

lt gq58x1_GetNextAny(lt x0,lt x1,lt N64,lt N0){ //N=2^64*N64+N0+1
  lt i,xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N0; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gq58x1_GetNextN(xp,xn,shift);
      xnnew=gq58x1_GetNextN(xn,gq58x1_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  i=N64; shift=64; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gq58x1_GetNextN(xp,xn,shift);
      xnnew=gq58x1_GetNextN(xn,gq58x1_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;                       // returns x_N, where N=2^64*N64+N0+1
}

void gq58x1_skipahead_(gq58x1_state* state, lt offset64, lt offset0){ // offset=offset64*2^64+offset0+1
  lt xn,xp,j; 
  for(j=0;j<32;j++){
    xp=gq58x1_GetNextAny(state->xP[j],state->xN[j],offset64,offset0);
    xn=gq58x1_GetNextAny(state->xP[j],state->xN[j],offset64,offset0+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

void gq58x1_init_(gq58x1_state* state){
  lt x0=100152853817629549ULL,x1=132388305121829306ULL,xp,xn,j;
  for(j=0;j<32;j++){
    xp=gq58x1_GetNextAny(x0,x1,0,9003035591327124ULL);
    xn=gq58x1_GetNextAny(x0,x1,0,9003035591327125ULL);
    state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
  }
}

void gq58x1_init_short_sequence_(gq58x1_state* state,unsigned SequenceNumber){
  gq58x1_init_(state);                     // 0 <= SequenceNumber < 10^8;   length of each sequence <= 8*10^7
  gq58x1_skipahead_(state,0,82927047ULL*(unsigned long long)SequenceNumber);
}

void gq58x1_init_medium_sequence_(gq58x1_state* state,unsigned SequenceNumber){
  gq58x1_init_(state);                     // 0 <= SequenceNumber < 10^6;   length of each sequence <= 8*10^9
  gq58x1_skipahead_(state,0,8799201913ULL*(unsigned long long)SequenceNumber);
}

void gq58x1_init_long_sequence_(gq58x1_state* state,unsigned SequenceNumber){
  gq58x1_init_(state);                     // 0 <= SequenceNumber < 10^4;   length of each sequence <= 8*10^11
  gq58x1_skipahead_(state,0,828317697521ULL*(unsigned long long)SequenceNumber);
}

unsigned int gq58x1_generate_(gq58x1_state* state){
  unsigned int sum=0,bit=1; int i; lt temp;
  for(i=0;i<32;i++){
    temp=gq58x1_mod_g((state->xN[i]+6*(gq58x1_g-state->xP[i]))<<3);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp<gq58x1_halfg)?0:bit); bit*=2;
  }
  return sum;
}

float gq58x1_generate_uniform_float_(gq58x1_state* state){
  unsigned int sum=0,bit=1; int i; lt temp;
  for(i=0;i<32;i++){
    temp=gq58x1_mod_g((state->xN[i]+6*(gq58x1_g-state->xP[i]))<<3);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp<gq58x1_halfg)?0:bit); bit*=2;
  }
  return (float) sum/(float) 4294967296;
}

void gq58x1_print_state_(gq58x1_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%llu",state->xN[i]%gq58x1_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%llu",state->xP[i]%gq58x1_g); printf((i<31)?",":"}\n\n");}
}

void gq58x1_print_sse_state_(gq58x1_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%llu",state->xN[i]%gq58x1_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%llu",state->xP[i]%gq58x1_g); printf((i<31)?",":"}\n\n");}
}
