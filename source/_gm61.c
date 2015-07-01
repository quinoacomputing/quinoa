// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gm61_k 24
#define gm61_q 74
#define gm61_g     2305843009213693951ULL
#define gm61_sixg  13835058055282163706ULL
#define gm61_halfg 1152921504606846976ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[32] __attribute__ ((aligned(16))),
     xP[32] __attribute__ ((aligned(16)));
} gm61_state;

typedef gm61_state gm61_sse_state;

lt gm61_Consts[2] __attribute__ ((aligned(16))) = {2305843009213693951UL,2305843009213693951UL};

unsigned int gm61_sse_generate_(gm61_sse_state* state){
  unsigned output1,output2;
  asm volatile(

      "movaps (%3),%%xmm0\n" \
      "movaps (%2),%%xmm3\n" \
      "movaps %%xmm3,(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm4\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm4\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,(%2)\n" \

      "movaps 16(%3),%%xmm0\n" \
      "movaps 16(%2),%%xmm3\n" \
      "movaps %%xmm3,16(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm5\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,16(%2)\n" \

      "movaps 32(%3),%%xmm0\n" \
      "movaps 32(%2),%%xmm3\n" \
      "movaps %%xmm3,32(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm6\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,32(%2)\n" \

      "movaps 48(%3),%%xmm0\n" \
      "movaps 48(%2),%%xmm3\n" \
      "movaps %%xmm3,48(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm7\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,48(%2)\n" \



      "psrlq  $60,%%xmm4\n" \
      "psrlq  $60,%%xmm5\n" \
      "psrlq  $60,%%xmm6\n" \
      "psrlq  $60,%%xmm7\n" \
      "packssdw  %%xmm5,%%xmm4\n" \
      "packssdw  %%xmm7,%%xmm6\n" \
      "packssdw  %%xmm6,%%xmm4\n" \


      "movaps 64(%3),%%xmm0\n" \
      "movaps 64(%2),%%xmm3\n" \
      "movaps %%xmm3,64(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm5\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,64(%2)\n" \


      "movaps 80(%3),%%xmm0\n" \
      "movaps 80(%2),%%xmm3\n" \
      "movaps %%xmm3,80(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm6\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,80(%2)\n" \

      "movaps 96(%3),%%xmm0\n" \
      "movaps 96(%2),%%xmm3\n" \
      "movaps %%xmm3,96(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm7\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,96(%2)\n" \

      "movaps 112(%3),%%xmm0\n" \
      "movaps 112(%2),%%xmm3\n" \
      "movaps %%xmm3,112(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm0\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm0\n" \
      "movaps %%xmm0,112(%2)\n" \


      "psrlq  $60,%%xmm5\n" \
      "psrlq  $60,%%xmm6\n" \
      "psrlq  $60,%%xmm7\n" \
      "psrlq  $60,%%xmm0\n" \
      "packssdw  %%xmm6,%%xmm5\n" \
      "packssdw  %%xmm0,%%xmm7\n" \
      "packssdw  %%xmm7,%%xmm5\n" \
      "packsswb  %%xmm5,%%xmm4\n" \
      "psllw  $7,%%xmm4\n" \
      "pmovmskb  %%xmm4,%0\n" \


      "movaps 128(%3),%%xmm0\n" \
      "movaps 128(%2),%%xmm3\n" \
      "movaps %%xmm3,128(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm4\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm4\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,128(%2)\n" \

      "movaps 144(%3),%%xmm0\n" \
      "movaps 144(%2),%%xmm3\n" \
      "movaps %%xmm3,144(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm5\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,144(%2)\n" \

      "movaps 160(%3),%%xmm0\n" \
      "movaps 160(%2),%%xmm3\n" \
      "movaps %%xmm3,160(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm6\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,160(%2)\n" \

      "movaps 176(%3),%%xmm0\n" \
      "movaps 176(%2),%%xmm3\n" \
      "movaps %%xmm3,176(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm7\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,176(%2)\n" \



      "psrlq  $60,%%xmm4\n" \
      "psrlq  $60,%%xmm5\n" \
      "psrlq  $60,%%xmm6\n" \
      "psrlq  $60,%%xmm7\n" \
      "packssdw  %%xmm5,%%xmm4\n" \
      "packssdw  %%xmm7,%%xmm6\n" \
      "packssdw  %%xmm6,%%xmm4\n" \



      "movaps 192(%3),%%xmm0\n" \
      "movaps 192(%2),%%xmm3\n" \
      "movaps %%xmm3,192(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm5\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm5\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,192(%2)\n" \


      "movaps 208(%3),%%xmm0\n" \
      "movaps 208(%2),%%xmm3\n" \
      "movaps %%xmm3,208(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm6\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm6\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm6\n" \
      "movaps %%xmm6,208(%2)\n" \

      "movaps 224(%3),%%xmm0\n" \
      "movaps 224(%2),%%xmm3\n" \
      "movaps %%xmm3,224(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm7\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm7\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm7\n" \
      "movaps %%xmm7,224(%2)\n" \

      "movaps 240(%3),%%xmm0\n" \
      "movaps 240(%2),%%xmm3\n" \
      "movaps %%xmm3,240(%3)\n" \
      "movaps (%4),%%xmm1\n" \
      "psubq  %%xmm0,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm0,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "paddq  %%xmm3,%%xmm3\n" \
      "paddq  %%xmm3,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm2\n" \
      "paddq  %%xmm0,%%xmm2\n" \
      "movaps %%xmm2,%%xmm1\n" \
      "psrlq  $61,%%xmm2\n" \
      "paddq  %%xmm2,%%xmm1\n" \
      "psllq  $61,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm1\n" \
      "movaps %%xmm1,%%xmm0\n" \
      "psrlq  $61,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm0\n" \
      "psllq  $61,%%xmm1\n" \
      "psubq  %%xmm1,%%xmm0\n" \
      "movaps %%xmm0,240(%2)\n" \


      "psrlq  $60,%%xmm5\n" \
      "psrlq  $60,%%xmm6\n" \
      "psrlq  $60,%%xmm7\n" \
      "psrlq  $60,%%xmm0\n" \
      "packssdw  %%xmm6,%%xmm5\n" \
      "packssdw  %%xmm0,%%xmm7\n" \
      "packssdw  %%xmm7,%%xmm5\n" \
      "packsswb  %%xmm5,%%xmm4\n" \
      "psllw  $7,%%xmm4\n" \
      "pmovmskb  %%xmm4,%1\n" \


      "shll $16,%1\n" \
      "addl %0,%1\n" \

      "":"=&r"(output1),"=&r"(output2):"r"(state->xN),"r"(state->xP),"r"(gm61_Consts));
      return output2;
}

void gm61_get_sse_state_(gm61_state* state,gm61_sse_state* sse_state){
  int i; for(i=0;i<32;i++) {sse_state->xN[i]=state->xN[i]; sse_state->xP[i]=state->xP[i];}
}

lt gm61_mod_g(lt x){ // returns x (mod g)
  lt F,G; F = (x>>61); G = x-(F<<61)+F;
  return ((G>=gm61_g) ? (G-gm61_g) : G);
}

lt gm61_CNext(lt N,lt P){ // returns (24N-74P) (mod g)
  lt curr1,curr2,curr3;
  curr1=P+P+P; curr1=gm61_mod_g(N+N+gm61_sixg-curr1-curr1);
  curr2=curr1+curr1+curr1; curr3=gm61_mod_g(curr2+curr2+gm61_g-P);
  return gm61_mod_g(curr3+curr3);
}

lt gm61_MyMult(lt A,lt B){ // returns AB (mod gm61_g), where it is implied that A,B<gm61_g;
  lt A1,A0,B1,B0,curr,x,m;
  A1=A>>32; B1=B>>32; A0=A-(A1<<32); B0=B-(B1<<32);
  curr=A1*B0+B1*A0; m=curr>>29; x=curr-(m<<29);
  curr=((x<<32)+m)+(8*A1*B1)+(gm61_mod_g(A0*B0));
  return gm61_mod_g(curr);
}

lt gm61_CNext2(lt N,lt P,lt myk,lt myq){   // returns (myk*N-myq*P) (mod gm61_g)
  lt curr1,curr2;
  curr1=gm61_MyMult(myk,N); curr2=gm61_MyMult(myq,P);
  if(curr1>=curr2) return (curr1-curr2); else return (gm61_g+curr1-curr2);
}

lt gm61_GetNextN(lt x0,lt x1,lt n){  // returns x_{2^n}
  lt myk=gm61_k,myq=gm61_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gm61_CNext2(x,x0,myk,myq);
    myk=gm61_CNext2(myk,2,myk,myq);
    myq=gm61_CNext2(myq,0,myq,0);
  }
  return x;
}

lt gm61_GetNextAny(lt x0,lt x1,lt N64,lt N0){ //N=2^64*N64+N0+1
  lt i,xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N0; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm61_GetNextN(xp,xn,shift);
      xnnew=gm61_GetNextN(xn,gm61_CNext2(xn,xp,gm61_k,gm61_q),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  i=N64; shift=64; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm61_GetNextN(xp,xn,shift);
      xnnew=gm61_GetNextN(xn,gm61_CNext2(xn,xp,gm61_k,gm61_q),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;                       // returns x_N, where N=2^64*N64+N0+1
}

void gm61_skipahead_(gm61_state* state, lt offset64, lt offset0){ // offset=offset64*2^64+offset0+1
  lt xn,xp,j; 
  for(j=0;j<32;j++){
    xp=gm61_GetNextAny(state->xP[j],state->xN[j],offset64,offset0);
    xn=gm61_GetNextAny(state->xP[j],state->xN[j],offset64,offset0+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

void gm61_init_(gm61_state* state){
  lt x0=463729159099860932ULL,x1=2069432424142811205ULL,xp,xn,j;
  for(j=0;j<32;j++){
    xp=gm61_GetNextAny(x0,x1,9007199254740991ULL,7371121435597082533ULL);
    xn=gm61_GetNextAny(x0,x1,9007199254740991ULL,7371121435597082534ULL);
    state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
  }
}

void gm61_init_sequence_(gm61_state* state,lt SequenceNumber){ // 0 <= SequenceNumber < 1.8*10^19
  lt n1,n2;                                                              // length of each sequence < 10^10
  gm61_init_(state);
  n1=SequenceNumber/892447987; n2=SequenceNumber%892447987;
  gm61_skipahead_(state,n1,n1*4193950067);
  gm61_skipahead_(state,0,n2*20669825409); // thus we are skipping ahead (SequenceNumber*20669825409) numbers
}

void gm61_init_long_sequence_(gm61_state* state,lt SequenceNumber){ // 0 <= SequenceNumber < 4*10^9
  gm61_init_(state);                                                      // length of each sequence  < 3*10^25
  gm61_skipahead_(state,2000000*SequenceNumber,2699204111*SequenceNumber);
}

unsigned int gm61_generate_(gm61_state* state){
  int i; lt temp; unsigned int sum=0,bit=1;
  for(i=0;i<32;i++){
    temp=gm61_CNext(state->xN[i],state->xP[i]);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum += ((temp<gm61_halfg)?0:bit); bit*=2;
  }
  return sum;
}

float gm61_generate_uniform_float_(gm61_state* state){
  int i; lt temp,sum=0,bit=1;
  for(i=0;i<32;i++){
    temp=gm61_CNext(state->xN[i],state->xP[i]);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum += ((temp<gm61_halfg)?0:bit); bit*=2;
  }
  return (float) sum/(float) 4294967296;
}

void gm61_print_state_(gm61_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%llu",state->xN[i]%gm61_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%llu",state->xP[i]%gm61_g); printf((i<31)?",":"}\n\n");}
}

void gm61_print_sse_state_(gm61_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%llu",state->xN[i]%gm61_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%llu",state->xP[i]%gm61_g); printf((i<31)?",":"}\n\n");}
}
