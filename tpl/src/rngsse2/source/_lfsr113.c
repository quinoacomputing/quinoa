// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

typedef unsigned long long lt;

typedef struct{
  unsigned z1,z2,z3,z4;
} lfsr113_state;

typedef struct{
  unsigned z[4] __attribute__ ((aligned(16)));
} lfsr113_sse_state;

unsigned lfsr113_a[4] __attribute__ ((aligned(16)))={4294967294U,4294967288U,4294967280U,4294967168U};
int lfsr113_b[4] __attribute__ ((aligned(16)))={262144,4,128,8192};
int lfsr113_c[4] __attribute__ ((aligned(16)))={64,4,8192,8};

unsigned int lfsr113_sse_generate_(lfsr113_sse_state* state){ // here SSE4 instruction pblendw is used
  unsigned output;
  asm volatile("movaps (%1),%%xmm1\n" \
      "movaps (%2),%%xmm2\n" \
      "movaps (%4),%%xmm0\n" \
      "pand   %%xmm1,%%xmm2\n" \
      "pmulld (%3),%%xmm2\n" \
      "pmulld %%xmm1,%%xmm0\n" \
      "pxor   %%xmm0,%%xmm1\n" \
      "psrld  $12,%%xmm1\n" \
      "pblendw $192,%%xmm1,%%xmm3\n" \
      "psrld  $1,%%xmm1\n" \
      "pblendw $3,%%xmm1,%%xmm3\n" \
      "psrld  $8,%%xmm1\n" \
      "pblendw $48,%%xmm1,%%xmm3\n" \
      "psrld  $6,%%xmm1\n" \
      "pblendw $12,%%xmm1,%%xmm3\n" \
      "pxor   %%xmm2,%%xmm3\n" \
      "movaps %%xmm3,(%1)\n" \
      "pshufd $255,%%xmm3,%%xmm0\n" \
      "pshufd $170,%%xmm3,%%xmm1\n" \
      "pshufd $85,%%xmm3,%%xmm2\n" \
      "pxor %%xmm0,%%xmm3\n" \
      "pxor %%xmm1,%%xmm2\n" \
      "pxor %%xmm2,%%xmm3\n" \
      "pextrd $0,%%xmm3,%0\n" \
      "":"=&r"(output):"r"(state->z),"r"(lfsr113_a),"r"(lfsr113_b),"r"(lfsr113_c));
  return output;
}

void lfsr113_get_sse_state_(lfsr113_state* state,lfsr113_sse_state* sse_state){
  sse_state->z[0]=state->z1;
  sse_state->z[1]=state->z2;
  sse_state->z[2]=state->z3;
  sse_state->z[3]=state->z4;
}

unsigned lfsr113_SkipAheadRoundSingleBit(unsigned state,unsigned bit,unsigned p,unsigned q,unsigned s,int n){
  char arr0[64]; char arr[128][32];                              // this function skips ahead 2^n*s bits
  unsigned e,i,j,k1=p,k2,l1,l2,ki,TwoInE;
  for(j=0;j<p; j++) arr0[j] = (state>>(31-j)) & 1;
  for(j=p;j<64;j++) arr0[j] = arr0[j-p]^arr0[j-p+q];
  if(((1<<n)*s)<32){
    j=(1<<n)*s;  return (arr0[j+bit]<<(31-bit));
  }
  else{
    i=0; ki=s; while(ki<p){ ki*=2; i++;}
    TwoInE=1;
    for(e=0;e<=(n-i);e++){
      for(j=0;j<p;j++){
        if(j==0) l1=128; else {k1=j; l1=0; while(k1<p) {k1*=2; l1++;}}
        k2=j+q; l2=0; while(k2<p) {k2*=2; l2++;}
        arr[e][j]=(e<l1 ? arr0[j*TwoInE+bit] : arr[e-l1][k1-p])^( e<l2 ? arr0[(j+q)*TwoInE+bit] : arr[e-l2][k2-p]);
      }
      TwoInE*=2;
    }
    return (arr[n-i][ki-p]<<(31-bit));
  }
}

unsigned lfsr113_SkipAheadRoundSingle(unsigned state,unsigned p,unsigned q,unsigned s,int n){
  char arr0[64]; char arr[128][32];                              // this function skips ahead 2^n*s bits
  unsigned bit,e,i,j,k1=p,k2,l1,l2,ki,TwoInE; unsigned output=0;
  for(j=0;j<p; j++) arr0[j] = (state>>(31-j)) & 1;
  for(j=p;j<64;j++) arr0[j] = arr0[j-p]^arr0[j-p+q];
  if(((1<<n)*s)<32){
    j=(1<<n)*s; for(i=0;i<32;i++) output+= (arr0[j+i]<<(31-i)); return output;
  }
  else{
    i=0; ki=s; while(ki<p){ ki*=2; i++;}
    for(bit=0;bit<32;bit++){
      TwoInE=1;
      for(e=0;e<=(n-i);e++){
        for(j=0;j<p;j++){
          if(j==0) l1=128; else {k1=j; l1=0; while(k1<p) {k1*=2; l1++;}}
          k2=j+q; l2=0; while(k2<p) {k2*=2; l2++;}
          arr[e][j]=(e<l1 ? arr0[j*TwoInE+bit] : arr[e-l1][k1-p])^( e<l2 ? arr0[(j+q)*TwoInE+bit] : arr[e-l2][k2-p]);
        }
        TwoInE*=2;
      }
      output+=(arr[n-i][ki-p]<<(31-bit));
    }
    return output;
  }
}

void lfsr113_SkipAheadRound(lfsr113_state* state,int n){ // Skips Ahead 2^n
  state->z1=lfsr113_SkipAheadRoundSingle(state->z1,31,6,18,n);
  state->z2=lfsr113_SkipAheadRoundSingle(state->z2,29,2,2,n);
  state->z3=lfsr113_SkipAheadRoundSingle(state->z3,28,13,7,n);
  state->z4=lfsr113_SkipAheadRoundSingle(state->z4,25,3,13,n);
}

void lfsr113_skipahead_(lfsr113_state* state,unsigned long long offset64,unsigned long long offset0){
  unsigned long long i=offset0; int shift=0;
  while(i>0){
    if(i%2==1) lfsr113_SkipAheadRound(state,shift);
    i/=2; shift++;
  }
  i=offset64; shift=64;
  while(i>0){
    if(i%2==1) lfsr113_SkipAheadRound(state,shift);
    i/=2; shift++;
  }
}

void lfsr113_init_(lfsr113_state* state){
  state->z1=state->z2=state->z3=state->z4=12345;
}

void lfsr113_init_sequence_(lfsr113_state* state,lt SequenceNumber){ 
  lt n1,n2;                // 0 <= SequenceNumber < 3.8*10^18,  length of each sequence < 10^10
  lfsr113_init_(state);
  n1=SequenceNumber/892447987; n2=SequenceNumber%892447987;
  lfsr113_skipahead_(state,n1,n1*4193950067); // 20669825409*892447987 = 2^64 + 4193950067
  lfsr113_skipahead_(state,0,n2*20669825409); // thus we are skipping ahead (SequenceNumber*20669825409) numbers
}

void lfsr113_init_long_sequence_(lfsr113_state* state,lt SequenceNumber){
  lfsr113_init_(state);     // 0 <= SequenceNumber < 4*10^9.  length of each sequence  < 10^24
  lfsr113_skipahead_(state,100000*SequenceNumber,2699204111*SequenceNumber);
}

unsigned int lfsr113_generate_(lfsr113_state* state){
   unsigned b;
   b = ((state->z1 <<  6) ^ state->z1) >> 13;
   state->z1 = ((state->z1 & 4294967294U) << 18) ^ b;
   b = ((state->z2 <<  2) ^ state->z2) >> 27;
   state->z2 = ((state->z2 & 4294967288U) <<  2) ^ b;
   b = ((state->z3 << 13) ^ state->z3) >> 21;
   state->z3 = ((state->z3 & 4294967280U) <<  7) ^ b;
   b = ((state->z4 <<  3) ^ state->z4) >> 12;
   state->z4 = ((state->z4 & 4294967168U) << 13) ^ b;
   return (state->z1 ^ state->z2 ^ state->z3 ^ state->z4);
}

float lfsr113_generate_uniform_float_(lfsr113_state* state){
   unsigned b;
   b = ((state->z1 <<  6) ^ state->z1) >> 13;
   state->z1 = ((state->z1 & 4294967294U) << 18) ^ b;
   b = ((state->z2 <<  2) ^ state->z2) >> 27;
   state->z2 = ((state->z2 & 4294967288U) <<  2) ^ b;
   b = ((state->z3 << 13) ^ state->z3) >> 21;
   state->z3 = ((state->z3 & 4294967280U) <<  7) ^ b;
   b = ((state->z4 <<  3) ^ state->z4) >> 12;
   state->z4 = ((state->z4 & 4294967168U) << 13) ^ b;
   return (state->z1 ^ state->z2 ^ state->z3 ^ state->z4) *  2.3283064365386963e-10;
}

void lfsr113_print_state_(lfsr113_state* state){
    printf("Generator State: z1=%u, z2=%u, z3=%u, z4=%u\n",state->z1,state->z2,state->z3,state->z4);
}
