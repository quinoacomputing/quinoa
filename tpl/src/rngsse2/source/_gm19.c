// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gm19_g 524287
#define gm19_halfg 262144
#define gm19_k 15
#define gm19_q 28
#define gm19_qg 14680036

unsigned gm19_PPPP[12] __attribute__ ((aligned(16))) =  
    {1048574,1048574,1048574,1048574,524286,524286,524286,524286,524287,524287,524287,524287};

typedef struct{
  unsigned xN[32] __attribute__ ((aligned(16))),
           xP[32] __attribute__ ((aligned(16)));
} gm19_state;

typedef gm19_state gm19_sse_state;

unsigned int gm19_sse_generate_(gm19_sse_state* state){
  unsigned output1,output2;
  asm volatile("movaps (%4),%%xmm7\n" \
      "movaps 32(%4),%%xmm6\n" \

      "movaps (%2),%%xmm0\n" \
      "movaps (%3),%%xmm1\n" \
      "movaps %%xmm0,(%3)\n" \
      "movaps %%xmm0,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm0\n" \
      "psubd  %%xmm1,%%xmm0\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "psubd  %%xmm2,%%xmm0\n" \
      "paddd  %%xmm1,%%xmm0\n" \

      "movaps %%xmm0,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm0\n" \
      "paddd  %%xmm1,%%xmm0\n" \
      "movaps %%xmm0,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm0\n" \
      "movaps %%xmm0,(%2)\n" \

      "movaps 16(%2),%%xmm3\n" \
      "movaps 16(%3),%%xmm1\n" \
      "movaps %%xmm3,16(%3)\n" \
      "movaps %%xmm3,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm3\n" \
      "psubd  %%xmm1,%%xmm3\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "psubd  %%xmm2,%%xmm3\n" \
      "paddd  %%xmm1,%%xmm3\n" \

      "movaps %%xmm3,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm3\n" \
      "paddd  %%xmm1,%%xmm3\n" \
      "movaps %%xmm3,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm3\n" \
      "movaps %%xmm3,16(%2)\n" \

      "movaps 32(%2),%%xmm4\n" \
      "movaps 32(%3),%%xmm1\n" \
      "movaps %%xmm4,32(%3)\n" \
      "movaps %%xmm4,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm4\n" \
      "psubd  %%xmm1,%%xmm4\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "psubd  %%xmm2,%%xmm4\n" \
      "paddd  %%xmm1,%%xmm4\n" \


      "movaps %%xmm4,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm4\n" \
      "paddd  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,32(%2)\n" \

      "movaps 48(%2),%%xmm5\n" \
      "movaps 48(%3),%%xmm1\n" \
      "movaps %%xmm5,48(%3)\n" \
      "movaps %%xmm5,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm5\n" \
      "psubd  %%xmm1,%%xmm5\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "psubd  %%xmm2,%%xmm5\n" \
      "paddd  %%xmm1,%%xmm5\n" \

      "movaps %%xmm5,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm5\n" \
      "paddd  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,48(%2)\n" \

      "psrld  $18,%%xmm0\n" \
      "psrld  $18,%%xmm3\n" \
      "psrld  $18,%%xmm4\n" \
      "psrld  $18,%%xmm5\n" \
      "packssdw %%xmm3,%%xmm0\n" \
      "packssdw %%xmm5,%%xmm4\n" \
      "packsswb %%xmm4,%%xmm0\n" \
      "psllw  $7,%%xmm0\n" \
      "pmovmskb %%xmm0,%0\n" \

      "movaps 64(%2),%%xmm0\n" \
      "movaps 64(%3),%%xmm1\n" \
      "movaps %%xmm0,64(%3)\n" \
      "movaps %%xmm0,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm0\n" \
      "psubd  %%xmm1,%%xmm0\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "paddd  %%xmm0,%%xmm0\n" \
      "psubd  %%xmm2,%%xmm0\n" \
      "paddd  %%xmm1,%%xmm0\n" \

      "movaps %%xmm0,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm0\n" \
      "paddd  %%xmm1,%%xmm0\n" \
      "movaps %%xmm0,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm0\n" \
      "movaps %%xmm0,64(%2)\n" \

      "movaps 80(%2),%%xmm3\n" \
      "movaps 80(%3),%%xmm1\n" \
      "movaps %%xmm3,80(%3)\n" \
      "movaps %%xmm3,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm3\n" \
      "psubd  %%xmm1,%%xmm3\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "paddd  %%xmm3,%%xmm3\n" \
      "psubd  %%xmm2,%%xmm3\n" \
      "paddd  %%xmm1,%%xmm3\n" \

      "movaps %%xmm3,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm3\n" \
      "paddd  %%xmm1,%%xmm3\n" \
      "movaps %%xmm3,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm3\n" \
      "movaps %%xmm3,80(%2)\n" \

      "movaps 96(%2),%%xmm4\n" \
      "movaps 96(%3),%%xmm1\n" \
      "movaps %%xmm4,96(%3)\n" \
      "movaps %%xmm4,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm4\n" \
      "psubd  %%xmm1,%%xmm4\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "paddd  %%xmm4,%%xmm4\n" \
      "psubd  %%xmm2,%%xmm4\n" \
      "paddd  %%xmm1,%%xmm4\n" \

      "movaps %%xmm4,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm4\n" \
      "paddd  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm4\n" \
      "movaps %%xmm4,96(%2)\n" \

      "movaps 112(%2),%%xmm5\n" \
      "movaps 112(%3),%%xmm1\n" \
      "movaps %%xmm5,112(%3)\n" \
      "movaps %%xmm5,%%xmm2\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm7,%%xmm5\n" \
      "psubd  %%xmm1,%%xmm5\n" \
      "paddd  %%xmm1,%%xmm1\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "paddd  %%xmm5,%%xmm5\n" \
      "psubd  %%xmm2,%%xmm5\n" \
      "paddd  %%xmm1,%%xmm5\n" \

      "movaps %%xmm5,%%xmm1\n" \
      "psrld  $19,%%xmm1\n" \
      "pand   %%xmm6,%%xmm5\n" \
      "paddd  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "pcmpgtd 16(%4),%%xmm1\n" \
      "pand   %%xmm6,%%xmm1\n" \
      "psubd  %%xmm1,%%xmm5\n" \
      "movaps %%xmm5,112(%2)\n" \

      "psrld  $18,%%xmm0\n" \
      "psrld  $18,%%xmm3\n" \
      "psrld  $18,%%xmm4\n" \
      "psrld  $18,%%xmm5\n" \
      "packssdw %%xmm3,%%xmm0\n" \
      "packssdw %%xmm5,%%xmm4\n" \
      "packsswb %%xmm4,%%xmm0\n" \
      "psllw  $7,%%xmm0\n" \
      "pmovmskb %%xmm0,%1\n" \
      "shll $16,%1\n" \
      "addl %1,%0\n" \
      "":"=&r"(output1),"=&r"(output2):"r"(state->xN),"r"(state->xP),"r"(gm19_PPPP));
  return output1;
}

void gm19_get_sse_state_(gm19_state* state,gm19_sse_state* sse_state){
  int i; for(i=0;i<32;i++) {sse_state->xN[i]=state->xN[i]; sse_state->xP[i]=state->xP[i];}
}

unsigned gm19_CNext(unsigned N,unsigned P){
  return (gm19_qg+gm19_k*N-gm19_q*P)%gm19_g;
}

unsigned gm19_CNext2(unsigned N,unsigned P,unsigned myk,unsigned myq){
  unsigned long long curr1,curr2,curr3;
  curr1=(unsigned long long)myk*(unsigned long long)N;
  curr2=(unsigned long long)myq*(unsigned long long)P;
  curr3=((unsigned long long)myq*(unsigned long long)gm19_g+curr1-curr2)%gm19_g;
  return curr3;
}

unsigned gm19_GetNextN(unsigned x0,unsigned x1,unsigned n){ // returns x_{2^n}
  unsigned myk=gm19_k,myq=gm19_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gm19_CNext2(x,x0,myk,myq);
    myk=gm19_CNext2(myk,2,myk,myq);
    myq=gm19_CNext2(myq,0,myq,0);
  }
  return x;
}

unsigned gm19_GetNextAny(unsigned x0,unsigned x1,unsigned long long N){ // returns x_N
  unsigned long long i; unsigned xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm19_GetNextN(xp,xn,shift);
      xnnew=gm19_GetNextN(xn,gm19_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;
}

void gm19_skipahead_(gm19_state* state, unsigned long long offset){
  unsigned xn,xp,j; 
  for(j=0;j<32;j++){
    xp=gm19_GetNextAny(state->xP[j],state->xN[j],offset);
    xn=gm19_GetNextAny(state->xP[j],state->xN[j],offset+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

void gm19_init_(gm19_state* state){
   unsigned x0=514932,x1=127293,xp,xn,j;
   for(j=0;j<32;j++){
     xp=gm19_GetNextAny(x0,x1,8382841959ULL);
     xn=gm19_GetNextAny(x0,x1,8382841960ULL);
     state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
   }
}

void gm19_init_sequence_(gm19_state* state,unsigned SequenceNumber){
  gm19_init_(state);                     // 0 <= SequenceNumber < 1000;   length of each sequence <= 6*10^6
  gm19_skipahead_(state,6927047ULL*(unsigned long long)SequenceNumber);
}

unsigned int gm19_generate_(gm19_state* state){
  int i; unsigned temp,sum=0,bit=1;
  for(i=0;i<32;i++){
    temp=gm19_CNext(state->xN[i],state->xP[i]);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum += ((temp<gm19_halfg)?0:bit); bit*=2;
  }
  return sum;
}

float gm19_generate_uniform_float_(gm19_state* state){
  int i; unsigned temp,sum=0,bit=1;
  for(i=0;i<32;i++){
    temp=gm19_CNext(state->xN[i],state->xP[i]);
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum += ((temp<gm19_halfg)?0:bit); bit*=2;
  }
  return (float) sum/(float) 4294967296;
}

void gm19_print_state_(gm19_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%u",state->xN[i]%gm19_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%u",state->xP[i]%gm19_g); printf((i<31)?",":"}\n\n");}
}

void gm19_print_sse_state_(gm19_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<32;i++) {printf("%u",state->xN[i]%gm19_g); printf((i<31)?",":"}\nxP={");}
    for(i=0;i<32;i++) {printf("%u",state->xP[i]%gm19_g); printf((i<31)?",":"}\n\n");}
}
