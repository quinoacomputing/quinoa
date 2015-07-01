// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#include<stdio.h>

#define gm55_k 256
#define gm55_q 176
#define gm55_g       36028797018961904ULL
#define gm55_gdiv16  2251799813685119ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[8] __attribute__ ((aligned(16))), 
     xP[8] __attribute__ ((aligned(16)));
} gm55_state;

typedef gm55_state gm55_sse_state;

lt gm55_sse_Consts[8] __attribute__ ((aligned(16))) = {396316767208580944ULL,396316767208580944ULL,
   2064ULL,2064ULL,36028792732385279ULL,36028792732385279ULL,36028797018961904ULL,36028797018961904ULL};

unsigned int gm55_sse_generate_(gm55_sse_state* state){
  unsigned output;
  asm volatile("movaps (%3),%%xmm0\n" \

      "movaps (%2),%%xmm1\n" \
      "movaps (%1),%%xmm4\n" \
      "movaps %%xmm4,(%2)\n" \
      "psllq  $4,%%xmm4\n" \
      "paddq  %%xmm0,%%xmm4\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm4\n" \
      "movaps %%xmm4,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm4\n" \
      "paddq  %%xmm3,%%xmm4\n" \
      "psllq  $4,%%xmm4\n" \
      "movaps %%xmm4,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm4\n" \
      "movaps %%xmm4,(%1)\n" \

      "movaps 16(%2),%%xmm1\n" \
      "movaps 16(%1),%%xmm5\n" \
      "movaps %%xmm5,16(%2)\n" \
      "psllq  $4,%%xmm5\n" \
      "paddq  %%xmm0,%%xmm5\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm5\n" \
      "movaps %%xmm5,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm5\n" \
      "paddq  %%xmm3,%%xmm5\n" \
      "psllq  $4,%%xmm5\n" \
      "movaps %%xmm5,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm5\n" \
      "movaps %%xmm5,16(%1)\n" \

      "movaps 32(%2),%%xmm1\n" \
      "movaps 32(%1),%%xmm6\n" \
      "movaps %%xmm6,32(%2)\n" \
      "psllq  $4,%%xmm6\n" \
      "paddq  %%xmm0,%%xmm6\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm6\n" \
      "movaps %%xmm6,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm6\n" \
      "paddq  %%xmm3,%%xmm6\n" \
      "psllq  $4,%%xmm6\n" \
      "movaps %%xmm6,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm6\n" \
      "movaps %%xmm6,32(%1)\n" \

      "movaps 48(%2),%%xmm1\n" \
      "movaps 48(%1),%%xmm7\n" \
      "movaps %%xmm7,48(%2)\n" \
      "psllq  $4,%%xmm7\n" \
      "paddq  %%xmm0,%%xmm7\n" \
      "movaps %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psllq  $3,%%xmm1\n" \
      "paddq  %%xmm1,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm7\n" \
      "movaps %%xmm7,%%xmm2\n" \
      "psrlq  $51,%%xmm2\n" \
      "movaps %%xmm2,%%xmm3\n" \
      "psllq  $7,%%xmm3\n" \
      "paddq  %%xmm2,%%xmm3\n" \
      "psllq  $51,%%xmm2\n" \
      "psubq  %%xmm2,%%xmm7\n" \
      "paddq  %%xmm3,%%xmm7\n" \
      "psllq  $4,%%xmm7\n" \
      "movaps %%xmm7,%%xmm1\n" \
      "paddq  16(%3),%%xmm1\n" \
      "pshufd $245,%%xmm1,%%xmm3\n" \
      "pcmpgtd 32(%3),%%xmm3\n" \
      "pand    48(%3),%%xmm3\n" \
      "psubq   %%xmm3,%%xmm7\n" \
      "movaps %%xmm7,48(%1)\n" \

      "psrlq  $51,%%xmm4\n" \
      "psrlq  $51,%%xmm5\n" \
      "psrlq  $51,%%xmm6\n" \
      "psrlq  $51,%%xmm7\n" \
      "packssdw  %%xmm5,%%xmm4\n" \
      "packssdw  %%xmm7,%%xmm6\n" \
      "packssdw  %%xmm6,%%xmm4\n" \
      "packsswb  %%xmm4,%%xmm4\n" \
      "movaps  %%xmm4,%%xmm0\n" \
      "psrldq   $4,%%xmm0\n" \
      "pslld    $4,%%xmm0\n" \
      "pxor    %%xmm0,%%xmm4\n"
      "movd    %%xmm4,%0\n" \
      "":"=&r"(output):"r"(state->xN),"r"(state->xP),"r"(gm55_sse_Consts));
      return output;
}

void gm55_get_sse_state_(gm55_state* state,gm55_sse_state* sse_state){
  int i; for(i=0;i<8;i++) {sse_state->xN[i]=state->xN[i]; sse_state->xP[i]=state->xP[i];}
}

lt gm55_mod_g(lt x){ // returns x (mod g)
  lt F,G; F = (x>>55); G = x-(F<<55)+(2064*F);
  return ((G>=gm55_g) ? (G-gm55_g) : G);
}

lt gm55_MyMult(lt A,lt B){ // returns AB (mod gm55_g), where it is implied that A,B<gm55_g;
  lt A1,A0,B1,B0,curr,x,m;
  A1=A>>28; B1=B>>27; A0=A-(A1<<28); B0=B-(B1<<27);
  curr=2*A1*B0+B1*A0; m=curr>>28; x=curr-(m<<28);
  curr=(x<<27)+2064*m+(gm55_mod_g(129*A1*B1)<<4)+A0*B0;
  return gm55_mod_g(curr);
}

lt gm55_CNext2(lt N,lt P,lt myk,lt myq){   // returns (myk*N-myq*P) (mod gm55_g)
  lt curr1,curr2;
  curr1=gm55_MyMult(myk,N); curr2=gm55_MyMult(myq,P);
  if(curr1>=curr2) return (curr1-curr2); else return (gm55_g+curr1-curr2);
}

lt gm55_CNext(lt N,lt P){ // returns (256*N-176*P) (mod gm55_g)
  return gm55_mod_g((N<<8)+176*(gm55_g-P));
}

lt gm55_GetNextN(lt x0,lt x1,unsigned int n){ //returns x_{2^n}
  lt myk=gm55_k,myq=gm55_q,i,x=x1;
  for(i=0;i<n;i++){
    x=gm55_CNext2(x,x0,myk,myq);
    myk=gm55_CNext2(myk,2,myk,myq);
    myq=gm55_CNext2(myq,0,myq,0);
  }
  return x;
}

lt gm55_GetNextAny(lt x0,lt x1,lt N64,lt N0){ //N=2^64*N64+N0+1
  lt i,xp=x0,xn=x1,xpnew,xnnew,shift=0;
  i=N0; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm55_GetNextN(xp,xn,shift);
      xnnew=gm55_GetNextN(xn,gm55_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  i=N64; shift=64; while(i>0){
    if(i%2==1){                        // xp,xn ----> 2^shift
      xpnew=gm55_GetNextN(xp,xn,shift);
      xnnew=gm55_GetNextN(xn,gm55_CNext(xn,xp),shift);
      xp=xpnew; xn=xnnew;
    }
    i/=2; shift++;
  }
  return xp;                       // returns x_N, where N=2^64*N64+N0+1
}

void gm55_skipahead_(gm55_state* state, lt offset64, lt offset0){ // offset=offset64*2^64+offset0+1
  lt xn,xp,j; 
  for(j=0;j<8;j++){
    xp=gm55_GetNextAny(state->xP[j],state->xN[j],offset64,offset0);
    xn=gm55_GetNextAny(state->xP[j],state->xN[j],offset64,offset0+1);
    state->xP[j]=xp; state->xN[j]=xn;
  }
}

void gm55_init_(gm55_state* state){
  lt x0=gm55_mod_g(100152853817629549ULL),x1=gm55_mod_g(132388305121829306ULL),xp,xn,j;
  for(j=0;j<8;j++){
    xp=gm55_GetNextAny(x0,x1,7730941120ULL,2741045636569588180ULL);
    xn=gm55_GetNextAny(x0,x1,7730941120ULL,2741045636569588181ULL);
    state->xP[j]=xp; state->xN[j]=xn; x0=xp; x1=xn;
  }
}

void gm55_init_short_sequence_(gm55_state* state,lt SequenceNumber){ // 0 <= SequenceNumber < 10^18
  lt n1,n2;                                                                // length of each sequence < 10^10
  gm55_init_(state);
  n1=SequenceNumber/892447987; n2=SequenceNumber%892447987;
  gm55_skipahead_(state,n1,n1*4193950067);
  gm55_skipahead_(state,0,n2*20669825409); // thus we are skipping ahead (SequenceNumber*20669825409) numbers
}

void gm55_init_long_sequence_(gm55_state* state,lt SequenceNumber){ // 0 <= SequenceNumber < 4*10^9
  gm55_init_(state);                                                      // length of each sequence  <   10^20
  gm55_skipahead_(state,8*SequenceNumber,2699204111*SequenceNumber);
}

unsigned int gm55_generate_(gm55_state* state){
  unsigned int sum=0; int i; lt temp;
  for(i=0;i<8;i++){ 
    temp=gm55_mod_g(((state->xN[i])<<8)+gm55_q*(gm55_g-state->xP[i]));
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp/gm55_gdiv16)<<((i<4)?(8*i):(8*i-28)));
  }
  return sum;
}

float gm55_generate_uniform_float_(gm55_state* state){
  unsigned int sum=0; int i; lt temp;
  for(i=0;i<8;i++){ 
    temp=gm55_mod_g(((state->xN[i])<<8)+gm55_q*(gm55_g-state->xP[i]));
    state->xP[i]=state->xN[i]; state->xN[i]=temp;
    sum+= ((temp/gm55_gdiv16)<<((i<4)?(8*i):(8*i-28)));
  }
  return (float) sum/(float) 4294967296;
}

void gm55_print_state_(gm55_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<8;i++) {printf("%llu",state->xN[i]%gm55_g); printf((i<7)?",":"}\nxP={");}
    for(i=0;i<8;i++) {printf("%llu",state->xP[i]%gm55_g); printf((i<7)?",":"}\n\n");}
}

void gm55_print_sse_state_(gm55_sse_state* state){int i;
    printf("Generator State:\nxN={");
    for(i=0;i<8;i++) {printf("%llu",state->xN[i]%gm55_g); printf((i<7)?",":"}\nxP={");}
    for(i=0;i<8;i++) {printf("%llu",state->xP[i]%gm55_g); printf((i<7)?",":"}\n\n");}
}
