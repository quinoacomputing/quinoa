// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#ifndef mrg32k3a_h
#define mrg32k3a_h

typedef struct{
  unsigned x0,x1,x2,y0,y1,y2;
} mrg32k3a_state;

typedef struct{
  unsigned s[20] __attribute__ ((aligned(16))); 
} mrg32k3a_sse_state;

void mrg32k3a_skipahead_(mrg32k3a_state* state, unsigned long long offset128,
                                      unsigned long long offset64,unsigned long long offset0);

void mrg32k3a_init_(mrg32k3a_state* state);

void mrg32k3a_init_sequence_(mrg32k3a_state* state,unsigned long long SequenceNumber);

unsigned int mrg32k3a_generate_(mrg32k3a_state* state);

float mrg32k3a_generate_uniform_float_(mrg32k3a_state* state);

unsigned int mrg32k3a_sse_generate_(mrg32k3a_sse_state* state);

void mrg32k3a_get_sse_state_(mrg32k3a_state* state,mrg32k3a_sse_state* sse_state);

void mrg32k3a_print_state_(mrg32k3a_state* state);

#endif mrg32k3a_h
