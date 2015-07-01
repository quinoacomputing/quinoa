// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#define gm19_g 524287

typedef struct{
  unsigned xN[32] __attribute__ ((aligned(16))),
           xP[32] __attribute__ ((aligned(16)));
} gm19_state;

typedef gm19_state gm19_sse_state;

void gm19_skipahead_(gm19_state* state, unsigned long long offset);

void gm19_init_(gm19_state* state);

void gm19_init_sequence_(gm19_state* state,unsigned SequenceNumber);

unsigned int gm19_generate_(gm19_state* state);

float gm19_generate_uniform_float_(gm19_state* state);

unsigned int gm19_sse_generate_(gm19_sse_state* state);

void gm19_get_sse_state_(gm19_state* state,gm19_sse_state* sse_state);

void gm19_print_state_(gm19_state* state);

void gm19_print_sse_state_(gm19_sse_state* state);
