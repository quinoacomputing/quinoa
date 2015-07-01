// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#define gm31_g     2147483647

typedef struct{
  unsigned xN[32],xP[32];
} gm31_state;

typedef struct{
  unsigned xN[64] __attribute__ ((aligned(16))),
           xP[64] __attribute__ ((aligned(16)));
} gm31_sse_state;

void gm31_skipahead_(gm31_state* state, unsigned long long offset);

void gm31_init_(gm31_state* state);

void gm31_init_short_sequence_(gm31_state* state,unsigned SequenceNumber);

void gm31_init_medium_sequence_(gm31_state* state,unsigned SequenceNumber);

void gm31_init_long_sequence_(gm31_state* state,unsigned SequenceNumber);

unsigned int gm31_generate_(gm31_state* state);

float gm31_generate_uniform_float_(gm31_state* state);

unsigned int gm31_sse_generate_(gm31_sse_state* state);

void gm31_get_sse_state_(gm31_state* state,gm31_sse_state* sse_state);

void gm31_print_state_(gm31_state* state);

void gm31_print_sse_state_(gm31_sse_state* state);
