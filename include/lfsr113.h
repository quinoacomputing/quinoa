// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

typedef unsigned long long lt;

typedef struct{
  unsigned z1,z2,z3,z4;
} lfsr113_state;

typedef struct{
  unsigned z[4] __attribute__ ((aligned(16)));
} lfsr113_sse_state;

void lfsr113_skipahead_(lfsr113_state* state,unsigned long long offset64,unsigned long long offset0);

void lfsr113_init_(lfsr113_state* state);

void lfsr113_init_sequence_(lfsr113_state* state,lt SequenceNumber);

void lfsr113_init_long_sequence_(lfsr113_state* state,lt SequenceNumber);

unsigned int lfsr113_generate_(lfsr113_state* state);

float lfsr113_generate_uniform_float_(lfsr113_state* state);

unsigned int lfsr113_sse_generate_(lfsr113_sse_state* state); // use this function only if CPU supports SSE4

void lfsr113_get_sse_state_(lfsr113_state* state,lfsr113_sse_state* sse_state);

void lfsr113_print_state_(lfsr113_state* state);

