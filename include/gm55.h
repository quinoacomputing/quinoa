// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#ifndef gm55_h
#define gm55_h

#define gm55_g       36028797018961904ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[8] __attribute__ ((aligned(16))), 
     xP[8] __attribute__ ((aligned(16)));
} gm55_state;

typedef gm55_state gm55_sse_state;

void gm55_skipahead_(gm55_state* state, lt offset64, lt offset0);

void gm55_init_(gm55_state* state);

void gm55_init_short_sequence_(gm55_state* state,lt SequenceNumber);

void gm55_init_long_sequence_(gm55_state* state,lt SequenceNumber);

unsigned int gm55_generate_(gm55_state* state);

float gm55_generate_uniform_float_(gm55_state* state);

unsigned gm55_sse_generate_(gm55_sse_state* state);

void gm55_get_sse_state_(gm55_state* state,gm55_sse_state* sse_state);

void gm55_print_state_(gm55_state* state);

void gm55_print_sse_state_(gm55_sse_state* state);

#endif // gm55_h
