// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#define gm61_g     2305843009213693951ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[32] __attribute__ ((aligned(16))),
     xP[32] __attribute__ ((aligned(16)));
} gm61_state;

typedef gm61_state gm61_sse_state;

void gm61_skipahead_(gm61_state* state, lt offset64, lt offset0);

void gm61_init_(gm61_state* state);

void gm61_init_sequence_(gm61_state* state,lt SequenceNumber);

void gm61_init_long_sequence_(gm61_state* state,lt SequenceNumber);

unsigned int gm61_generate_(gm61_state* state);

float gm61_generate_uniform_float_(gm61_state* state);

unsigned int gm61_sse_generate_(gm61_sse_state* state);

void gm61_get_sse_state_(gm61_state* state,gm61_sse_state* sse_state);

void gm61_print_state_(gm61_state* state);

void gm61_print_sse_state_(gm61_sse_state* state);
