// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#define gq58x1_g       288230374541099008ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[32] __attribute__ ((aligned(16))),
     xP[32] __attribute__ ((aligned(16)));
} gq58x1_state;

typedef gq58x1_state gq58x1_sse_state;

void gq58x1_skipahead_(gq58x1_state* state, lt offset64, lt offset0);

void gq58x1_init_(gq58x1_state* state);

void gq58x1_init_short_sequence_(gq58x1_state* state,unsigned SequenceNumber);

void gq58x1_init_medium_sequence_(gq58x1_state* state,unsigned SequenceNumber);

void gq58x1_init_long_sequence_(gq58x1_state* state,unsigned SequenceNumber);

unsigned int gq58x1_generate_(gq58x1_state* state);

float gq58x1_generate_uniform_float_(gq58x1_state* state);

unsigned int gq58x1_sse_generate_(gq58x1_sse_state* state);

void gq58x1_get_sse_state_(gq58x1_state* state,gq58x1_sse_state* sse_state);

void gq58x1_print_state_(gq58x1_state* state);

void gq58x1_print_sse_state_(gq58x1_sse_state* state);
