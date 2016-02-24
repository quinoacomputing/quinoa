// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#ifndef gq58x4_h
#define gq58x4_h

#define gq58x4_g       288230374541099008ULL

typedef unsigned long long lt;

typedef struct{
  lt xN[8] __attribute__ ((aligned(16))),
     xP[8] __attribute__ ((aligned(16)));
} gq58x4_state;

typedef gq58x4_state gq58x4_sse_state;

void gq58x4_skipahead_(gq58x4_state* state, lt offset64, lt offset0);

void gq58x4_init_(gq58x4_state* state);

void gq58x4_init_short_sequence_(gq58x4_state* state,unsigned SequenceNumber);

void gq58x4_init_medium_sequence_(gq58x4_state* state,unsigned SequenceNumber);

void gq58x4_init_long_sequence_(gq58x4_state* state,unsigned SequenceNumber);

unsigned int gq58x4_generate_(gq58x4_state* state);

float gq58x4_generate_uniform_float_(gq58x4_state* state);

unsigned int gq58x4_sse_generate_(gq58x4_sse_state* state);

void gq58x4_get_sse_state_(gq58x4_state* state,gq58x4_sse_state* sse_state);

void gq58x4_print_state_(gq58x4_state* state);

void gq58x4_print_sse_state_(gq58x4_sse_state* state);

#endif // gq58x4_h
