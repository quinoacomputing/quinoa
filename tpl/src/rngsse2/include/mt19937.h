// (c) Copyright 2013 Lev Barash, Landau Institute for Theoretical Physics, Russian Academy of Sciences
// This is supplement to the paper:
// L.Yu. Barash, L.N. Shchur, "RNGSSELIB: Program library for random number generation. More generators, parallel streams of random numbers and Fortran compatibility".
// e-mail: barash @ itp.ac.ru (remove space)

#define mt19937_N 624

typedef struct{
  unsigned mt[mt19937_N];
  int mti;
} mt19937_state;

typedef struct{
  unsigned mt_aligned[3*mt19937_N+5] __attribute__ ((aligned(16)));
  unsigned out[3*mt19937_N+5] __attribute__ ((aligned(16)));
  unsigned *mt;
  int mti;
} mt19937_sse_state;

void mt19937_skipahead_(mt19937_state* state, unsigned long long offset0, unsigned offset_log);

void mt19937_init_(mt19937_state* state);

void mt19937_init_sequence_(mt19937_state* state, unsigned long long SequenceNumber);

unsigned mt19937_generate_(mt19937_state* state);

float mt19937_generate_uniform_float_(mt19937_state* state);

unsigned int mt19937_sse_generate_(mt19937_sse_state* state);

void mt19937_get_sse_state_(mt19937_state* state, mt19937_sse_state* sse_state);

void mt19937_print_state_(mt19937_state* state);
