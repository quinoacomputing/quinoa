#ifndef NDIFF_H
#define NDIFF_H

/*
 o---------------------------------------------------------------------o
 |
 | Ndiff
 |
 | Copyright (c) 2012+ laurent.deniau@cern.ch
 | Gnu General Public License
 |
 o---------------------------------------------------------------------o
  
   Purpose:
     numerical diff of files
     provides the main ndiff loop
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>
#include "types.h"

// ----- types

struct utest;
struct ndiff;
struct context;
struct constraint;

// ----- constrant

enum ndiff_options {
  ndiff_recycle_left = -1, ndiff_norecycle = 0, ndiff_recycle_right = 1
};

// ----- interface

#define T struct ndiff
#define C struct constraint

T*    ndiff_alloc    (FILE *lhs, FILE *rhs, struct context*, int n_, int r_);
void  ndiff_clear    (T*);
void  ndiff_free     (T*);
void  ndiff_option   (T*, const int *keep_, const int *blank_, const int *check_, const int *recycle_);
void  ndiff_result   (T*, FILE *lhs, FILE *rhs);

// high level API
void  ndiff_loop     (T*);

// low level API
int   ndiff_outLine  (T*);
int   ndiff_skipLine (T*);
int   ndiff_readLine (T*);
int   ndiff_fillLine (T*, const char *lhs, const char *rhs);

int   ndiff_gotoLine (T*, const C*);
int   ndiff_gotoNum  (T*, const C*);

int   ndiff_nextNum  (T*, const C*); // return 0 if no number is found
int   ndiff_testNum  (T*, const C*);

void  ndiff_getInfo  (const T*, int *row_, int *col_, int *cnt_, long *num_);
int   ndiff_feof     (const T*, int both);
int   ndiff_isempty  (const T*);

#undef T
#undef C

// ----- testsuite

#ifndef NTEST

void ndiff_utest (struct utest*);

#endif // NTEST
#endif
