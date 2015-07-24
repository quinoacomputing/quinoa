#ifndef UTEST_H
#define UTEST_H

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
     manage unit tests
     display results
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>

// ----- types

struct utest;

// ----- interface

#define T struct utest

T*   utest_alloc (FILE*);
void utest_free  (T*);
void utest_title (T*, const char*);
void utest_init  (T*, const char*);
void utest_fini  (T*);
void utest_stat  (T*);

// return the number of failed tests, prefer the UTEST macro
int  utest_test  (T*, int pass, const char *cond, const char *file, int line);

#undef T

// ----- macros

// assume T* to be named "utest"
#define UTEST(cond) utest_test(utest, (cond), #cond, __FILE__, __LINE__)

#endif
