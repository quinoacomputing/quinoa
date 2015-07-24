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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "utils.h"
#include "utest.h"
#include "error.h"

// ----- constants

#ifndef MAXUTESTKEEP
#define MAXUTESTKEEP 25
#endif

enum { UTEST_KEEP = MAXUTESTKEEP };

// ----- types

struct utest {
  // output file
  FILE*       out;

  // history of failed tests
  const char* fail_cond[UTEST_KEEP];
  const char* fail_file[UTEST_KEEP];
  int         fail_line[UTEST_KEEP];
  int         fail_n;

  // current test info
  const char* test_str;
  int         pass_cnt;
  int         fail_cnt;
  double      test_t0;

  // tests stat
  int         total_pass;
  int         total_fail;
  double      total_time;
};

// ----- public

struct utest*
utest_alloc(FILE *out)
{
  struct utest *ut = malloc(sizeof *ut);
  ensure(ut, "out of memory");

  ut->out = out ? out : stdout;
  ut->fail_n     = 0;
  ut->total_pass = 0;
  ut->total_fail = 0;
  ut->total_time = 0;

  return ut;
}

void
utest_free(struct utest *ut)
{
  assert(ut);
  free(ut);
}

void
utest_title(struct utest *ut, const char *str)
{
  fprintf(ut->out, " [ %s ]\n", str);
}

void
utest_init(struct utest *ut, const char *str)
{
  assert(ut);

  fprintf(ut->out, " + %-50s ", str);
  ut->test_str  = str;
  ut->pass_cnt  = 0;
  ut->fail_cnt  = 0;
  ut->fail_n    = 0;
  ut->test_t0   = clock();
}

int
utest_test(struct utest *ut, int pass, const char *cond, const char *file, int line)
{
  assert(ut);

  if (pass)
    ++ut->pass_cnt;
  else {
    ++ut->fail_cnt;
    if (ut->fail_n < UTEST_KEEP) {
      const char *p = strrchr(file, '/');
      if (p) file = p+1;
      ut->fail_cond[ut->fail_n  ] = cond;
      ut->fail_file[ut->fail_n  ] = file;
      ut->fail_line[ut->fail_n++] = line;
    }
  }

  return ut->fail_cnt;
}

void
utest_fini(struct utest *ut)
{
  assert(ut);

  double t1 = clock();
  double t = (t1 - ut->test_t0) / CLOCKS_PER_SEC;
  
  fprintf(ut->out,
          "(%.2f s) - %3d/%3d : %s\n",
          t,
          ut->pass_cnt,
          ut->pass_cnt+ut->fail_cnt,
          ut->fail_cnt ? fail_str : pass_str);

  for (int i = 0; i < ut->fail_n; i++)
    fprintf(ut->out,
            "   - (%s,%d) %s\n",
            ut->fail_file[i],
            ut->fail_line[i],
            ut->fail_cond[i]);

  ut->total_time += t;
  ut->total_pass += ut->pass_cnt;
  ut->total_fail += ut->fail_cnt;
}

void
utest_stat(struct utest *ut)
{
  assert(ut);

  fprintf(ut->out,
  " = %5d total, %5d passed, %5d failed            (%.2f s)             %s\n",
    ut->total_pass+ut->total_fail,
    ut->total_pass,
    ut->total_fail,
    ut->total_time,
    ut->total_fail ? fail_str : pass_str);
}

