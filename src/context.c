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
     manage contexts of constraints
     print, scan contexts from file
 
 o---------------------------------------------------------------------o
*/

#include <stdlib.h>
#include <assert.h>
#include <float.h>

#include "context.h"
#include "constraint.h"

#define T struct context
#define C struct constraint

// ----- types

struct context {
  // heaps of contraints
  const C **fut; // future, sorted
  const C **act; // active, sorted
  const C **row; // active, sorted
  int fut_n, act_n, row_n;

  // current status
  int row_u, row_i, col_i;
  bool sorted;

  // storage
  int dat_n, dat_sz;
  C dat[];
};

// ----- forward decl

static void
ut_trace(const T *cxt, int i, int j, const C* cst1, const C* cst2);

// ----- private (sort helpers)

static inline int
cmpCst (const void *cst1_, const void *cst2_)
{
  const C *cst1 = *(const void* const*)cst1_; 
  const C *cst2 = *(const void* const*)cst2_;
  const struct slice *row1 = &cst1->row;
  const struct slice *row2 = &cst2->row;
  const struct slice *col1 = &cst1->col;
  const struct slice *col2 = &cst2->col;

  // sorted by (> row start, > col start, > idx)
  return slice_first(row1) < slice_first(row2) ?  1
       : slice_first(row1) > slice_first(row2) ? -1
       : slice_first(col1) < slice_first(col2) ?  1
       : slice_first(col1) > slice_first(col2) ? -1
       : cst1 < cst2 ? 1 : -1;
}

static inline int
cmpRow (const C *cst1, const C *cst2)
{
  const struct slice *row1 = &cst1->row;
  const struct slice *row2 = &cst2->row;

  // sorted by (> row last, > idx)
  return slice_last(row1) < slice_last(row2) ?  1
       : slice_last(row1) > slice_last(row2) ? -1
       : cst1 < cst2 ? 1 : -1;
}

// ----- private (ctor & dtor helpers)

static inline void
context_setup (T *cxt)
{
  cxt->fut = malloc(3 * cxt->dat_n * sizeof *cxt->fut);
  ensure(cxt->fut, "out of memory");

  cxt->act = cxt->fut + cxt->dat_n;
  cxt->row = cxt->act + cxt->dat_n;
  cxt->fut_n = cxt->dat_n;
  cxt->act_n = cxt->row_n = 0;
  cxt->row_i = cxt->col_i = 0;

  for (int i = 0; i < cxt->fut_n; i++)
    cxt->fut[i] = cxt->dat+i;

  qsort(cxt->fut, cxt->fut_n, sizeof *cxt->fut, cmpCst);
  cxt->sorted = true;  
}

static void
context_teardown (T *cxt)
{
  free(cxt->fut);

  *cxt = (T) {
      .dat_n  = cxt->dat_n,
      .dat_sz = cxt->dat_sz
    };  
}

static inline T*
context_grow (T *cxt, int n)
{
  // enlarge on need
  if (n > cxt->dat_sz) {
    if (cxt->sorted == true)
      context_teardown(cxt);

    cxt = realloc(cxt, sizeof *cxt + n * sizeof *cxt->dat);
    ensure(cxt, "out of memory");
    cxt->dat_sz = n;
  }

  return cxt;  
}

static void
context_add0(T *cxt)
{
  // add rule #0: "* * abs=DBL_MIN"
  const C c = constraint_init(slice_initAll(), slice_initAll(), eps_init(eps_abs, DBL_MIN), -1, 0);

  context_add(cxt, &c);
}

// ----- private (eps helpers)

static inline void
context_updateAct (T *cxt, int row_i)
{
  trace("->updateAct row %d", row_i);
  int na = cxt->act_n;

  // remove obsolete constraints
  for (; cxt->act_n; --cxt->act_n) {
    const C *act = cxt->act[cxt->act_n-1];
    uint i = slice_last(&act->row);

    if (i >= (uint)row_i) {
      if (i < (uint)cxt->row_u) cxt->row_u = i;
      break;
    }
  }

  trace("%d obsolete constraints removed", na -= cxt->act_n);

  // select future constraints
  for (; cxt->fut_n; --cxt->fut_n) {
    const C *fut = cxt->fut[cxt->fut_n-1];
    uint i = slice_first(&fut->row);

    if (i > (uint)row_i) {         // not yet active
      if (i < (uint)cxt->row_u) cxt->row_u = i;
      break;
    }

    if (slice_last(&fut->row) < (uint)row_i) continue; // already obsolete

    // insert future constraint
    if (!cxt->act_n) *cxt->act = fut;
    else {
      const C **act = cxt->act+cxt->act_n-1;
      for (; act >= cxt->act; --act) {
        if (cmpRow(fut, *act) >= 0) break;
        act[1] = act[0];
      }
      act[1] = fut;
    }
    ++cxt->act_n;
  }

  trace("%d future constraints added", cxt->act_n-na);
  trace("<-updateAct row %d", row_i);
}

static inline void
context_setupRow (T *cxt, int row_i)
{
  trace("->setupRow row %d", row_i);
  cxt->row_n = 0;

  // select active constraints for this row
  for (int i = 0; i < cxt->act_n; ++i) {
    const C *act = cxt->act[i];
    if (!slice_isEnum(&act->row, row_i)) continue; // not active

    // action always dominates, unless hidden...
    if (act->eps.cmd >= eps_skip && !(act->eps.cmd & eps_alt)) {
      cxt->row[0] = act;
      cxt->row_n  = 1;
      break;
    }

    // add active constraint
    cxt->row[cxt->row_n++] = act;
  }

  trace("%d active constraints selected ([0] #%d, line %d)",
        cxt->row_n, cxt->row[0]->idx, cxt->row[0]->line);
  trace("<-setupRow row %d", row_i);
}

static inline const C*
context_setupCol (T *cxt, int col_i)
{
  trace("->setupCol col %d", col_i);
  const C *cst = 0;

  // select last-added active constraint for this col
  for (int i = 0; i < cxt->row_n; ++i) {
    const C *act = cxt->row[i];
    if (act > cst && !(act->eps.cmd & eps_alt) && 
        (act->eps.cmd >= eps_skip || slice_isElem(&act->col, col_i))) cst = act;
  }

  trace("constraint #%d (line %d) selected", cst->idx, cst->line);
  trace("<-setupCol col %d", col_i);

  return cst;
}

static inline const C*
context_getIncCst (T *cxt, int row_i, int col_i)
{
  const C *cst = 0;

  ensure(row_i >= cxt->row_i, "obsolete row");

  if (row_i > cxt->row_i) {
    if (row_i >= cxt->row_u)
      context_updateAct(cxt, row_i);       // update active constraints

    context_setupRow(cxt, row_i);          // setup constraints for this row

    cxt->row_i = row_i;
    cxt->col_i = 0;
  }

  ensure(col_i >= cxt->col_i, "obsolete column");

  cst = context_setupCol(cxt, col_i);      // setup constraints for this col

  cxt->col_i = col_i;

  return cst;
}

static inline const C*
context_getAtCst (T *cxt, int row_i, int col_i)
{
  const C *cur = cxt->dat+cxt->dat_n-1;
  const C *cst = 0;

  // select last-added active constraint, brute force...
  for (; cur >= cxt->dat; --cur)
    if (!(cur->eps.cmd & eps_alt) && slice_isElem(&cur->row, row_i)) {
      if (cur->eps.cmd >= eps_skip) return cur;
      if (slice_isElem(&cur->col, col_i)) { cst = cur--; break; }
    }

  // check for pending actions
  for (; cur >= cxt->dat; --cur)
    if (cur->eps.cmd >= eps_skip && !(cur->eps.cmd & eps_alt) && slice_isElem(&cur->row, row_i))
      return cur;

  return cst;
}

// -----------------------------------------------------------------------------
// ----- interface
// -----------------------------------------------------------------------------

T*
context_alloc (int n)
{
  enum { min_alloc = 128 };

  if (n < min_alloc) n = min_alloc;

  T *cxt = malloc(sizeof *cxt + n * sizeof *cxt->dat);
  ensure(cxt, "out of memory");

  *cxt = (T) { .dat_sz = n };

  context_add0(cxt);

  return cxt;
}

void
context_clear (T *cxt)
{
  assert(cxt);
  context_teardown(cxt);
  cxt->dat_n = 0;
  context_add0(cxt);
}

void
context_free (T *cxt)
{
  assert(cxt);
  context_teardown(cxt);
  free(cxt);
}

T*
context_add (T *cxt, const C *cst)
{
  assert(cxt && cst);

  // check if in use
  if (cxt->sorted == true)
    context_teardown(cxt);

  // check for storage space
  if (cxt->dat_n == cxt->dat_sz)
    cxt = context_grow(cxt, 1.75*cxt->dat_sz); // +75%

  // add constraint and setup index
  cxt->dat[cxt->dat_n] = *cst;
  cxt->dat[cxt->dat_n].idx = cxt->dat_n;

  // check for alternate qualifier, set onfail on previous rule
  if (cst->eps.cmd & eps_alt) {
    assert(cxt->dat_n > 0);
    int cmd = cxt->dat[cxt->dat_n-1].eps.cmd | eps_onfail;
    cxt->dat[cxt->dat_n-1].eps.cmd = (enum eps_cmd)cmd;
  }

  cxt->dat_n++;

  return cxt;
}

void
context_onfail(T *cxt, const C* cst)
{
  assert(cst->idx > 0);

  // clear alt qualifier of previous rule
  int cmd = cxt->dat[cst->idx-1].eps.cmd & ~eps_alt;
  cxt->dat[cst->idx-1].eps.cmd = (enum eps_cmd)cmd;
}

const C*
context_getAt (T *cxt, int row, int col)
{
  assert(cxt);
  ensure(row > 0, "null row");

  // check if ready for use
  if (cxt->sorted == false)
    context_setup(cxt);

  return context_getAtCst(cxt, row, col);
}

const C*
context_getInc (T *cxt, int row, int col)
{
  assert(cxt);
  ensure(row > 0, "null row");

  // check if ready for use
  if (cxt->sorted == false)
    context_setup(cxt);

  return context_getIncCst(cxt, row, col);
}

const C*
context_getIdx (const T *cxt, int idx)
{
  assert(cxt);
  return idx < cxt->dat_n ? cxt->dat+idx : 0;
}

int
context_findIdx (const T *cxt, const C *cst)
{
  assert(cxt && cst);
  return cst->idx >= 0   && cst->idx < cxt->dat_n     ? cst->idx :
         cst >= cxt->dat && cst < cxt->dat+cxt->dat_n ? cst-cxt->dat : -1;
}

int
context_findLine (const T *cxt, const C *cst)
{
  assert(cxt && cst);
  return cxt->dat_n > 0 && cst->line > 0 &&
         cst->line <= cxt->dat[cxt->dat_n-1].line ? cst->line : -1;
}

T*
context_scan(T *cxt, FILE *fp)
{
  int row = 1;
  C cst;

  while(!feof(fp)) {
    constraint_scan(&cst, fp, &row);
    if (cst.eps.cmd != eps_invalid)
      cxt = context_add(cxt, &cst);
  }

  return cxt;
}

void
context_print(const T *cxt, FILE *fp)
{
  const C *c;

  for (int i = 0; (c = context_getIdx(cxt, i)) != 0; i++) {
    fprintf(fp,"[#%d:%d] ", i, c->line);
    constraint_print(c, fp);
    putc('\n', fp);
  }
}

#undef T
#undef C

// -----------------------------------------------------------------------------
// ----- testsuite
// -----------------------------------------------------------------------------

#ifndef NTEST

#include "utest.h"

#define T struct context
#define C struct constraint

enum { NROW = 5, NCOL = 5 };

// ----- debug

static void
ut_trace(const T *cxt, int i, int j, const C* cst1, const C* cst2)
{
  fprintf(stderr, "(%d,%d)\n", i, j);
  if (cst1) {
    fprintf(stderr, "[%d].1: ", context_findIdx(cxt, cst1));
    constraint_print(cst1, stderr);
    putc('\n', stderr);
  }
  if (cst2) {
    fprintf(stderr, "[%d].2: ", context_findIdx(cxt, cst2));
    constraint_print(cst2, stderr);
    putc('\n', stderr);
  }
  fprintf(stderr, "{F} ");
  for(int k = 0; k < cxt->fut_n; k++)
    fprintf(stderr, "%d ", context_findIdx(cxt, cxt->fut[k]));

  fprintf(stderr, "\n{A} ");
  for(int k = 0; k < cxt->act_n; k++)
    fprintf(stderr, "%d ", context_findIdx(cxt, cxt->act[k]));

  fprintf(stderr, "\n{R} ");
  for(int k = 0; k < cxt->row_n; k++)
    fprintf(stderr, "%d ", context_findIdx(cxt, cxt->row[k]));

  putc('\n', stderr);
}

// ----- teardown

static T*
ut_teardown(T *cxt)
{
  context_clear(cxt);
  return cxt;
}

// ----- test

#if 0
/* tests for speed:
   - context_getAt  is much simpler but has bad complexity
   - context_getInc is (almost) always faster and much more stable
*/

static void
ut_testAt(struct utest *utest, T* cxt, int i, int j)
{
  const C* cst = context_getAt(cxt, i, j);
  UTEST(cst || !cst);
}

static void
ut_testInc(struct utest *utest, T* cxt, int i, int j)
{
  const C* cst = context_getInc(cxt, i, j);
  UTEST(cst || !cst);
}
#endif

static void
ut_testNul(struct utest *utest, T* cxt, int i, int j)
{
  const C* cst1 = context_getAt (cxt, i, j);
  const C* cst2 = context_getInc(cxt, i, j);

  UTEST(cst1 == cst2 &&
        cst1 == context_getIdx(cxt,0) &&
           0 == context_findIdx(cxt, cst1));

  if (cst1 != cst2)
    ut_trace(cxt, i, j, cst1, cst2);
}

static void
ut_testEqu(struct utest *utest, T* cxt, int i, int j)
{
  const C* cst1 = context_getAt (cxt, i, j);
  const C* cst2 = context_getInc(cxt, i, j);

  UTEST( (cst1 == cst2 ||
          (cst1 && (cst1->eps.cmd & eps_skip) &&
           cst2 && (cst2->eps.cmd & eps_skip)) ) );

  if (cst1 != cst2)
    ut_trace(cxt, i, j, cst1, cst2);
}

// ----- setup

static T*
ut_setup1(T *cxt)
{
  C cst;
  struct eps eps = eps_init(eps_dig, 1);

  // 3  2
  cst = constraint_init(slice_init(3), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);

  return cxt;
}

static T*
ut_setup2(T *cxt)
{
  C cst;
  struct eps eps = eps_init(eps_dig, 2);

  // 1-2  2
  cst = constraint_init(slice_initLast(1, 2), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 3    2
  cst = constraint_init(slice_init(3), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 4-5  2
  cst = constraint_init(slice_initSize(4, 2), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);

  return cxt;
}

static T*
ut_setup3(T *cxt)
{
  C cst;
  struct eps eps = eps_init(eps_dig, 3);

  // 2  1-2
  cst = constraint_init(slice_init(2), slice_initSize(1, 2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 2  3
  cst = constraint_init(slice_init(2), slice_init(3), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 2  4-5
  cst = constraint_init(slice_init(2), slice_initSize(4, 2), eps, -1, 0);
  cxt = context_add(cxt, &cst);

  return cxt;
}

static T*
ut_setup4(T *cxt)
{
  C cst;
  struct eps eps = eps_init(eps_dig , 4);
  struct eps skp = eps_init(eps_skip, 4);

  // 2    4-5
  cst = constraint_init(slice_init(2), slice_initSize(4, 2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 1-2  2
  cst = constraint_init(slice_initSize(1, 2), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 2    3
  cst = constraint_init(slice_init(2), slice_init(3), skp, -1, 0);
  cxt = context_add(cxt, &cst);
  // 3    2
  cst = constraint_init(slice_init(3), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 2    1-2
  cst = constraint_init(slice_init(2), slice_initSize(1, 2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 4-5  2
  cst = constraint_init(slice_initSize(4, 2), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 1-3  2-4
  cst = constraint_init(slice_initSize(1, 3), slice_initSize(2, 3), eps, -1, 0);
  cxt = context_add(cxt, &cst);

  return cxt;
}

static T*
ut_setup5(T *cxt)
{
  C cst;
  struct eps eps = eps_init(eps_dig , 5);
  struct eps skp = eps_init(eps_skip, 5);

  // 2-4/2    4-5
  cst = constraint_init(slice_initSizeStride(2, 2, 2), slice_initSize(4, 2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 4-5      2
  cst = constraint_init(slice_initSize(4, 2), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 1-3      2-4
  cst = constraint_init(slice_initSize(1, 3), slice_initSize(2, 3), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 1-6/2    2
  cst = constraint_init(slice_initSizeStride(1, 2, 3), slice_init(2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 2        3-7/2
  cst = constraint_init(slice_init(2), slice_initSizeStride(3, 2, 2), eps, -1, 0);
  cxt = context_add(cxt, &cst);
  // 3        2
  cst = constraint_init(slice_init(3), slice_init(2), skp, -1, 0);
  cxt = context_add(cxt, &cst);
  // 2        1-5/2
  cst = constraint_init(slice_init(2), slice_initSizeStride(1, 2, 2), eps, -1, 0);
  cxt = context_add(cxt, &cst);

  return cxt;
}

static T*
ut_setup6(T *cxt)
{
  C cst;
  struct eps eps = eps_init(eps_dig , 6);

  ut_setup3(cxt);
  ut_setup2(cxt);

  // 1-5    1-5
  cst = constraint_init(slice_initSize(1, 5), slice_initSize(1, 5), eps, -1, 0);
  cxt = context_add(cxt, &cst);

  return cxt;
}

static T*
ut_setup7(T *cxt)
{
  ut_setup3(cxt);
  ut_setup5(cxt);
  ut_setup2(cxt);
  ut_setup4(cxt);
  ut_setup1(cxt);
  return cxt;
}

// ----- unit tests

static struct spec {
  const char *name;
  T*        (*setup)   (T*);
  void      (*test )   (struct utest*, T*, int, int);
  T*        (*teardown)(T*);
} spec[] = {
  { "no constraint",                        0        , ut_testNul, ut_teardown },
  { "single constraint",                    ut_setup1, ut_testEqu, ut_teardown },
  { "multiple row constraints",             ut_setup2, ut_testEqu, ut_teardown },
  { "multiple column constraints",          ut_setup3, ut_testEqu, ut_teardown },
  { "overlapping constraints",              ut_setup4, ut_testEqu, ut_teardown },
  { "overlapping strided constraints",      ut_setup5, ut_testEqu, ut_teardown },
  { "sparse mixed strided constraints",     ut_setup6, ut_testEqu, ut_teardown },
  { "many mixed strided constraints",       ut_setup7, ut_testEqu, ut_teardown },
  { "no constraint (after use)",            0        , ut_testNul, ut_teardown }
};
enum { spec_n = sizeof spec/sizeof *spec };

// ----- interface

void
context_utest(struct utest *ut)
{
  assert(ut);
  T *cxt = context_alloc(0);

  utest_title(ut, "Context");

  for (int k = 0; k < spec_n; k++) {
    utest_init(ut, spec[k].name);
    if (spec[k].setup) cxt = spec[k].setup(cxt);

    spec[k].test(ut, cxt, 1, 1); // idempotent

    for (int i = 1; i <= NROW; i++)
    for (int j = 1; j <= NCOL; j++)
      spec[k].test(ut, cxt, i, j);
 
    spec[k].test(ut, cxt, NROW, NCOL); // idempotent

    if (spec[k].teardown) cxt = spec[k].teardown(cxt);
    utest_fini(ut);
  }

  context_free(cxt);
}

#endif
