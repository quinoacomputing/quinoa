#ifndef SLICE_H
#define SLICE_H

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
     manage slices
 
 o---------------------------------------------------------------------o
*/

#include "types.h"
#include "error.h"

// ----- types

struct slice {
  uint first, last, stride;
};

// ----- interface

#define T struct slice

static inline T
slice_initZero(void)
{
  return (T){ 0, 0, 1 };
}

static inline T
slice_initOne(void)
{
  return (T){ 1, 1, 1 };
}

static inline T
slice_initMax(void)
{
  return (T){ UINT_MAX, UINT_MAX, 1 };
}

static inline T
slice_initAll(void)
{
  return (T){ 0, UINT_MAX, 1 };
}

static inline T
slice_init(uint first)
{
  return (T){ first, first, 1 };
}

static inline T
slice_initSize(uint first, uint size)
{
  ensure(size, "invalid slice size");

  ullong last = first + (ullong)(size-1);
  if (size == UINT_MAX || last > UINT_MAX) last = UINT_MAX;

  return (T){ first, last, 1 };
}

static inline T
slice_initSizeStride(uint first, uint size, uint stride)
{
  ensure(size  , "invalid slice size"  );
  ensure(stride, "invalid slice stride");

  ullong last = first + (ullong)(size-1) * stride;
  if (size == UINT_MAX || last > UINT_MAX) last = UINT_MAX;

  return (T){ first, last, stride };
}

static inline T
slice_initLast(uint first, uint last)
{
  ensure(first <= last, "invalid range bounds");
  return (T){ first, last, 1 };
}

static inline T
slice_initLastStride(uint first, uint last, uint stride)
{
  ensure(first <= last, "invalid range bounds");
  ensure(stride       , "invalid range stride");

  return (T){ first, last, stride };
}

static inline uint
slice_first(const T* s)
{
  return s->first;
}

static inline uint
slice_stride(const T* s)
{
  return s->stride;
}

static inline uint
slice_last(const T* s)
{
  return s->last;
}

static inline uint
slice_end(const T* s)
{
  return s->last > (UINT_MAX-s->stride) ? UINT_MAX : s->last+s->stride;
}

static inline uint
slice_get(const T* s, uint n)
{
  return s->first + n * s->stride;
}

static inline uint
slice_sget(const T* s, uint n)
{
  ullong at = s->first + (ullong)n * s->stride;
  ensure(at <= s->last, "index out of range");
  return at;
}

static inline uint
slice_size(const T* s)
{
  uint size = (s->last - s->first)/s->stride + 1;
  return size ? size : UINT_MAX;
}

static inline uint
slice_width(const T* s)
{
  uint width = s->last - s->first + 1;
  return width ? width : UINT_MAX;
}

static inline bool
slice_isDense(const T* s)
{
  return s->stride == 1;
}

static inline bool
slice_isUnit(const T* s)
{
  return s->first == s->last;
}

static inline bool
slice_isInfinite(const T* s)
{
  return s->last == UINT_MAX;
}

static inline bool
slice_isFull(const T* s)
{
  return slice_first(s) == 0 && slice_isInfinite(s) && slice_isDense(s);
}

static inline bool
slice_isWithin(const T* s, uint n)
{
  return n >= s->first && n <= s->last;
}

static inline bool
slice_isEnum(const T* s, uint n)
{
  return slice_isDense(s) || (n - s->first) % s->stride == 0;
}

static inline bool
slice_isElem(const T* s, uint n)
{
  return slice_isWithin(s, n) && slice_isEnum(s, n);
}

#undef T

#endif
