#ifndef UTILS_H
#define UTILS_H

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
     provides utilities
 
 o---------------------------------------------------------------------o
*/

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "args.h"

#define MkString(a) MkString_(a)
#define MkString_(a) #a

#define MkConcat(a,b) MkConcat_(a,b)
#define MkConcat_(a,b) a ## b

// extern constants

extern const char* pass_str;
extern const char* fail_str;

// extern functions

FILE*  open_file(const char* str, FILE **res_fp, int *idx, const char *ext, int optext, int required);
void  close_file(FILE *fp, int zip);

void  accum_summary(int total, int failed, long lines, long numbers);

// inline functions

#if !__STDC__ || __STDC_VERSION__ < 199901L
static inline int
isblank(int c)
{
  return c == ' ' || c == '\t';
}

static inline double
fmin(double a, double b)
{
  return a < b ? a : b;
}

static inline double
fmax(double a, double b)
{
  return a > b ? a : b;
}
#endif

static inline int
imin (int a, int b)
{
  return a<b ? a : b;
}

static inline int
imax (int a, int b)
{
  return a>b ? a : b;
}

static inline double 
pow10(int i)
{
  extern const double *const pow10_table99;
  return -100 < i && i < 100 ? pow10_table99[i] : pow(10, i);
}

// ----- public (read helpers)

static inline int
isComment(const char *buf)
{
  if (isspace(*option.cchr))
    while(isspace(*buf)) ++buf;

  return strchr(option.cchr, *buf) != 0;
}

static inline int
skipSpace (FILE *fp, int *i_)
{
  int c = 0, i = 0;

  while ((c = getc(fp)) != EOF) {
    if (!isspace(c)) break;
    i++;
  }

  if (i_) *i_ = i;

  return c; 
}

static inline int
skipLine (FILE *fp, int *i_)
{
  int c = 0, i = 0;

  while ((c = getc(fp)) != EOF) {
    if (c == '\n') break;                 // \n   : Unix, Linux, MacOSX
    if (c == '\r') {
      if ((c = getc(fp)) != '\n')         // \r\n : Windows
        ungetc(c, fp);                    // \r   : Mac (old)
      c = '\n'; break;
    }
    i++;
  }

  if (i_) *i_ = i;

  return c;
}

static inline int
readLine (FILE *fp, char *buf, int n, int *i_)
{
  int c = 0, i = 0;

  while (i < n-1 && (c = getc(fp)) != EOF) {
    if (c == '\n') break;                 // \n   : Unix, Linux, MacOSX
    if (c == '\r') {
      if ((c = getc(fp)) != '\n')         // \r\n : Windows
        ungetc(c, fp);                    // \r   : Mac (old)
      c = '\n'; break;
    }
    buf[i++] = c;
  }
  buf[i] = 0;

  if (isComment(buf)) buf[i=0] = 0;

  if (i_) *i_ += i;

  return c;
}

#endif
