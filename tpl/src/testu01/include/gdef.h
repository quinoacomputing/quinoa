 
/* gdef.h  for ANSI C */
#ifndef GDEF_H
#define GDEF_H

#include <gdefconf.h>
#include <limits.h>

#ifdef HAVE_LONG_LONG
#define USE_LONGLONG
#else
#undef USE_LONGLONG
#endif

#ifdef HAVE_GMP_H
#define USE_GMP
#else
#undef USE_GMP
#endif
 

#define FALSE 0
#define TRUE 1

typedef int  lebool;



#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#ifndef HAVE_UINT32_T
#if UINT_MAX >= 4294967295UL
   typedef unsigned int  uint32_t;
#else
   typedef unsigned long  uint32_t;
#endif
#endif

#ifndef HAVE_UINT8_T
   typedef unsigned char  uint8_t;
#endif



#define USE_LONGLONG



#ifdef USE_LONGLONG
   typedef long long  longlong;
   typedef unsigned long long  ulonglong;
#define PRIdLEAST64  "lld"
#define PRIuLEAST64  "llu"
#endif



#if ULONG_MAX == 4294967295UL
#define IS_ULONG32
#endif



#undef USE_ANSI_CLOCK



#define DIR_SEPARATOR "/"



#undef USE_GMP



#undef HAVE_MATHEMATICA


void gdef_GetHostName (char machine[], int n);



void gdef_WriteHostName (void);

#endif

