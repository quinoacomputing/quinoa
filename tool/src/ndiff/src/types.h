#ifndef TYPES_H
#define TYPES_H

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
     provides portable types and limits
 
 o---------------------------------------------------------------------o
*/

#include <limits.h>

typedef unsigned char          bool;
typedef unsigned int           uint;
typedef unsigned long int      ulong;
typedef unsigned long long int ullong;

enum { false, true };

#endif
