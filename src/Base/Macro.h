//******************************************************************************
/*!
  \file      src/Base/Macro.h
  \author    J. Bakosi
  \date      Thu Aug 29 14:54:18 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Macro definitions
  \details   Macro definitions
*/
//******************************************************************************
#ifndef Macro_h
#define Macro_h

#include <sys/time.h>
#include <ctime>
#include <cstdlib>

namespace quinoa {

#define SWAP(a,b,tmp) {tmp=a; a=b; b=tmp;}
#define IGNORE(expr) (static_cast<void>(expr))

// macros for fine-grained profiling
#define STARTTIME \
struct timeval START_TIME, END_TIME; \
int total_usecs; \
gettimeofday(&START_TIME, (struct timezone*)0);

#define ENDTIME \
gettimeofday(&END_TIME, (struct timezone*)0); \
total_usecs = (END_TIME.tv_sec-START_TIME.tv_sec) * 1000000 + (END_TIME.tv_usec-START_TIME.tv_usec); \
printf("Total time was %d uSec.\n", total_usecs);

} // namespace quinoa

#endif // Macro_h
