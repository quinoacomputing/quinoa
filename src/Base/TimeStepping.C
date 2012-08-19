// -----------------------------------------------------------------------------
// \file    src/Base/TimeStepping.C
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Time stepping
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

#include <sys/time.h>
#include <ctime>
#include "Const.h"
#include "TimeStepping.h"

void calctimes(struct timeval *start_time, int it, int it0, double t, double dt,
	       long int *hrs2beg, long int *mins2beg, long int *secs2beg,
	       long int *hrs2end, long int *mins2end, long int *secs2end )
// -----------------------------------------------------------------------------
// Routine: calctimes - Calculate elapsed and estimated time remaining
//                      in hours:mins:seconds
// Author : J. Bakosi
// -----------------------------------------------------------------------------
{
  struct timeval cur_time;
  long int secs_elapsed;

  gettimeofday( &cur_time, (struct timezone*)0 );
  secs_elapsed = ((cur_time.tv_sec-start_time->tv_sec) * 1000000 +
                  (cur_time.tv_usec-start_time->tv_usec))/1000000;

  // calculate elapsed time...
  *secs2beg = secs_elapsed;
  *mins2beg = (*secs2beg)/60;
  *hrs2beg = (*mins2beg)/60;
  if ( *secs2beg >= 60 ) *secs2beg %= 60;
  if ( *secs2beg >= 60 ) *secs2beg %= 60;
  if ( *mins2beg >= 60 ) *mins2beg %= 60;

  // estimate time remaining...
  if (it-it0) *secs2end = (long int)(secs_elapsed*(MAXTIME-t)/(dt*(it-it0)));
  else *secs2end = 0;
  *mins2end = (*secs2end)/60;
  *hrs2end = (*mins2end)/60;
  if ( *secs2end >= 60 ) *secs2end %= 60;
  if ( *secs2end >= 60 ) *secs2end %= 60;
  if ( *mins2end >= 60 ) *mins2end %= 60;
}
