//  ------------------------------------------------------------------------------------------------------------
//
//  Copyright 2007 Jozsef Bakosi
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  ------------------------------------------------------------------------------------------------------------
//
//  A diagnostic function for estimating times needed for run.
//  For more info see main.cc
//

#include <sys/time.h>
#include <time.h>
#include "const.h"
#include "timestepping.h"







void calctimes( struct timeval *start_time, int it, int it0, double t, double dt,
	       long int *hrs2beg, long int *mins2beg, long int *secs2beg,
	       long int *hrs2end, long int *mins2end, long int *secs2end )
//
// calculates elapsed time and estimated time remaining in hours:mins:seconds
//
{
  struct timeval cur_time;
  long int secs_elapsed;


  gettimeofday( &cur_time, (struct timezone*)0 );
  secs_elapsed = ((cur_time.tv_sec-start_time->tv_sec) * 1000000 + (cur_time.tv_usec-start_time->tv_usec))/1000000;

  // calculate elapsed time...
  *secs2beg = secs_elapsed;
  *mins2beg = (*secs2beg)/60;
  *hrs2beg = (*mins2beg)/60;
  if ( *secs2beg >= 60 ) *secs2beg %= 60;
  if ( *secs2beg >= 60 ) *secs2beg %= 60;
  if ( *mins2beg >= 60 ) *mins2beg %= 60;

  // estimate time remaining...
  if (it-it0) *secs2end = (long int)(secs_elapsed*(MAXTIME-t)/(dt*(it-it0))); else *secs2end = 0;
  *mins2end = (*secs2end)/60;
  *hrs2end = (*mins2end)/60;
  if ( *secs2end >= 60 ) *secs2end %= 60;
  if ( *secs2end >= 60 ) *secs2end %= 60;
  if ( *mins2end >= 60 ) *mins2end %= 60;
}
