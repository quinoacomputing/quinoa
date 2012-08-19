// -----------------------------------------------------------------------------
// \file    src/Base/TimeStepping.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Time stepping
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

void calctimes( struct timeval *start_time, int it, int it0, double t, double dt,
	       long int *hrs2beg, long int *mins2beg, long int *secs2beg,
	       long int *hrs2end, long int *mins2end, long int *secs2end );
