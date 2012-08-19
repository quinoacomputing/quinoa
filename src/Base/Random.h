// -----------------------------------------------------------------------------
// \file    src/Base/Macros.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Utilities for random numbers
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

void preprng_tables( int ngr, int nthreads, int restarted, int samenthreads,
                     #ifndef WALLFUNCTIONS
                     int nur, double **ru,
		     #endif
		     double **rg
		   );

void destroyrng_tables( int nthreads,
                        #ifndef WALLFUNCTIONS
                        double **ru,
			#endif
			double **rg
                      );

void regenrng_tables( int nthreads,
                      #ifndef WALLFUNCTIONS
                      double *ru,
		      #endif
		      double *rg );

void preprng_streams( int nthreads, VSLStreamStatePtr **stream
                      #ifndef WALLFUNCTIONS
                      , int restarted, int samenthreads
		      #endif
		    );

void destroyrng_streams( int nthreads, VSLStreamStatePtr **stream );

void saverng_streams( int nthreads
                      #ifndef WALLFUNCTIONS
                      , VSLStreamStatePtr *stream
                      #endif
                    );
