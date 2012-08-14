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
//  externally callable functions defined in random.cc
//


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
