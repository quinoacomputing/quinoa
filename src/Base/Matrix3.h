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
//  externally callable functions defined in matrix3.cc
//

#include <stdio.h>


double m3det( double *m );
void m3inv( double *m, double *invm );
void m3mult( double *A, double *B, double *mult );
void m3exp( double *m, double *expm );
void m3trans( double *m, double *trans );
void m3cholesky( double *m );
void m3out( double *m, char *name, FILE *ofile );

