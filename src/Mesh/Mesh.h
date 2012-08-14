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
//  externally callable functions defined in mesh.cc
//

void prepmsh( int *npoin, int *nbpoin, int *nelem, double **coord, int **bpg, int **binpoel,
	      int **inpoel, int **esup1, int **esup2, int **psup1, int **psup2, int **bpsup1,
	      int **bpsup2, int **esupel1, int **esupel2, int **esuel, int**bptags, int **betags,
	      double **Ae, double **dNx, double **dNy, double **dete, double **sqrtAp, double *minsqrtAp );
