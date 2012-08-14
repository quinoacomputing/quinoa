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


// structure prototype for sparse matrices
struct sparsemat {
      int size;			// size: 1,		matrix of (dof x size) x (dof x size)
      int rsize;		// size: 1,		real size of matrix (rsize = dof x size)
      int dof;			// size: 1,		number of unknowns/node (degree of freedom)
      int nnz;			// size: 1,		total number of nonzeros
      int *rnz;			// size: size,		number of nonzeros of each row
      int *ia;			// size: size*dof+1,	row pointers
      int *ja;			// size: nnz,		column indices
      double *a;		// size: nnz,		nonzero values
      double *r, *p, *z, *q,	// sizes: rsize,	auxiliary vectors for conjugate gradients
             *d, *u;
      int nlev;			// size: 1,		number of levels for level scheduling
      int *lev1, *lev2;		// sizes: rsize,nlev+1	linked lists for level scheduling
};


void create_sparsemat( sparsemat *m, int DOF, int size, int *psup1, int *psup2 );
void destroy_sparsemat( sparsemat *m );
void addmr( sparsemat *M, int row, int column, int i, double value );
void addma( sparsemat *M, int row, int column, double value );
void insmr( sparsemat *M, int row, int column, int i, double value );
void insma( sparsemat *M, int row, int column, double value );
double getmr( sparsemat *M, int row, int column, int i );
double getma( sparsemat *M, int row, int column );
void dpzero( double *ptr, int size, int nthreads );
void ipzero( int *ptr, int size, int nthreads );
void printmat_as_stored( sparsemat *M );
void printmat_as_structure( sparsemat *M );
void printmat_as_matrix( sparsemat *M );
void printmat_as_matlab( sparsemat *M );
void printmat2file_as_stored( sparsemat *M, char *filename, double *rhs );
