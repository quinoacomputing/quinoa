// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Notes:
//   This file deserves a Teuchos Parameter List parameter.
//   However, that would require making the triutils package dependent on Teuchos


#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "Trilinos_Util_CountTriples.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

int Trilinos_Util_ReadTriples2Epetra( char *data_file,
				      bool symmetric, 
				      const Epetra_Comm  &comm, 
				      Epetra_Map *& map, 
				      Epetra_CrsMatrix *& A, 
				      Epetra_Vector *& x, 
				      Epetra_Vector *& b,
				      Epetra_Vector *&xexact,
				      bool NonUniformMap=false,
				      bool TimDavisHeader=false,
				      bool ZeroBased=false ) {
  FILE *in_file ;
  int N_rows = 0, nnz = 0; 
  
  const int BUFSIZE = 800 ; 
  char buffer[BUFSIZE] ; 
  std::vector<int> non_zeros;   // Number of non-zeros in each row
  Trilinos_Util_CountTriples( data_file, symmetric, non_zeros, N_rows, nnz, comm, TimDavisHeader, ZeroBased ) ;

#if 0
  cout << " Trilinos_Util_ReadTriples2Epetra.cpp::" << __LINE__ << "  N_rows = " << N_rows << endl << "non_zeros = " ; 
  for (int i= 0; i<N_rows; i++ )  cout << non_zeros[i] ; 
  cout << endl ; 
  comm.Barrier();
#endif
 
  int NumNonZeroRows=0;
  std::vector<int> MapIndices;
  if(comm.MyPID() == 0)  { 
    if ( NonUniformMap ) {
      for (int i= 0; i<N_rows; i++ )  if( non_zeros[i] > 0 ) NumNonZeroRows++;
      MapIndices.resize(NumNonZeroRows);
      int IndexNumber =0;
      for (int i= 0; i<N_rows; i++ ) if( non_zeros[i] > 0 ) MapIndices[IndexNumber++]=i;
    }
  }
  //
  //  Copy MapIndices[] to all processes
  //
  if ( NonUniformMap ) {
    //    int value[1];
    //    values[0] = 
    comm.Broadcast( &NumNonZeroRows, 1, 0 ) ; 
    MapIndices.resize(NumNonZeroRows);
    comm.Broadcast( &MapIndices[0], NumNonZeroRows, 0 ) ; 
  }

  std::vector<int> ptrs(N_rows+1) ; // Pointers into inds and vals for the start of each row
  std::vector<int> inds(nnz);     //  Column Indices
  std::vector<double> vals(nnz);  //  Matrix values

  if(comm.MyPID() == 0)  { 
    //  ptrs, inds and vals together constitute a compressed row storage of the matrix.


    in_file = fopen( data_file, "r");
    assert (in_file != NULL) ;  // Checked in T_U_CountTriples() 

    ptrs[0] = 0 ; 
    for( int i=0; i< N_rows; i++ ) { 
      ptrs[i+1] = ptrs[i] + non_zeros[i]; 
    }

    std::vector<int> iptrs = ptrs ; //  Current pointers into inds and vals for each row

    if ( TimDavisHeader ) fgets( buffer, BUFSIZE, in_file ); // Throw away the Tim Davis Header Line 
    while ( fgets( buffer, BUFSIZE, in_file ) ) { 
      int i, j; 
      double val ; 
      sscanf( buffer, "%d %d %lg", &i, &j, &val ) ; 
      const int i_index = ( ZeroBased?i:i-1 );
      const int j_index = ( ZeroBased?j:j-1 );

      int iptr = iptrs[i_index] ; 
      iptrs[i_index]++ ;
      vals[iptr] = val ; 
      inds[iptr] = j_index ; 
      //
      //  If this is a symmetric matrix, we need to enter the entry 
      //  for the other triangular half
      //
      if (symmetric && i != j ) {
	iptr = iptrs[j_index] ; 
	iptrs[j_index]++;
	vals[iptr] = val ; 
	inds[iptr] = i_index ; 
      }
    } 
    fclose(in_file);
  }

  int nlocal = 0;
  if ( NonUniformMap ) {
    if (comm.MyPID()==0) nlocal = NumNonZeroRows;
    map = new Epetra_Map(NumNonZeroRows, nlocal, &MapIndices[0], 0, comm); // Create map with all elements on PE 0
  } else { 
    if (comm.MyPID()==0) nlocal = N_rows;
    map = new Epetra_Map(N_rows, nlocal, 0, comm); // Create map with all elements on PE 0
  }
  
  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0

  if (comm.MyPID()==0)
    for (int i=0; i<N_rows; i++) {
      if ( ptrs[i+1]>ptrs[i]) A->InsertGlobalValues(i, ptrs[i+1]-ptrs[i], &vals[ptrs[i]], &inds[ptrs[i]]);
    }
  A->FillComplete();

  std::vector<double> hbx(N_rows);
  std::vector<double> hbb(N_rows);
  std::vector<double> hbxexact(N_rows);

  
  x = new Epetra_Vector(Copy, *map, &hbx[0]);
  b = new Epetra_Vector(Copy, *map, &hbb[0]);
  xexact = new Epetra_Vector(Copy, *map, &hbxexact[0]);

  EPETRA_CHK_ERR( x->PutScalar( 0.0 ) );
  EPETRA_CHK_ERR( xexact->Random( ) ) ; 
  EPETRA_CHK_ERR( A->Multiply( false, *xexact, *b ) ); 

  return 0;
}
