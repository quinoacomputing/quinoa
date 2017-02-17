/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef CALL_EPETRA_H
#define CALL_EPETRA_H

/* This is a true C header file to be used to call Epetra functions from 
   the Euclid solver. They all pass the matrix as a void-pointer. */

#ifdef __cplusplus
extern "C"
{
#endif

  int ExtractIndicesView (void *A, int GlobalRow, int *NumEntries,
			  int **Indices);

  int ExtractValuesView (void *A, int GlobalRow, int *NumEntries,
			 double **Values);

  int MinMaxMyGID (void *A, bool Row, bool min);

  int NumGlobalRowCol (void *A, bool Row);

  int NumMyRowEntries (void *A, int Row, int *numEntries);

#ifdef __cplusplus
}
#endif


#endif
