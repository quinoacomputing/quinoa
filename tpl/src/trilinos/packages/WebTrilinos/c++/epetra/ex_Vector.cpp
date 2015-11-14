// @HEADER
// ***********************************************************************
//
//               WebTrilinos: A Web Interface to Trilinos
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ***********************************************************************
// @HEADER

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

  // Total number of elements in the vector
  int NumElements = 10;

  // Construct a Map with NumElements and index base of 0
  Epetra_Map Map(NumElements, 0, Comm);

  // Create x as an Epetra_vector
  Epetra_Vector x(Map);
 
  // get the local size of the vector
  int MyLength = x.MyLength();

  /* First way to define the vector:   */
  /* use the [] operator on the object */

  for (int i = 0 ; i < MyLength ; ++i) x[i] = 1.0 * i;

  // need a double pointer because this works with multi-vectors
  double * pointer;
  
  x.ExtractView( &pointer );

  for (int i = 0; i < MyLength; ++i)
    cout << "on proc " << Comm.MyPID() << ", x["
	 << i << "] = " << pointer[i] << endl;
  
  // now modify the values
  for (int i = 0; i < MyLength; ++i)
    pointer[i] *= 10;

  // this affects the object x
  cout << x;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
} /* main */
