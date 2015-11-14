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
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList GaleriList;

  // The problem is defined on a 2D grid, global size is nx * nx.
  int nx = 30; 
  GaleriList.set("n", nx * nx);
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx);
  Epetra_Map* Map = Galeri::CreateMap("Linear", Comm, GaleriList);
  Epetra_RowMatrix* A = Galeri::CreateCrsMatrix("Laplace2D", Map, GaleriList);

  // defines LHS and RHS
  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorDomainMap());

  // solution is constant
  LHS.PutScalar(1.0);
  // now build corresponding RHS
  A->Apply(LHS,RHS);

  // now randomize the solution
  RHS.Random();

  // need an Epetra_LinearProblem to define AztecOO solver
  Epetra_LinearProblem Problem(A,&LHS,&RHS);

  // now we can allocate the AztecOO solver
  AztecOO Solver(Problem);

  // specify options
  Solver.SetAztecOption(AZ_solver,AZ_gmres);
  Solver.SetAztecOption(AZ_output,32);
  Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut); 
  Solver.SetAztecOption(AZ_graph_fill, 3);
  Solver.SetAztecOption(AZ_overlap, 0);

  // .. and here we solve
  Solver.Iterate(1550,1e-8);

  // delete memory
  delete Map;
  delete A;
  
#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}
