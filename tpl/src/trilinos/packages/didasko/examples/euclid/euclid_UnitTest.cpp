// @HEADER
// ***********************************************************************
// 
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
// 
// ***********************************************************************
// @HEADER

#include "Ifpack.h"
#include "Ifpack_Euclid.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Epetra_MultiVector.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "euclid_Helpers.hpp"
#include "Teuchos_Array.hpp"
#include <string>
#include <stdio.h>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;

const double tol = 1E-7;
const int numVec = 1;

TEUCHOS_UNIT_TEST( Ifpack_Hypre, Euclid){
  RCP<Epetra_CrsMatrix> Matrix = rcp(newCrsMatrix(27));
  //cout << endl << *Matrix << endl;
  TEST_EQUALITY(Matrix->RowMap().LinearMap(), true);
  Ifpack_Euclid preconditioner(Matrix.get());
  TEST_EQUALITY(preconditioner.Initialize(),0);
  TEST_EQUALITY(preconditioner.SetParameter("setmem", 0),0);
  TEST_EQUALITY(preconditioner.SetParameter("setstats", 0), 0);
  TEST_EQUALITY(preconditioner.SetParameter("setlevel", 2), 0);
  TEST_EQUALITY(preconditioner.Compute(),0);
  cout << endl << preconditioner << endl;
  Epetra_MultiVector KnownX(Matrix->DomainMap(), numVec);
  KnownX.Random();

  Epetra_MultiVector B(Matrix->RangeMap(), numVec);
  TEST_EQUALITY(Matrix->Apply(KnownX, B), 0);

  Epetra_MultiVector X(Matrix->DomainMap(), numVec);

  AztecOO Solver(Matrix.get(), &X, &B);
  TEST_EQUALITY(Solver.SetPrecOperator(&preconditioner),0);
  Solver.Iterate(1000, 1E-9);
  //TEST_EQUALITY(preconditioner.ApplyInverse(B,X),0);
  TEST_EQUALITY(EquivalentVectors(X, KnownX, tol), true);
  cout << endl << preconditioner << endl;
}

