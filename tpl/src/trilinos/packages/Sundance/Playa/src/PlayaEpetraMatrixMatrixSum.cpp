/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#include "PlayaEpetraMatrix.hpp"
#include "PlayaEpetraMatrixMatrixSum.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaExceptions.hpp"
#include "EpetraExt_MatrixMatrix.h"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif



namespace Playa
{
using namespace Teuchos;


LinearOperator<double> epetraMatrixMatrixSum(
  const LinearOperator<double>& A,
  const LinearOperator<double>& B)
{
  /* Extract the underlying Epetra matrix for A. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);

  /* Extract the underlying Epetra matrix for A. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> B_crs = EpetraMatrix::getConcretePtr(B);
  
  bool transA = false;
  bool transB = false;

  TEUCHOS_TEST_FOR_EXCEPTION(A.range() != B.range(), RuntimeError,
    "incompatible operand ranges in epetraMatrixMatrixSum()"
    << std::endl << "A.range()=" << A.range()
    << std::endl << "B.range()=" << B.range()
    );
  

  TEUCHOS_TEST_FOR_EXCEPTION(A.domain() != B.domain(), RuntimeError,
    "incompatible operand domains in epetraMatrixMatrixSum()"
    << std::endl << "A.domain()=" << A.domain()
    << std::endl << "B.domain()=" << B.domain()
    );
  

  /* Get the row map from A. We will need this to build the target matrix C */
  const Epetra_Map* rowmap 
    = transA ? &(A_crs->DomainMap()) : &(A_crs->RowMap());

  /* make the target matrix */
  RCP<Epetra_CrsMatrix> C = rcp(new Epetra_CrsMatrix(Copy, *rowmap, 1));
  Epetra_CrsMatrix* CPtr = C.get();

  /* Carry out the multiplication */
  int ierr 
    = EpetraExt::MatrixMatrix::Add(
      *A_crs, transA, 1.0, 
      *B_crs, transB, 1.0, CPtr);
  TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, RuntimeError,
    "EpetraExt Matrix-matrix add failed with error code ierr=" << ierr);

  /* Need to call FillComplete on the result */
  C->FillComplete();

  /* Prepare an operator object for the added matrix */
  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(C, B.domain(), A.range()));
  return rtn;
  
}

}
