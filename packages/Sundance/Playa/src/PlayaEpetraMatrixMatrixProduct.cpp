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
#include "PlayaEpetraMatrixMatrixProduct.hpp"
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


LinearOperator<double> epetraLeftScale(
  const Vector<double>& d,
  const LinearOperator<double>& A)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);
  
  /* Make a deep copy of A */
  RCP<Epetra_CrsMatrix> mtxCopy = rcp(new Epetra_CrsMatrix(*A_crs));

  /* Extract the underlying Epetra vector. Type checking is done
   * internally, so we need no error check here. */
  const Epetra_Vector& epv = EpetraVector::getConcrete(d);
  
  /* Scale the copy */
  mtxCopy->LeftScale(epv);

  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(mtxCopy, A.domain(), A.range()));
  return rtn;
  
}

LinearOperator<double> epetraRightScale(
  const LinearOperator<double>& A,
  const Vector<double>& d)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);
  
  /* Make a deep copy of A */
  RCP<Epetra_CrsMatrix> mtxCopy = rcp(new Epetra_CrsMatrix(*A_crs));

  /* Extract the underlying Epetra vector. Type checking is done
   * internally, so we need no error check here. */
  const Epetra_Vector& epv = EpetraVector::getConcrete(d);
  
  /* Scale the copy */
  mtxCopy->RightScale(epv);

  /* Prepare an operator object for the scaled matrix */
  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(mtxCopy, A.domain(), A.range()));
  return rtn;
  
}


LinearOperator<double> epetraMatrixMatrixProduct(
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
  

  /* Get the row map from A. We will need this to build the target matrix C */
  const Epetra_Map* rowmap 
    = transA ? &(A_crs->DomainMap()) : &(A_crs->RowMap());

  /* make the target matrix */
  RCP<Epetra_CrsMatrix> C = rcp(new Epetra_CrsMatrix(Copy, *rowmap, 1));

  /* Carry out the multiplication */
  int ierr 
    = EpetraExt::MatrixMatrix::Multiply(*A_crs, transA, *B_crs, transB, *C);
  TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, RuntimeError,
    "EpetraExt Matrix-matrix multiply failed with error code ierr=" << ierr);

  /* Prepare an operator object for the scaled matrix */
  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(C, B.domain(), A.range()));
  return rtn;
  
}

}
