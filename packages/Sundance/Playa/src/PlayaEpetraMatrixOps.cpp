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
#include "PlayaEpetraMatrixFactory.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"  // changed from Impl
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"  // changed from Impl

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;




namespace Playa
{

Vector<double> getEpetraDiagonal(const LinearOperator<double>& A)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);

  VectorSpace<double> rowSpace = A.domain();
  Vector<double> rtn = rowSpace.createMember();

  Epetra_Vector* xPtr = EpetraVector::getConcretePtr(rtn);
  A_crs->ExtractDiagonalCopy(*xPtr);

  return rtn;
}


LinearOperator<double> makeEpetraDiagonalMatrix(const Vector<double>& d)
{
  VectorSpace<double> space = d.space();
  RCP<const EpetraVectorSpace> eps 
    = rcp_dynamic_cast<const EpetraVectorSpace>(space.ptr());

  EpetraMatrixFactory mf(eps, eps);

  int nLocal = space.numLocalElements();
  int offset = space.baseGlobalNaturalIndex();
  for (int i=0; i<nLocal; i++)
  {
    int rowIndex = offset + i;
    mf.initializeNonzerosInRow(rowIndex, 1, &rowIndex);
  }

  mf.finalize();
  LinearOperator<double> rtn = mf.createMatrix();

  RCP<EpetraMatrix> epm = rcp_dynamic_cast<EpetraMatrix>(rtn.ptr());
  epm->zero();
  
  for (int i=0; i<nLocal; i++)
  {
    int rowIndex = offset + i;
    double val = d[i];
    epm->addToRow(rowIndex, 1, &rowIndex, &val);
  }
  
  return rtn;
}

}
