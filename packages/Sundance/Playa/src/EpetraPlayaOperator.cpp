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
#include "PlayaEpetraVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"  // changed from Impl
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"  // changed from Impl
#include "Teuchos_Array.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaIfpackILUOperator.hpp"
#include "PlayaGenericLeftPreconditioner.hpp"
#include "PlayaGenericRightPreconditioner.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_getConst.hpp"
#include "EpetraPlayaOperator.hpp"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;

namespace Epetra
{

Epetra_PlayaOperator::Epetra_PlayaOperator(const LinearOperator<double>& A,
  const LinearSolver<double>& solver)
  : A_(A), solver_(solver), useTranspose_(false), comm_(), domain_(), range_(),
    isNativeEpetra_(false), isCompoundEpetra_(false), label_(A.description())
{
  const EpetraMatrix* em = dynamic_cast<const EpetraMatrix*>(A.ptr().get());
  const EpetraVectorSpace* ed = dynamic_cast<const EpetraVectorSpace*>(A.domain().ptr().get());
  const EpetraVectorSpace* er = dynamic_cast<const EpetraVectorSpace*>(A.range().ptr().get());

  if (em)
  {
    isNativeEpetra_ = true;
    const Epetra_CrsMatrix* crs = em->crsMatrix();
    domain_ = rcp(new Epetra_Map(crs->OperatorDomainMap()));
    range_ = rcp(new Epetra_Map(crs->OperatorRangeMap()));
    useTranspose_ = crs->UseTranspose();
    comm_ = rcp(crs->Comm().Clone());
  }
  else if (er != 0 && ed != 0)
  {
    domain_ = ed->epetraMap();
    range_ = er->epetraMap();
    comm_ = rcp(domain_->Comm().Clone());
    isCompoundEpetra_ = true;
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }
}


int Epetra_PlayaOperator::Apply(const Epetra_MultiVector& in, Epetra_MultiVector& out) const
{
  if (isNativeEpetra_)
  {
    const EpetraMatrix* em = dynamic_cast<const EpetraMatrix*>(A_.ptr().get());
    return em->crsMatrix()->Multiply(useTranspose_, in, out);
  }
  else if (isCompoundEpetra_)
  {
    const Epetra_Vector* cevIn = dynamic_cast<const Epetra_Vector*>(&in);
    Epetra_Vector* evIn = const_cast<Epetra_Vector*>(cevIn);
    Epetra_Vector* evOut = dynamic_cast<Epetra_Vector*>(&out);
    TEUCHOS_TEST_FOR_EXCEPTION(evIn==0, std::runtime_error, "Epetra_PlayaOperator::Apply "
      "cannot deal with multivectors");
    TEUCHOS_TEST_FOR_EXCEPTION(evOut==0, std::runtime_error, "Epetra_PlayaOperator::Apply "
      "cannot deal with multivectors");


    RCP<VectorBase<double> > vpIn 
      = rcp(new EpetraVector(A_.domain(), rcp(evIn, false)));
    RCP<VectorBase<double> > vpOut 
      = rcp(new EpetraVector(A_.range(), rcp(evOut, false)));
    Vector<double> vIn = vpIn;
    Vector<double> vOut = vpOut;

    A_.apply(vIn, vOut);
    out = EpetraVector::getConcrete(vOut);
    return 0;
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPT(true);
    return -1; // -Wall
  }
}

int Epetra_PlayaOperator::ApplyInverse(const Epetra_MultiVector& in, Epetra_MultiVector& out) const
{
  
  TEUCHOS_TEST_FOR_EXCEPTION(solver_.ptr().get()==0, std::runtime_error,
    "no solver provided for Epetra_PlayaOperator::ApplyInverse");
  TEUCHOS_TEST_FOR_EXCEPTION(!isNativeEpetra_ && !isCompoundEpetra_, std::runtime_error,
    "Epetra_PlayaOperator::ApplyInverse expects either "
    "a native epetra operator or a compound operator with "
    "Epetra domain and range spaces");
  const Epetra_Vector* cevIn = dynamic_cast<const Epetra_Vector*>(&in);
  Epetra_Vector* evIn = const_cast<Epetra_Vector*>(cevIn);
  Epetra_Vector* evOut = dynamic_cast<Epetra_Vector*>(&out);

  TEUCHOS_TEST_FOR_EXCEPTION(evIn==0, std::runtime_error, "Epetra_PlayaOperator::Apply "
    "cannot deal with multivectors");
  TEUCHOS_TEST_FOR_EXCEPTION(evOut==0, std::runtime_error, "Epetra_PlayaOperator::Apply "
    "cannot deal with multivectors");

  RCP<VectorBase<double> > vpIn 
    = rcp(new EpetraVector(A_.domain(),
        rcp(evIn, false)));
  RCP<VectorBase<double> > vpOut 
    = rcp(new EpetraVector(A_.range(),
        rcp(evOut, false)));
  Vector<double> vIn = vpIn;
  Vector<double> vOut = vpOut;
  
  SolverState<double> state = solver_.solve(A_, vIn, vOut);

  if (state.finalState() == SolveCrashed) return -1;
  else if (state.finalState() == SolveFailedToConverge) return -2;
  else out = EpetraVector::getConcrete(vOut);

  return 0;
}




double Epetra_PlayaOperator::NormInf() const 
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return -1; // -Wall
}

const char* Epetra_PlayaOperator::Label() const 
{
  return label_.c_str();
}

}
