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

#include "PlayaIfpackICCOperator.hpp"
#include "Teuchos_Array.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaEpetraVector.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

namespace Playa
{

  using namespace Teuchos;

  IfpackICCOperator::IfpackICCOperator(const EpetraMatrix* A,
				       int fillLevels,
				       int overlapFill,
				       double dropTol,
				       double relaxationValue,
				       double relativeThreshold,
				       double absoluteThreshold)
    : LinearOpWithSpaces<double>(A->domain(), A->range()),
      precond_()
  {
    const Epetra_CrsMatrix* matrix = A->crsMatrix();
  
    Ifpack_ICT* precond = new Ifpack_ICT(matrix);
    precond_ = rcp(precond);

    Teuchos::ParameterList List;
  
    List.set("fact: ict level-of-fill", (double) fillLevels);
    List.set("fact: absolute threshold", absoluteThreshold);
    List.set("fact: relative threshold", relativeThreshold);
    List.set("fact: relax value", relaxationValue);
    List.set("fact: drop tolerance", dropTol);

    int ierr = precond->SetParameters(List);
    TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
			       "IfpackOperator ctor: "
			       "precond->SetParameters() failed with ierr="
			       << ierr);
  
    ierr = precond->Compute();

    TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
			       "IfpackOperator ctor: "
			       "precond->Compute() failed with ierr="
			       << ierr);
  }
  
  void IfpackICCOperator::apply(Teuchos::ETransp transApplyType,
				const Vector<double>& in,
				Vector<double> out) const
  {
    /* grab the epetra vector objects underlying the input and output vectors */
    const Epetra_Vector& epIn = EpetraVector::getConcrete(in);
    Epetra_Vector& epOut = EpetraVector::getConcrete(out);

    /* ifpack's solve is logically const but declared non-const because
     *  internal data changes. So, do a const_cast. */
    Ifpack_ICT* p = const_cast<Ifpack_ICT*>(precond_.get());

    int ierr;

    /* do the solve (or transpose solve) */
    if (transApplyType==NO_TRANS)
      {
	ierr = p->H().Solve(false,false,false, epIn, epOut);
      }
    else
      {
	ierr = p->H().Solve(false,true,false, epIn, epOut);
      }

    TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
			       "Playa::IfpackICCOperator apply: "
			       "p->H().Solve() (transpose state="
			       << transApplyType << ") failed with ierr="
			       << ierr);
  }
  
  void IfpackICCOperator::print(std::ostream& os) const 
  {
    precond_->Print(os);
  }

  string IfpackICCOperator::description() const 
  {
    std::string rtn = "ICC Preconditioner";
    return rtn;
  }

}
