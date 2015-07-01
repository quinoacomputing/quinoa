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

#ifndef EPETRA_PLAYA_OPERATOR_HPP
#define EPETRA_PLAYA_OPERATOR_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "Epetra_Operator.h"

namespace Epetra
{
using namespace Teuchos;
using namespace Playa;
  
  

/** */
class Epetra_PlayaOperator : public Epetra_Operator
{
public:
  /** */
  Epetra_PlayaOperator(const LinearOperator<double>& A,
    const LinearSolver<double>& solver=LinearSolver<double>());
    
  /** */
  int SetUseTranspose(bool useTrans) {useTranspose_ = useTrans; return 0;}

  /** */
  int Apply(const Epetra_MultiVector& in, Epetra_MultiVector& out) const ;

  /** */
  int ApplyInverse(const Epetra_MultiVector& in, Epetra_MultiVector& out) const ;

  /** */
  double NormInf() const ;

  /** */
  const char* Label() const ;

  /** */
  bool UseTranspose() const {return useTranspose_;}

  /** */
  bool HasNormInf() const {return false;}

  /** */
  const Epetra_Comm& Comm() const {return *comm_;}

  /** */
  const Epetra_Map& OperatorDomainMap() const {return *domain_;}

  /** */
  const Epetra_Map& OperatorRangeMap() const {return *range_;}

    

private:
  LinearOperator<double> A_;
  LinearSolver<double> solver_;
  bool useTranspose_;
  RCP<Epetra_Comm> comm_;
  RCP<const Epetra_Map> domain_;
  RCP<const Epetra_Map> range_;
  bool isNativeEpetra_;
  bool isCompoundEpetra_;
  std::string label_;
};
}

#endif
