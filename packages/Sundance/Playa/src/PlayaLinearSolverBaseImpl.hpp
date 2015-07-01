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

#ifndef PLAYA_LINEARSOLVERBASEIMPL_HPP
#define PLAYA_LINEARSOLVERBASEIMPL_HPP

#include "PlayaLinearSolverBaseDecl.hpp"
#include "PlayaPreconditioner.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "Teuchos_ParameterList.hpp"


using namespace Teuchos;


namespace Playa
{

template <class Scalar> inline
const ParameterList& LinearSolverBase<Scalar>::parameters() const 
{return params_;}

template <class Scalar> inline
LinearSolverBase<Scalar>::LinearSolverBase(const ParameterList& params)
  : ObjectWithVerbosity(), params_(params) 
{
  if (this->parameters().isParameter(this->verbosityParam()))
  {
    this->setVerb(this->parameters().template get<int>(this->verbosityParam()));
  }
}

template <class Scalar> inline
ParameterList& LinearSolverBase<Scalar>::parameters() {return params_;}


template <class Scalar> inline
string LinearSolverBase<Scalar>::verbosityParam() const {return "Verbosity";}

template <class Scalar>
template <typename T> inline
void LinearSolverBase<Scalar>::setParameter(const ParameterList& params,
  T* dataPtr,
  const std::string& name)
{
  if (!params.isParameter(name)) return;

  TEUCHOS_TEST_FOR_EXCEPTION(!params.template isType<T>(name), std::runtime_error,
    "invalid type for parameter [" << name << "]"); 

  *dataPtr = params.template get<T>(name);
}

template <class Scalar> inline
void LinearSolverBase<Scalar>::setUserPrec(const PreconditionerFactory<Scalar>& pf)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "User-defined preconditioning not allowed for generic "
    "linear solver subtypes");
}

template <class Scalar> inline
void LinearSolverBase<Scalar>::setUserPrec(const LinearOperator<Scalar>& P,
  const LinearSolver<Scalar>& pSolver)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "User-defined preconditioning not allowed for generic "
    "linear solver subtypes");
}

}

#endif
