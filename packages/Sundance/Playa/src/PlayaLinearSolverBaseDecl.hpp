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

#ifndef PLAYA_LINEARSOLVERBASEDECL_HPP
#define PLAYA_LINEARSOLVERBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaSolverState.hpp"
#include "PlayaObjectWithVerbosity.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  using namespace Teuchos;
  template <class Scalar>
  class LinearOperator;

  template <class Scalar>
  class Preconditioner;

  template <class Scalar>
  class PreconditionerFactory;

  template <class Scalar>
  class Vector;
  

  /** */
  template <class Scalar>
  class LinearSolverBase : public ObjectWithVerbosity
  {
  public:
    /** */
    LinearSolverBase(const ParameterList& params);

    /** */
    virtual ~LinearSolverBase(){;}

    /** */
    virtual SolverState<Scalar> solve(const LinearOperator<Scalar>& op,
                                      const Vector<Scalar>& rhs,
                                      Vector<Scalar>& soln) const = 0;

    /** Change the convergence tolerance. Default does nothing. */
    virtual void updateTolerance(const double& tol) {;}

    /** Set a user-defined preconditioning operator. Default is an error. */
    virtual void setUserPrec(const PreconditionerFactory<Scalar>& pf);

    /** Set a user-defined preconditioning operator. Default is an error. */
    virtual void setUserPrec(const LinearOperator<Scalar>& P,
      const LinearSolver<Scalar>& pSolver);

    /** */
    const ParameterList& parameters() const ;

    /** */
    ParameterList& parameters();

    /** */
    std::string verbosityParam() const ;

    /** */
    template <typename T>
    static void setParameter(const ParameterList& params,
                             T* valuePtr, 
                             const std::string& paramName);
  private:
    ParameterList params_;
  };
}

#endif
