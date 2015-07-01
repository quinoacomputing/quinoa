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

#ifndef PLAYA_ITERATIVESOLVER_HPP
#define PLAYA_ITERATIVESOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   *
   */
  template <class Scalar>
  class IterativeSolver : public LinearSolverBase<Scalar>
  {
  public:
    /** */
    IterativeSolver(const ParameterList& params = ParameterList());

    /** */
    virtual ~IterativeSolver(){;}
    
    /** */
    int getMaxiters() const 
    {return this->parameters().template get<int>(maxitersParam());}

    /** */
    Scalar getTol() const 
    {return this->parameters().template get<double>(tolParam());}

    /** Change the convergence tolerance. */
    virtual void updateTolerance(const double& tol)
    {getParameter<double>(this->parameters(), tolParam())=tol;}

    /** */
    static std::string maxitersParam() {return "Max Iterations";}

    /** */
    static std::string tolParam() {return "Tolerance";}

    /** */
    static int defaultMaxiters() {return 500;}

    /** */
    static Scalar defaultTol() {return 1.0e-10;}
  };

  
  template <class Scalar> inline
  IterativeSolver<Scalar>::IterativeSolver(const ParameterList& params)
    : LinearSolverBase<Scalar>(params)
  {;}
  
}

#endif
