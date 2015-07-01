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

#ifndef PLAYA_BICGSTABSOLVER_DECL_HPP
#define PLAYA_BICGSTABSOLVER_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaKrylovSolver.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "PlayaOut.hpp"
#include "Teuchos_Describable.hpp"


namespace Playa
{
using namespace Teuchos;

/**
 *
 */
template <class Scalar>
class BICGSTABSolver : public KrylovSolver<Scalar>,
                       public Playa::Handleable<LinearSolverBase<Scalar> >,
                       public Printable,
                       public Describable
{
public:
  /** */
  BICGSTABSolver(const ParameterList& params = ParameterList());

  /** */
  BICGSTABSolver(const ParameterList& params,
    const PreconditionerFactory<Scalar>& precond);

  /** */
  virtual ~BICGSTABSolver(){;}

  /** \name Printable interface */
  //@{
  /** Write to a stream  */
  void print(std::ostream& os) const ;
  //@}
    
  /** \name Describable interface */
  //@{
  /** Write a brief description */
  std::string description() const {return "BICGSTABSolver";}
  //@}

  /** \name Handleable interface */
  //@{
  /** Return a ref count pointer to a newly created object */
  virtual RCP<LinearSolverBase<Scalar> > getRcp() 
    {return rcp(this);}
  //@}
    
protected:

  /** */
  virtual SolverState<Scalar> solveUnprec(const LinearOperator<Scalar>& op,
    const Vector<Scalar>& rhs,
    Vector<Scalar>& soln) const ;

    
};


}

#endif
