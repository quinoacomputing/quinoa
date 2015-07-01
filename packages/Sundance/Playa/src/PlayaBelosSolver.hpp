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

#ifndef PLAYA_BELOSSOLVER_HPP
#define PLAYA_BELOSSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "PlayaBelosAdapter.hpp"
#include "Teuchos_Describable.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosTFQMRSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   *
   */
  class BelosSolver : public LinearSolverBase<double>,
                      public Handleable<LinearSolverBase<double> >,
                      public Printable,
                      public Describable
  {
  public:
    /** */
    BelosSolver(const Teuchos::ParameterList& params);

    /** */
    virtual ~BelosSolver(){;}

    /** Set the preconditioning operator */
    void setUserPrec(const PreconditionerFactory<double>& pf) {pf_=pf;}

    /** \name Printable interface */
    //@{
    /** Write to a stream  */
    void print(std::ostream& os) const 
    {
      os << description() << std::endl;
    }
    //@}
    
    /** \name Describable interface */
    //@{
    /** Write a brief description */
    std::string description() const {return "BelosSolver";}
    //@}

    

    /** */
    virtual SolverState<double> solve(const LinearOperator<double>& op,
                                      const Vector<double>& rhs,
                                      Vector<double>& soln) const ;

    /** \name Handleable interface */
    //@{
    /** Return a ref count pointer to a newly created object */
    virtual RCP<LinearSolverBase<double> > getRcp() 
    {return rcp(this);}
    //@}


  protected:

  private:
    
    /** */
    PreconditionerFactory<double> pf_;
    /** */
    mutable RCP<Belos::SolverManager<double,Anasazi::SimpleMV, LinearOperator<double> > > solver_ ;
    /** */
    mutable bool hasSolver_;
  };
  
}

#endif
