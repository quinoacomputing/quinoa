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

#ifndef PLAYA_ICCPRECONDITIONERFACTORY_HPP
#define PLAYA_ICCPRECONDITIONERFACTORY_HPP

#include "PlayaDefs.hpp"
#include "PlayaPreconditionerFactoryBase.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaICCFactorizableOp.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"

namespace Playa
{
  using namespace Teuchos;

  /**
   * 
   */
  template <class Scalar>
  class ICCPreconditionerFactory
    : public PreconditionerFactoryBase<Scalar>
  {
  public:
    /** Magnitude type */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    /** Construct with a parameter list */
    //const ParameterList& params inside
    ICCPreconditionerFactory():
    	fillLevels_(1),
        overlapFill_(0),
	dropTolerance_(1.0e-4),
        relaxationValue_(0.0),
        relativeThreshold_(1.0),
        absoluteThreshold_(0.0){;}
        
    ICCPreconditionerFactory(int fillLevels, 
			     int overlapFill, 
			     ScalarMag dropTolerance,
			     ScalarMag relaxationValue, 
			     ScalarMag relativeThreshold, 
			     ScalarMag absoluteThreshold)
    {
    	fillLevels_=fillLevels;
    	overlapFill_=overlapFill;
	dropTolerance_ = dropTolerance;
    	relaxationValue_=relaxationValue;
    	relativeThreshold_=relativeThreshold;
    	absoluteThreshold_=absoluteThreshold;
    }

    /** virtual dtor */
    virtual ~ICCPreconditionerFactory(){;}

    
    /** */
    virtual Preconditioner <Scalar>
    createPreconditioner(const LinearOperator<Scalar>& A) const 
    {
      /* In order for ICC factorization to work, the operator A must
       * implement the ICCFactorizableOp interface. We cast A's pointer
       * to a ICCFactorizableOp ptr. If the cast fails, throw a spoke. */
      
      const ICCFactorizableOp<Scalar>* fop 
        = dynamic_cast<const ICCFactorizableOp<Scalar>*>(A.ptr().get());

      TEUCHOS_TEST_FOR_EXCEPTION(fop==0, std::runtime_error,
                         "ICCPreconditionerFactory attempted to "
                         "create an ICC preconditioner for an operator type "
                         "that does not implement the ICCFactorizableOp "
                         "interface. The op is " << A.description());

      
      /* Now we can delegate the construction of the ICC factors to 
      * the factorizable op. */
      Preconditioner<Scalar> P;
      fop->getICCPreconditioner(fillLevels_,
				overlapFill_,
				dropTolerance_,
				relaxationValue_,
				relativeThreshold_,
				absoluteThreshold_,
				P);
      /* Return the preconditioner */
      return P;
    }

    /* Handleable boilerplate */
    GET_RCP(PreconditionerFactoryBase<Scalar>);
  private:

    int fillLevels_;
    int overlapFill_;
    ScalarMag dropTolerance_;
    ScalarMag relaxationValue_;
    ScalarMag relativeThreshold_;
    ScalarMag absoluteThreshold_;
  };


}

#endif
