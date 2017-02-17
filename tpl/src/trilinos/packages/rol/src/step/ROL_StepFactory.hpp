// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_STEPFACTORY_H
#define ROL_STEPFACTORY_H

#include "ROL_Types.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "ROL_Step.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_AugmentedLagrangianStep.hpp"
#include "ROL_MoreauYosidaPenaltyStep.hpp"
#include "ROL_BundleStep.hpp"
#include "ROL_InteriorPointStep.hpp"

namespace ROL {

  template<class Real>
  class StepFactory {
    public:
    ~StepFactory(void){}

    Teuchos::RCP<Step<Real> > getStep(const std::string &type,
                                      Teuchos::ParameterList &parlist) const {
      EStep els = StringToEStep(type);
      switch(els) {
        case STEP_AUGMENTEDLAGRANGIAN: return Teuchos::rcp( new AugmentedLagrangianStep<Real>(parlist) );
        case STEP_BUNDLE:              return Teuchos::rcp( new BundleStep<Real>(parlist) );
        case STEP_COMPOSITESTEP:       return Teuchos::rcp( new CompositeStep<Real>(parlist) );
        case STEP_LINESEARCH:          return Teuchos::rcp( new LineSearchStep<Real>(parlist) );
        case STEP_MOREAUYOSIDAPENALTY: return Teuchos::rcp( new MoreauYosidaPenaltyStep<Real>(parlist) );
        case STEP_PRIMALDUALACTIVESET: return Teuchos::rcp( new PrimalDualActiveSetStep<Real>(parlist) );
        case STEP_TRUSTREGION:         return Teuchos::rcp( new TrustRegionStep<Real>(parlist) );
        case STEP_INTERIORPOINT:       return Teuchos::rcp( new InteriorPointStep<Real>(parlist) ); 
        default:                       return Teuchos::null;
      }
    }
  };

}

#endif
