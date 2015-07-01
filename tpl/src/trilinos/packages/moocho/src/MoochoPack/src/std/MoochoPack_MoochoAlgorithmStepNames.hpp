// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RSQP_ALGORITHM_STEP_NAMES_H
#define RSQP_ALGORITHM_STEP_NAMES_H

#include <string>

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** @name Names for MOOCHOsteps */
//@{

const std::string EvalNewPoint_name                 = "EvalNewPoint";
const std::string ReducedGradient_name              = "ReducedGradient";
const std::string ReducedHessian_name               = "ReducedHessian";
const std::string QuasiNormalStep_name              = "QuasiNormalStep";
const std::string TangentialStep_name               = "TangentialStep";
const std::string SearchDirec_name                  = "SearchDirec";
const std::string LineSearch_name                   = "LineSearch";
const std::string CheckConvergence_name             = "CheckConvergence";

const std::string CalcLambdaIndep_name              = "CalcLambdaIndep";
const std::string CalcReducedGradLagrangian_name    = "CalcReducedGradLagrangian";
const std::string CheckSkipBFGSUpdate_name          = "CheckSkipBFGSUpdate";
const std::string CalcDFromYPYZPZ_name				= "CalcDFromYPYZPZ";

//@}
}	// end namespace MoochoPack

#endif // RSQP_ALGORITHM_STEP_NAMES_H
