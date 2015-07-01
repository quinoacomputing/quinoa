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

#ifndef BFGS_HELPERS_H
#define BFGS_HELPERS_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief @name Functions to be used in BFGS updating.
 */
//@{

/** \brief Check that s'*y is sufficiently positive and print the result if it is not.
 *
 * @param  s           [in] DVector (size n): Secant update vector B*s=y.
 * @param  y           [in] DVector (size n): Secant update vector for B*s=y.
 * @param  sTy         [in] If sTy != NULL then *sTy must contain the value computed
 *                     from dot(s,y).  If sTy == NULL, then this value will be computed
 *                     internally.  This argument is included so that the same computation
 *                     does not have to be performed more than once by the client and
 *                     this function.
 * @param  out         [in/out] If out==NULL then no output will be printed if the
 *                     condition fails.  If out!=NULL and the function returns false
 *                     then a small amount of output will be sent to *out.
 * @param  func_name   [in] The name of the function this is being called from.
 *                     If the condition is not met then this name in included
 *                     in the printout to *out.  If func_name == NULL then this
 *                     string will obviously not be printed.
 *
 * @return If s'*y >= sqrt(mach_epsilon) * ||s||2 * ||y||2 then this function will return true.
 * Otherwise it will return false.
 */
bool BFGS_sTy_suff_p_d(
  const Vector    &s
  ,const Vector   &y
  ,const value_type     *sTy        = NULL
  ,std::ostream         *out        = NULL
  ,const char           func_name[] = NULL
  );

//@}

} // end namespace AbstractLinAlgPack

#endif // BFGS_HELPERS_H
