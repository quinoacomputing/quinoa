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

#ifndef MERIT_FUNC_PENALTY_PARAM_UPDATE_ADDED_STEP_H
#define MERIT_FUNC_PENALTY_PARAM_UPDATE_ADDED_STEP_H

#include "MoochoPack_Types.hpp"
#include "IterationPack_AlgorithmStep.hpp"

namespace MoochoPack {

/** \brief Base class for steps that update penalty parameters based on the
 * Lagrange multipliers lambda_k (or some approximation to them).
 *
 * This class contains methods for setting and querying values that
 * determine how the penalty parameters are updated.
 */
class MeritFunc_PenaltyParamUpdate_AddedStep
  : public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

  /** @name Pure virtual methods to be overridden by subclasses */
  //@{

  /** \brief Set the smallest value a penalty parameter is allowed to be.
    */
  virtual void small_mu( value_type small_mu ) = 0;
  /** \brief . */
  virtual value_type small_mu() const = 0;
  /** \brief Set the ratio of min(mu(i))/max(mu(i)) >= min_mu_ratio.
    *
    * If there is only one penalty parameter this is ignored.
    * The default implementation is to just return 1.
    */
  virtual void min_mu_ratio( value_type min_mu_ratio )
  {}
  /** \brief . */
  virtual value_type min_mu_ratio() const
  {	return 1.0; };
  /** \brief Set set the factor for mu = (1 + mult_factor) * abs(lambda_k(i)).
    *
    * Here it is expented that mult_factor will be very small
    * (i.e. 1e-4)
    */
  virtual void mult_factor( value_type mult_factor ) = 0;
  /** \brief . */
  virtual value_type mult_factor() const = 0;
  /** \brief Set the total KKT error ( max(||rGL||inf/max(1,||Gf||inf),||c||inf) )
    * above which the penalty parameter will be allowed to be reduced.
    *
    * If the KKT error is less than near_solution the penalty parameter will
    * only be increased.  This is usually needed to prove convergence.
    */
  virtual void kkt_near_sol( value_type kkt_near_sol ) = 0;
  /** \brief . */
  virtual value_type kkt_near_sol() const = 0;

  //@}

};	// end class MeritFunc_PenaltyParamUpdate_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_PENALTY_PARAM_UPDATE_ADDED_STEP_H
