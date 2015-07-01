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

#ifndef RSQP_ALGO_INTERFACE_H
#define RSQP_ALGO_INTERFACE_H

#include "MoochoPack_Types.hpp"
#include "MoochoPack_NLPSolverClientInterface.hpp"

namespace MoochoPack {

/** \brief Interface \c NLPAlgoContainer uses to access \c NLPAlgo.
 *
 * This interface helps avoid dangerous usage stategies for an \c NLPAlgo object.
 */
class NLPAlgoInterface {
public:

  /** \brief . */
  virtual ~NLPAlgoInterface() {}

  /** \brief Print the algorithm description.
   */
  virtual void interface_print_algorithm(std::ostream& out) const = 0;

  /** \brief Start the iterations.
    *
    * This function returns true if the solution was found and false
    * if the maximum number of iterations was reached before the
    * solution was found.
    */
  virtual NLPSolverClientInterface::EFindMinReturn dispatch() = 0;

  /** \brief Return the state object.
   */
  virtual const NLPAlgoState& retrieve_state() const = 0;

  /** @name Algorithm timing */
  //@{

  /** \brief . */
  virtual void interface_set_algo_timing( bool algo_timing ) = 0;
  /** \brief . */
  virtual bool interface_algo_timing() const = 0;
  /** \brief . */
  virtual void interface_print_algorithm_times( std::ostream& out ) const = 0;

  //@}

};	// end class NLPAlgoInterface

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_INTERFACE_H
