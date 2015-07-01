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
//

#ifndef EXAMPLE_NLP_FIRST_ORDER_INFO_RUN_H
#define EXAMPLE_NLP_FIRST_ORDER_INFO_RUN_H

#include <iosfwd>

#include "NLPInterfacePack_Types.hpp"
#include "MoochoPack_MoochoSolver.hpp"

namespace NLPInterfacePack {

/** \defgroup ExampleNLPFirstOrderRun_grp Helper function for ExampleNLPFirstOrder */
//@{

/** \brief Function accepts a VectorSpace object and then uses it to define
 * an example NLP and run <tt>MoochoPack::MoochoSolver</tt> on it.
 *
 * @param  vec_space   [in] The vector space object used to create all of the
 *                     needed vector spaces and vectors.  This vector space and
 *                     the vectors it creates will get a though testing.
 * @param  xo          [in] The initial starting point for unknown variables (before
 *                     they are forced in bounds).
 * @param  has_bounds  [in] If true, then the NLP will have bounds on the variables.
 * @param  dep_bounded [in] (valid only if has_bounds == true) If true, then
 *                     the dependent variables will be bounded, if false the
 *                     independent variables will be bounded.
 * @param  console_out [in/out] If != NULL then *console_out gets the output.
 * @param  error_out   [in/out] If != NULL then *eout gets minimal summary output.
 * @param  throw_solve_exception
 *                     [in] If true then solver will not throw exception (but other code may).
 * @param  algo_out    [in/out] If != NULL then it gets algo outptut, otherwise goes to 'MoochoAlgo.out'
 * @param  summary_out [in/out] If != NULL then it gets summary outptut, otherwise goes to 'MoochoSummary.out'
 * @param  journal_out [in/out] If != NULL then it gets journal outptut, otherwise goes to 'MoochoJournal.out'
 *
 * @returns Returns the return value from <tt>MoochoPack::rsqp_mama_jama_solve()</tt>
 * (see this function for most of the documentation).
 */
MoochoPack::MoochoSolver::ESolutionStatus
ExampleNLPFirstOrderRun(
  const VectorSpace&   vec_space
  ,value_type          xo
  ,bool                has_bounds
  ,bool                dep_bounded
  ,std::ostream*       console_out
  ,std::ostream*       error_out
  ,bool                throw_solve_exception = false
  ,std::ostream*       algo_out              = NULL
  ,std::ostream*       summary_out           = NULL
  ,std::ostream*       journal_out           = NULL
  );

//@}

} // end namespace NLPInterfacePack

#endif // EXAMPLE_NLP_FIRST_ORDER_INFO_RUN_H


