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

#ifndef ALGORITHM_STEP_H
#define ALGORITHM_STEP_H

#include <iosfwd>

#include "IterationPack_Types.hpp"

namespace IterationPack {

/** \brief Base type for all objects that perform steps in an <tt>Algorithm</tt>.
  */
class AlgorithmStep {
public:

  /** \brief . */
  typedef size_t poss_type;

  /** @name Pure virtual functions that must be overridden */
  //@{

  /** \brief Called by <tt>Algorithm</tt> to perform a main, pre or post step at step_poss and assoc_step_poss.
    *
    * @return Should return false if this step object has terminated the algorithm or
    * redirected control to another step.  In this case it is assumed that <tt>this</tt> called
    * <tt>algo\ref Algorithm::terminate ".terminate(...)"</tt>
    * or <tt>algo\ref Algorithm::do_step_next ".do_step_next(...)"</tt>.
    */
  virtual bool do_step(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    ) = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief . */
  virtual ~AlgorithmStep() {}

  /** \brief Called by <tt>Algorithm</tt> just before the algorithm is run.
   *
   * This allows step objects to reinitialize themselves just before
   * an algorithm is run.
   *
   * The default implementation does nothing.
   */
  virtual void initialize_step(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    )
    {}

  /** \brief Called by <tt>Algorithm</tt> to inform when a runtime configuration change
   * is finihed.
   *
   * This function is only called when the algorithm is already running
   * but the configuration has changed.
   *
   * The default implementation does nothing.
   */
  virtual void inform_updated(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    )
    {}

  /** \brief Called by <tt>Algorithm</tt> just after an algorithm is terminiated.
   *
   * This allows step objects to perform any final processing or cleanup
   * just after an algorithm is finished.
   *
   * The default implementation does nothing.
   */
  virtual void finalize_step(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    )
    {}

  /** \brief Called by <tt>Algorithm::print_algorithm()</tt> to print out what this step does in Matlab like format.
    *
    * The default does nothing.
    */
  virtual void print_step(
    const Algorithm&         algo
    ,poss_type               step_poss
    ,EDoStepType             type
    ,poss_type               assoc_step_poss
    , std::ostream&          out
    ,const std::string&      leading_str
    ) const
  {}

  //@}

};	// end class AlgorithmStep

}	// end namespace IterationPack 

#endif // ALGORITHM_STEP_H
