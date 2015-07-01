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

#ifndef RSQP_TRACK_SUMMARY_STD_H
#define RSQP_TRACK_SUMMARY_STD_H

#include "MoochoPack_quasi_newton_stats.hpp"
#include "MoochoPack_qp_solver_stats.hpp"
#include "MoochoPack_act_set_stats.hpp"
#include "IterationPack_AlgorithmTracker.hpp"

namespace MoochoPack {

/** \brief This class simply outputs the convergence information
  * for each iteration.
  */
class MoochoTrackerSummaryStd
  : public IterationPack::AlgorithmTracker
{
public:

  /** \brief . */
  enum EOptError { OPT_ERROR_REDUCED_GRADIENT_LAGR, OPT_ERROR_GRADIENT_LAGR };

  /// Construct with an output stream
  MoochoTrackerSummaryStd(
    const ostream_ptr_t      &o
    ,const ostream_ptr_t     &journal_out
    ,EOptError               opt_error = OPT_ERROR_REDUCED_GRADIENT_LAGR
    );

  /// Set the output stream for summary outputting
  void set_output_stream(const ostream_ptr_t& o);

  /// Get the output stream for summary outputting.
  const ostream_ptr_t& get_output_stream() const;

  /** \brief Output the total number of qp iterations back to and
    * the k=0 iteration.
    */
  int num_total_qp_iter() const
  {	return num_total_qp_iter_;	}

  /** @name Overridden from AlgorithmTracker */
  //@{

  /** \brief . */
  void output_iteration(const Algorithm& algo) const;
  /** \brief . */
  void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;
  
  //@}

protected:

  /// Print the header to the output
  void print_header(const NLPAlgoState &s) const;

  std::ostream& o() const
  {	return *o_; }

private:
  ostream_ptr_t                  o_;
  EOptError                      opt_error_;
  mutable int                     num_total_qp_iter_;
  quasi_newton_stats_iq_member    quasi_newton_stats_;
  qp_solver_stats_iq_member       qp_solver_stats_;
  act_set_stats_iq_member		    act_set_stats_;

  // Not defined and not to be called
  MoochoTrackerSummaryStd();

};	// end class MoochoTrackerSummaryStd

}	// end namespace MoochoPack 

#endif	// RSQP_TRACK_SUMMARY_STD_H
