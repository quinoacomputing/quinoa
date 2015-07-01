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

#include "MoochoPack_NLPSolverClientInterface.hpp"
#include "IterationPack_AlgorithmTracker.hpp"

MoochoPack::NLPSolverClientInterface::NLPSolverClientInterface(
  int                      max_iter
  ,double                  max_run_time
  ,value_type              opt_tol
  ,value_type              feas_tol
  ,value_type              comp_tol
  ,value_type              step_tol
  ,EJournalOutputLevel     journal_output_level
  ,EJournalOutputLevel     null_space_journal_output_level
  ,int                     journal_print_digits
  ,bool                    check_results
  ,bool                    calc_conditioning
  ,bool                    calc_matrix_norms
  ,bool                    calc_matrix_info_null_space_only
  )
  :max_iter_(max_iter)
  ,max_run_time_(max_run_time)
  ,opt_tol_(opt_tol)
  ,feas_tol_(feas_tol)
  ,comp_tol_(comp_tol)
  ,step_tol_(step_tol)
  ,journal_output_level_(journal_output_level)
  ,null_space_journal_output_level_(null_space_journal_output_level)
  ,journal_print_digits_(journal_print_digits)
  ,check_results_(check_results)
  ,calc_conditioning_(calc_conditioning)
  ,calc_matrix_norms_(calc_matrix_norms)
  ,calc_matrix_info_null_space_only_(calc_matrix_info_null_space_only)
{}
