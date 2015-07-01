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

#ifndef QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H
#define QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace ConstrainedOptPack {

/** \brief Set options for QPSolverRelaxedQPSchur from an
  * OptionsFromStream object.
  *
  * The default options group name is QPSolverRelaxedQPSchur.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group QPSolverRelaxedQPSchur {
  *    max_qp_iter_frac  = 10.0;   *** (+dbl) max_qp_itr = max_qp_itr_frac * (# variables)
  *    max_real_runtime  = 1e+20;  *** (+dbl) maximum runtime in minutes
  *    inequality_pick_policy = ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY;
  *    inequality_pick_policy = ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY; *** not supported yet!
  *    inequality_pick_policy = ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY;
  *    bounds_tol        = 1e-10;  *** (+dbl) feasibility tolerance for bound constriants
  *    inequality_tol    = 1e-10;  *** (+dbl) feasibility tolerance for general inequality constriants
  *    equality_tol      = 1e-10;  *** (+dbl) feasibility tolerance for general equality constriants
  *    loose_feas_tol    = 1e-9;   *** (+dbl) (Expert use only)
  *    dual_infeas_tol   = 1e-12;  *** (+dbl) allowable dual infeasiblity before error
  *    huge_primal_step  = 1e+20;  *** (+dbl) value of a near infinite primal step
  *    huge_dual_step    = 1e+20;  *** (+dbl) value of a near infinite dual step
  *    bigM              = 1e+10;  *** (+dbl) value or relaxation penalty in objective
  *    warning_tol   = 1e-10;  *** Testing warning tolerance
  *    error_tol     = 1e-5;   *** Testing error tolerance
  *    iter_refine_min_iter = 1;  *** (+int) Minimum number of iterative refinement iterations
  *    iter_refine_max_iter = 3;  *** (+int) Maximum number of iterative refinement iterations
  *    iter_refine_opt_tol  = 1e-12; *** (+dbl) Optimality tolarance for iterative refinement
  *    iter_refine_feas_tol = 1e-12; *** (+dbl) Feasibility tolerance for iterative refinement
  *    iter_refine_at_solution = true; *** (+dbl) If true then iterative refinement will always be used
  *    pivot_warning_tol        = 1e-6;  *** (+dbl) Relative warning tolerance for a pivot element in the schur complement
  *    pivot_singular_tol       = 1e-8;  *** (+dbl) Relative singularity tolerance for a pivot element in the schur complement
  *    pivot_wrong_iniertia_tol = 1e-10; *** (+dbl) Relative tolerance for a pivot element in the schur complement for wrong inertia
  *    print_level = USE_INPUT_ARG;  *** Use the input argument to solve_qp(...)
  *    print_level = NO_OUTPUT;
  *    print_level = OUTPUT_BASIC_INFO;
  *    print_level = OUTPUT_ITER_SUMMARY;
  *    print_level = OUTPUT_ITER_STEPS;
  *    print_level = OUTPUT_ACT_SET;
  *    print_level = OUTPUT_ITER_QUANTITIES;
  }
  \end{verbatim}
  */
class QPSolverRelaxedQPSchurSetOptions
  : public OptionsFromStreamPack::SetOptionsFromStreamNode 
    , public OptionsFromStreamPack::SetOptionsToTargetBase<
      QPSolverRelaxedQPSchur >
{
public:

  /** \brief . */
  QPSolverRelaxedQPSchurSetOptions(
      QPSolverRelaxedQPSchur* target = 0
    , const char opt_grp_name[] = "QPSolverRelaxedQPSchur" );

protected:

  /// Overridden from SetOptionsFromStreamNode
  void setOption( int option_num, const std::string& option_value );

};	// end class QPSolverRelaxedQPSchurSetOptions

}	// end namespace ConstrainedOptPack

#endif	// QP_SOLVER_RELAXED_QP_SCHUR_SET_OPTIONS_H
