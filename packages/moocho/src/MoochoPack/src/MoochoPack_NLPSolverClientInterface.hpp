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

#ifndef RSQP_SOLVER_CLIENT_INTERFACE_H
#define RSQP_SOLVER_CLIENT_INTERFACE_H

#include <stdexcept>

#include "MoochoPack_Types.hpp"
#include "NLPInterfacePack_NLP.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \brief This is the most basic interface that clients use to solve an NLP.
 *
 * ToDo: Finish documentaiton.
 */
class NLPSolverClientInterface {
public:

  /** @name Public Types */
  //@{

  /** \brief . */
  enum EFindMinReturn {
    SOLUTION_FOUND
    ,MAX_ITER_EXCEEDED
    ,MAX_RUN_TIME_EXCEEDED
    ,ALGORITHMIC_ERROR
  };

  /// Thrown if the setup is not valid
  class InvalidSetup : public std::logic_error
  {public: InvalidSetup(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}
  
  /** @name Solver Parameters */
  //@{

  /// Set the maximum number of iterations the rSQP algorithm can perform
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, max_iter );

  /** \brief Set the maximum run_time
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( double, max_run_time );

  /** \brief Set the termination tolerance for the relative (scaled) linear dependence of the
   * gradients part of the first order necessary optimality conditions.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_tol );

  /** \brief Set the termination tolerance for the (scaled) equality constraints ||c(x*)||inf
   * which is part of the first order necessary optimality conditions.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_tol );

  /** \brief Set the termination tolerance for the complementarity condition 
   *  for the (scaled) bound constraints
   *  which is part of the first order necessary optimality conditions.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, comp_tol );

  /** \brief Set the termination tolerance for the change in the estimate of the solution.
   *
   * The test is: <tt>|d(i)|/(1+|x(i)|) < step_tol</tt>. 
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, step_tol );

  /** \brief Determine the amount of output to a journal file.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EJournalOutputLevel, journal_output_level );

  /** \brief Determine the amount of output of the null space to a journal file.
   *
   * This option allows the user to perform a higher level of output
   * for quantities in the null space.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EJournalOutputLevel, null_space_journal_output_level );

  /** \brief Set the precesion of the journal output.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, journal_print_digits );

  /** \brief Set whether computations will be double checked or not.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_results );

  /** \brief Set whether the condition numbers of important matrics is
   * computed and printed or not.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, calc_conditioning );

  /** \brief Set whether or not matrix norms are computed and printed.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, calc_matrix_norms );

  /** \brief Set whether calc_conditioning and calc_matrix_norms apply to only
   * null space matrices.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, calc_matrix_info_null_space_only );

  //@}

  /** @name Constructors/initalizers */
  //@{

  /// <<std comp>> members for the nlp
  STANDARD_COMPOSITION_MEMBERS( NLP, nlp );

  /// <<std comp>> members for the track
  STANDARD_COMPOSITION_MEMBERS( AlgorithmTracker, track );

  /** \brief Construct with no references set to nlp or track objects.
   */
  NLPSolverClientInterface(
    int                    max_iter             = 10000
    ,double                max_run_time         = 1e+10 // run forever
    ,value_type            opt_tol              = 1e-6
    ,value_type            feas_tol             = 1e-6
    ,value_type            comp_tol             = 1e-6
    ,value_type            step_tol             = 1e-2
    ,EJournalOutputLevel   journal_output_level = PRINT_ALGORITHM_STEPS
    ,EJournalOutputLevel   null_space_journal_output_level = PRINT_ALGORITHM_STEPS
    ,int                   journal_print_digits = 6
    ,bool                  check_results        = false
    ,bool                  calc_conditioning    = false
    ,bool                  calc_matrix_norms    = false
    ,bool                  calc_matrix_info_null_space_only = false
    );

  /** \brief . */
  virtual ~NLPSolverClientInterface() {}

  //@}

  /** @name Solve the NLP*/
  //@{

  /** \brief Find the minimun of the set NLP.
   *
   * This function returns <tt>SOLUTION_FOUND</tt> if the NLP has been solved to the desired
   * tolerances.  In this case <tt>this->track().output_final(...,TERMINATE_TRUE)</tt>
   * and <tt>this->nlp().report_final_solution(...,true)</tt> are called before this function
   * returns..
   *
   * If the solution is not found, then <tt>this->nlp().report_final_solution(...,false)</tt> is
   * called and one of the following occurs:
   * <ul>
   * <li>  If the maximum number of iterations has been exceeded then <tt>MAX_ITER_EXCEEDED</tt> will
   * be returned.  In this case <tt>this->track().output_final(...,MAX_ITER_EXCEEDED)</tt> is called.
   * <li> If the maximum runtime has been exceeded then <tt>MAX_RUN_TIME_EXCEEDED</tt> will be
   * returned.  In this case <tt>this->track().output_final(...,MAX_RUN_TIME_EXCEEDED)</tt> is called.
   * <li> An exception is thrown.  The client should be prepaired to catch any exceptions thrown
   * from this function.
   * All of the purposefully thrown exceptions are derived from std::exception so the
   * client can check the what() function to see and description of the error.  If an
   * exception is thrown then <tt>this->track().output_final(...,TERMINATE_FALSE)</tt> will
   * be called before this exception is rethrown out of this functiion.  If the constraints
   * are found to be infeasible, then the exception <tt>InfeasibleConstraints</tt> will be thrown.
   * If a line search failure occurs then the exception <tt>LineSearchFailure</tt> will be thrown.
   * If some test failed then the exception <tt>TestFailed</tt> will be thrown.  Many other exceptions
   * may be thrown but these are the main ones that the SQP algorithms known about and will
   * purposefully generate.
   * </ul>
   *
   * Preconditions:<ul>
   * <li> <tt>this->nlp() != 0</tt> (throw InvalidSetup)
   * </ul>
   *
   * Postcondtions:<ul>
   * <li> Minimum of %NLP is found to opt_tol, max_iter was reached
   * or max_run_time reached (throw std::exection)
   * </ul>
   */
  virtual EFindMinReturn find_min() = 0;

  //@}

  /** @name Algorithm description */
  //@{

  /** \brief Prints a description of the algorithm.
    */
  virtual void print_algorithm(std::ostream& out) const = 0;

  //@}

  /** @name Algorithm timing */
  //@{

  /** \brief Causes algorithm to be timed.
   *
   * Call with <tt>algo_timing == true</tt> before calling <tt>find_min()</tt>
   * to have the algorithm timed.
   */
  virtual void set_algo_timing( bool algo_timing ) = 0;

  /** \brief . */
  virtual bool algo_timing() const = 0;

  /** \brief Outputs table of times for each step and the cummulative times.
   *
   * Call after <tt>find_min()</tt> has executed to get a table
   * of times.
   */
  virtual void print_algorithm_times( std::ostream& out ) const = 0;

  //@}

private:

#ifdef DOXYGEN_COMPILE // Strictly for doxygen diagrams
  /** \brief . */
  NLPInterfacePack::NLP                  *nlp;
  /** \brief . */
  IterationPack::AlgorithmTracker   *track;
#endif

}; // end class NLPSolverClientInterface

}  // end namespace MoochoPack

#endif	// RSQP_SOLVER_CLIENT_INTERFACE_H
