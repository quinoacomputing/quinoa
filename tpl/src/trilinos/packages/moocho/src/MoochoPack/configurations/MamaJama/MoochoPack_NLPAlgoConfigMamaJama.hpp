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

#ifndef RSQP_ALGO_CONFIG_MAMA_JAMA_H
#define RSQP_ALGO_CONFIG_MAMA_JAMA_H

#include "MoochoPack_NLPAlgoConfig.hpp"
#include "MoochoPack_NLPAlgo.hpp"
#include "MoochoPack_DecompositionSystemStateStepBuilderStd.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"

namespace MoochoPack {

/** \brief This is a do all configuration class for <tt>NLPAlgo</tt>.
 *
 * This class relies on the builder class <tt>DecompositionSystemStateStepBuilderStd</tt>
 * to perform many different tasks.
 *
 * Options specific for to this configuration class and the classes that
 * it works with that can be set through <tt>this->set_options()</tt>, see the files
 * <tt>\ref DecompositionSystemStateStepBuilderStd_opts "Moocho.opt.DecompositionSystemStateStepBuilderStd"</tt>.
 * and <tt>\ref NLPAlgoConfigMamaJama_opts "Moocho.opt.NLPAlgoConfigMamaJama"</tt>.
 *
 * Note that all built-in support for basis permutations and direct sparse solvers
 * can be left out if the macro MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS is defined
 * This will result in smaller executables for programs that don't need this
 * extra functionality..
 * 
 * ToDo: Finish documentation!
 */
class NLPAlgoConfigMamaJama : public NLPAlgoConfig {
public:

  /** \brief . */
  NLPAlgoConfigMamaJama();

  /** \brief . */
  ~NLPAlgoConfigMamaJama();

  /** Overridden from NLPAlgoConfig */
  //@{

  /** \brief Set the <tt>OptionsFromStream</tt> object that will be used for specifying the options.
   *
   *  @param  options
   *               [in] If \c NULL then no options will be set.  If <tt>!=NULL</tt> then
   *               this is the \c OptionsFromStream object that will be used to extract the
   *               options to use for the algorithm.  The state of this object must
   *               be maintained by the client until \c config_algo_cntr() is called
   *               and it is at this point that the options are read.
   *
   */
  void set_options( const options_ptr_t& options );
  /** \brief . */
  const options_ptr_t& get_options() const;
  /** \brief . */
  void config_algo_cntr(NLPAlgoContainer* algo_cntr, std::ostream* trase_out);
  /** \brief . */
  void init_algo(NLPAlgoInterface* algo);

  //@}

public:

  /** @name Enums for variaous options categories */
  //@{

  /** \brief . */
  enum EQuasiNewton {
    QN_AUTO, QN_BFGS, QN_PBFGS, QN_LBFGS, QN_LPBFGS };
  /** \brief . */
  enum EHessianInitialization {
    INIT_HESS_AUTO, INIT_HESS_SERIALIZE, INIT_HESS_IDENTITY, INIT_HESS_FIN_DIFF_SCALE_IDENTITY
    , INIT_HESS_FIN_DIFF_SCALE_DIAGONAL, INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS };
  /** \brief . */
  enum EQPSolverType {
    QP_AUTO, QP_QPSOL, QP_QPOPT, QP_QPKWIK, QP_QPSCHUR };
  /** \brief . */
  enum ELineSearchMethod {
    LINE_SEARCH_AUTO, LINE_SEARCH_NONE, LINE_SEARCH_DIRECT
    , LINE_SEARCH_2ND_ORDER_CORRECT, LINE_SEARCH_WATCHDOG
    , LINE_SEARCH_FILTER };
  /** \brief . */
  enum EMeritFunctionType {
    MERIT_FUNC_AUTO, MERIT_FUNC_L1, MERIT_FUNC_MOD_L1
    , MERIT_FUNC_MOD_L1_INCR };
  /** \brief . */
  enum EL1PenaltyParamUpdate {
    L1_PENALTY_PARAM_AUTO, L1_PENALTY_PARAM_WITH_MULT
    , L1_PENALTY_PARAM_MULT_FREE };

  //@}

  /** @name Struct for options values */
  //@{

  /** \brief . */
  struct SOptionValues {
    // Constructor (sets default values)
    SOptionValues();
    // Variable Reduction, Range/Null space decompositions
    value_type max_basis_cond_change_frac_; // If < , don't change default
    // Reduced Hessian Approximations
    bool exact_reduced_hessian_;
    EQuasiNewton quasi_newton_;
    int num_lbfgs_updates_stored_; // If < 0, don't change default
    bool lbfgs_auto_scaling_;
    EHessianInitialization hessian_initialization_;
    // QP subproblem solvers
    EQPSolverType qp_solver_type_;
    bool reinit_hessian_on_qp_fail_;
    // Line search methods
    ELineSearchMethod line_search_method_;
    EMeritFunctionType merit_function_type_;
    EL1PenaltyParamUpdate l1_penalty_param_update_;
    int full_steps_after_k_; // If < 0, do not use this option at all.
    value_type max_pz_norm_;
    int num_pz_damp_iters_;
  };
  
  //@}

private:

  /// Builder class for some common code
  DecompositionSystemStateStepBuilderStd decomp_sys_step_builder_;
  
  /// Smart pointer to options
  options_ptr_t options_;
  
  /// Options structs
  SOptionValues uov_; // options set by user
  SOptionValues cov_; // current option values actually used

  // ///////////////////////////////////////////////////////
  // Private member functions

  /// Read in the options from a stream
  static void readin_options(
    const OptionsFromStreamPack::OptionsFromStream& options
    , SOptionValues *option_values, std::ostream* trase_out );

  /// Set the defaults for options not set by the user
  static void set_default_options(
    const SOptionValues& user_option_values
    , SOptionValues *current_option_values
    , std::ostream* trase_out );

};	// end class NLPAlgoConfigMamaJama

/** \defgroup NLPAlgoConfigMamaJama_opts Options for NLPAlgoConfigMamaJama.
 *
 * The following is the contents of the file <tt>Moocho.opt.NLPAlgoConfigMamaJama</tt>
 * which are options specific to the class <tt>MoochoPack::NLPAlgoConfigMamaJama</tt>
 * and the class objects that it configures.
 *
 * \verbinclude configurations/MamaJama/Moocho.opt.NLPAlgoConfigMamaJama
 */

}	// end namespace MoochoPack 

#endif	// RSQP_ALGO_CONFIG_MAMA_JAMA_H
