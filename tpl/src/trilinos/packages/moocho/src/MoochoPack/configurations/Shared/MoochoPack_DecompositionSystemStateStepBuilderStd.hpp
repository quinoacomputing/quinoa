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

#ifndef DECOMPOSITION_SYSTEM_STATE_STEP_BUILDER_STD_H
#define DECOMPOSITION_SYSTEM_STATE_STEP_BUILDER_STD_H

#include "MoochoPack_Types.hpp"
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
#include "AbstractLinAlgPack_BasisSystemPerm.hpp"
#endif
#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"

namespace OptionsFromStreamPack {
  class OptionsFromStream;
}

namespace MoochoPack {

/** \brief Standard builder object for creating DecompositionSystem, EvalNewPoint Step and other objects
 * and setting up some of the state object.
 *
 * This class is designed to be used by <tt>NLPAlgoConfig</tt> subclasses based on SQP
 * and performs many different tasks that are common to all of these algorithms.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemStateStepBuilderStd
{
public:
  
  /** @name Enums for variaous options categories */
  //@{

  /** \brief . */
  enum ENullSpaceMatrixType {
    NULL_SPACE_MATRIX_AUTO, NULL_SPACE_MATRIX_EXPLICIT
    , NULL_SPACE_MATRIX_IMPLICIT };
  /** \brief . */
  enum ERangeSpaceMatrixType {
    RANGE_SPACE_MATRIX_AUTO, RANGE_SPACE_MATRIX_COORDINATE
    , RANGE_SPACE_MATRIX_ORTHOGONAL };

  //@}

  /** @name Struct for options values */
  //@{

  /** \brief . */
  struct SOptionValues {
    // Constructor (sets default values)
    SOptionValues();
    ENullSpaceMatrixType	null_space_matrix_type_;
    ERangeSpaceMatrixType	range_space_matrix_type_;
    int						max_dof_quasi_newton_dense_;    // If < 0, don't change default
  };

  //@}
  
  /** \brief . */
  typedef Teuchos::RCP<
    const OptionsFromStreamPack::OptionsFromStream>             options_ptr_t;

  /** \brief . */
  DecompositionSystemStateStepBuilderStd();

  /** \brief Set the options that will be used to configure the algorithmic objects.
   *
   *  @param  options
   *               [in] If \c NULL then no options will be set.  If <tt>!=NULL</tt> then
   *               this is the \c OptionsFromStream object that will be used to extract the
   *               options to use for the algorithm.  The state of this object must
   *               be maintained by the client until this object is no longer needed.
   */
  void set_options( const options_ptr_t& options );
  /** \brief . */
  const options_ptr_t& get_options() const;
  /** \brief Process the %NLP and process the options passed in from <tt>set_options()</tt>.
   * Postconditions:<ul>
   * <li> <tt>this->current_option_values()</tt> returns the options that will be
   *      used in all of the following method.
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  void process_nlp_and_options(
    std::ostream          *trase_out
    ,NLP                  &nlp
    ,NLPFirstOrder        **nlp_foi
    ,NLPSecondOrder       **nlp_soi
    ,NLPDirect            **nlp_fod
    ,bool                 *tailored_approach
    );
  /** \brief Create the decomposition system object.
   *
   * ToDo: Finish documentation!
   */
  void create_decomp_sys(
    std::ostream                                       *trase_out
    ,NLP                                               &nlp
    ,NLPFirstOrder                                     *nlp_foi
    ,NLPSecondOrder                                    *nlp_soi
    ,NLPDirect                                         *nlp_fod
    ,bool                                              tailored_approach
    ,Teuchos::RCP<DecompositionSystem>         *decomp_sys
    );
  /** \brief Add the common iteration quantities to the state object.
   * 
   * ToDo: Finish documentation!
   */
  void add_iter_quantities(
    std::ostream                                           *trase_out
    ,NLP                                                   &nlp
    ,NLPFirstOrder                                         *nlp_foi
    ,NLPSecondOrder                                        *nlp_soi
    ,NLPDirect                                             *nlp_fod
    ,bool                                                  tailored_approach
    ,const Teuchos::RCP<DecompositionSystem>       &decomp_sys
    ,const Teuchos::RCP<NLPAlgoState>              &state
    );

  /** \brief Create the EvalNewPoint step object and allocated objects.
   *
   * ToDo: Finish documentation!
   */
  void create_eval_new_point(
    std::ostream                                                 *trase_out
    ,NLP                                                         &nlp
    ,NLPFirstOrder                                               *nlp_foi
    ,NLPSecondOrder                                              *nlp_soi
    ,NLPDirect                                                   *nlp_fod
    ,bool                                                        tailored_approach
    ,const Teuchos::RCP<DecompositionSystem>             &decomp_sys
    ,Teuchos::RCP<IterationPack::AlgorithmStep>          *eval_new_point_step
    ,Teuchos::RCP<CalcFiniteDiffProd>                    *calc_fd_prod
    ,Teuchos::RCP<VariableBoundsTester>                  *bounds_tester
    ,Teuchos::RCP<NewDecompositionSelection_Strategy>    *new_decomp_selection_strategy
    );

  /** \brief Return the current option values being used.
   */
  SOptionValues& current_option_values();

private:
  
  // ///////////////////////////
  // Private data members

  /// Smart pointer to options
  options_ptr_t   options_;

  /// Options structs
  SOptionValues       uov_; // options set by user
  SOptionValues       cov_; // current option values actually used

#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
  Teuchos::RCP<BasisSystemPerm> basis_sys_perm_;
#endif

  // /////////////////////////
  // Private member functions

  /// Read in the options from a stream
  static void readin_options(
    const OptionsFromStreamPack::OptionsFromStream& options
    ,SOptionValues *option_values, std::ostream* trase_out
    );

  /// Set the defaults for options not set by the user
  static void set_default_options(
    const SOptionValues& user_option_values
    ,SOptionValues *current_option_values
    ,std::ostream* trase_out
    );

}; // end class DecompositionSystemStateStepBuilderStd

// ///////////////////////////////
// Inline members

inline
DecompositionSystemStateStepBuilderStd::SOptionValues&
DecompositionSystemStateStepBuilderStd::current_option_values()
{
  return cov_;
}

}  // end namespace MoochoPack

#endif // DECOMPOSITION_SYSTEM_STATE_STEP_BUILDER_STD_H
