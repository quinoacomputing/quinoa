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

#include "MoochoPack_DecompositionSystemStateStepBuilderStd.hpp"

// NLP Stuff

#include "NLPInterfacePack_NLPSecondOrder.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "NLPInterfacePack_NLPVarReductPerm.hpp"
#include "NLPInterfacePack_NLPDirectTester.hpp"
#include "NLPInterfacePack_NLPDirectTesterSetOptions.hpp"

// Basis system and direct sparse solvers

#include "AbstractLinAlgPack_BasisSystemTester.hpp"
#include "AbstractLinAlgPack_BasisSystemTesterSetOptions.hpp"
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
#include "ConstrainedOptPack_DecompositionSystemVarReductPermStd.hpp"
#endif

// Range/null decomposition

#include "AbstractLinAlgPack_MatrixSymIdent.hpp"
#include "MoochoPack_DecompositionSystemHandlerVarReductPerm_Strategy.hpp"
#include "MoochoPack_DecompositionSystemHandlerStd_Strategy.hpp"
#include "ConstrainedOptPack_DecompositionSystemTester.hpp"
#include "ConstrainedOptPack_DecompositionSystemTesterSetOptions.hpp"
#include "ConstrainedOptPack_DecompositionSystemCoordinate.hpp"
#include "ConstrainedOptPack_DecompositionSystemOrthogonal.hpp"

// Iteration quantities

#include "ConstrainedOptPack_MatrixIdentConcatStd.hpp"               // Y, Z
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"

// Eval new point

#include "MoochoPack_EvalNewPointStd_StepSetOptions.hpp"
#include "MoochoPack_EvalNewPointTailoredApproach_StepSetOptions.hpp"
#include "MoochoPack_EvalNewPointTailoredApproachCoordinate_Step.hpp"
#include "MoochoPack_EvalNewPointTailoredApproachOrthogonal_Step.hpp"

// Other classes

#include "MoochoPack_NLPAlgoState.hpp"
#include "MoochoPack_NewDecompositionSelectionStd_Strategy.hpp"
#include "ConstrainedOptPack_VariableBoundsTesterSetOptions.hpp"
#include "NLPInterfacePack_CalcFiniteDiffProdSetOptions.hpp"
#include "NLPInterfacePack_NLPFirstDerivTester.hpp"
#include "NLPInterfacePack_NLPFirstDerivTesterSetOptions.hpp"

// Common utilities
#include "OptionsFromStreamPack_StringToIntMap.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {
  const int DEFAULT_MAX_DOF_QUASI_NEWTON_DENSE = 200;
} // end namespace

namespace MoochoPack {

//
// Here is where we define the default values for the algorithm.  These
// should agree with what are in the Moocho.opt.NLPAlgoConfigMamaJama file.
//
DecompositionSystemStateStepBuilderStd::SOptionValues::SOptionValues()
  :null_space_matrix_type_(NULL_SPACE_MATRIX_AUTO)
  ,range_space_matrix_type_(RANGE_SPACE_MATRIX_AUTO)
  ,max_dof_quasi_newton_dense_(-1)
{}

DecompositionSystemStateStepBuilderStd::DecompositionSystemStateStepBuilderStd()
{}

void DecompositionSystemStateStepBuilderStd::set_options( const options_ptr_t& options )
{
  options_ = options;
}

const DecompositionSystemStateStepBuilderStd::options_ptr_t&
DecompositionSystemStateStepBuilderStd::get_options() const
{
  return options_;
}

void DecompositionSystemStateStepBuilderStd::process_nlp_and_options(
  std::ostream      *trase_out
  ,NLP              &nlp
  ,NLPFirstOrder    **nlp_foi
  ,NLPSecondOrder   **nlp_soi
  ,NLPDirect        **nlp_fod
  ,bool             *tailored_approach
  )
{

  //
  // Probe the NLP interfaces
  //

  // Get the dimensions of the NLP
  const size_type
    n   = nlp.n(),
    m   = nlp.m(),
    r   = m, // ToDo: Compute this for real!
    dof = n - r;

  if(trase_out)
    *trase_out << "\n*** Probing the NLP object for supported interfaces ...\n";

  // Determine which NLP interface is supported
  *nlp_foi = dynamic_cast<NLPFirstOrder*>(&nlp);	
  *nlp_soi = dynamic_cast<NLPSecondOrder*>(&nlp);	
  *nlp_fod = dynamic_cast<NLPDirect*>(&nlp);
  *tailored_approach = false;
  if( *nlp_foi ) {
    if(trase_out)
      *trase_out << "\nDetected that NLP object supports the NLPFirstOrder interface!\n";
    *tailored_approach = false;
  }
  else if( *nlp_fod ) {
    if(trase_out)
      *trase_out << "\nDetected that NLP object supports the NLPDirect interface!\n";
    *tailored_approach = true;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error
      ,"NLPAlgoConfigMamaJama::config_algo_cntr(...) : Error, "
      "the NLP object of type \'" << typeName(nlp) <<
      "\' does not support the NLPFirstOrder or NLPDirect "
      "interfaces!" );
  }
  if( *nlp_soi ) {
    if(trase_out)
      *trase_out << "\nDetected that NLP object also supports the NLPSecondOrder interface!\n";
  }
  
  //
  // Process the options
  //

  // Readin the options
  if(options_.get()) {
    readin_options( *options_, &uov_, trase_out );
  }
  else {
    if(trase_out) {
      *trase_out
        << "\n*** Warning, no OptionsFromStream object was set so a default set"
          " of options will be used!\n";
    }
  }	

  // Set default
  if( uov_.max_dof_quasi_newton_dense_ < 0 )
    cov_.max_dof_quasi_newton_dense_ = DEFAULT_MAX_DOF_QUASI_NEWTON_DENSE;
  else
    cov_.max_dof_quasi_newton_dense_ = uov_.max_dof_quasi_newton_dense_;

  // Decide what type of range-space matrix to use
  if( r && uov_.range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO ) {
    const bool use_orth = (dof*dof/r) <= cov_.max_dof_quasi_newton_dense_*cov_.max_dof_quasi_newton_dense_;
    if(trase_out)
      *trase_out
        << "\nrange_space_matrix == AUTO:"
        << "\n(n-r)^2/r = (" << dof << ")^2/r = " << (dof*dof/r)
        << ( use_orth ? " <= " : " > " ) << "max_dof_quasi_newton_dense^2 = ("
        << cov_.max_dof_quasi_newton_dense_ << ")^2 = "
        << cov_.max_dof_quasi_newton_dense_*cov_.max_dof_quasi_newton_dense_
        << ( use_orth
           ? "\nsetting range_space_matrix = ORTHOGONAL\n"
           : "\nsetting range_space_matrix = COORDINATE\n" );
    cov_.range_space_matrix_type_ =
      ( use_orth
        ? RANGE_SPACE_MATRIX_ORTHOGONAL
        : RANGE_SPACE_MATRIX_COORDINATE );
  }

  // Set the default options that where not already set yet
  set_default_options(uov_,&cov_,trase_out);

}

void DecompositionSystemStateStepBuilderStd::create_decomp_sys(
  std::ostream                                       *trase_out
  ,NLP                                               &nlp
  ,NLPFirstOrder                                     *nlp_foi
  ,NLPSecondOrder                                    *nlp_soi
  ,NLPDirect                                         *nlp_fod
  ,bool                                              tailored_approach
  ,Teuchos::RCP<DecompositionSystem>         *decomp_sys
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::RCP;

  const size_type
    m  = nlp.m();

  if( m == 0 ) {
    *decomp_sys = Teuchos::null;
    return;
  }

  RCP<BasisSystem> basis_sys = Teuchos::null;
  if(!tailored_approach) {
    // Set the default basis system if one is not set
    basis_sys = nlp_foi->basis_sys();
    if( basis_sys.get() == NULL ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,"\nA basis system object was not specified by the NLPFirstOrder object of type \'"
        << typeName(nlp) << "\' and we can not build a rSQP algorithm without one!" );

    }
    // Create the testing object for the basis system and set it up.
    RCP<BasisSystemTester>
      basis_sys_tester = Teuchos::rcp(new BasisSystemTester());
    if(options_.get()) {
      BasisSystemTesterSetOptions
        opt_setter(basis_sys_tester.get());
      opt_setter.set_options(*options_);
    }
    // Determine what type of null space matrix to use
    DecompositionSystemVarReduct::EExplicitImplicit
      D_imp;
    switch( cov_.null_space_matrix_type_ ) {
      case NULL_SPACE_MATRIX_AUTO:
        D_imp = DecompositionSystemVarReduct::MAT_IMP_AUTO;
        break;
      case NULL_SPACE_MATRIX_EXPLICIT:
        D_imp = DecompositionSystemVarReduct::MAT_IMP_EXPLICIT;
        break;
      case NULL_SPACE_MATRIX_IMPLICIT:
        D_imp = DecompositionSystemVarReduct::MAT_IMP_IMPLICIT;
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
    // See if the basis system object supports basis permutations
    basis_sys_perm_ = Teuchos::rcp_dynamic_cast<BasisSystemPerm>(basis_sys);
#endif
    // Create the DecompositionSystem implementation object
    typedef RCP<DecompositionSystemVarReductImp> decomp_sys_imp_ptr_t;
    decomp_sys_imp_ptr_t decomp_sys_imp;
    switch( cov_.range_space_matrix_type_ ) {
      case RANGE_SPACE_MATRIX_COORDINATE:
        decomp_sys_imp
          = Teuchos::rcp(new DecompositionSystemCoordinate(
            nlp.space_x()
            ,nlp.space_c()
            ,basis_sys  // Will have basis_sys_->var_dep().size() == 0 if permutation
            ,basis_sys_tester
            ,D_imp
            ) );
        break;
      case RANGE_SPACE_MATRIX_ORTHOGONAL: {
        decomp_sys_imp
          = Teuchos::rcp(new DecompositionSystemOrthogonal(
            nlp.space_x()
            ,nlp.space_c()
            ,basis_sys  // Will have basis_sys_->var_dep().size() == 0 if permutation
            ,basis_sys_tester
            ) );
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);	// only a local error
    }
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
    // Create the actual DecompositionSystem object being used
    if( basis_sys_perm_.get() != NULL ) {
      if(trase_out)
        *trase_out
          << "\nThe BasisSystem object with concreate type \'" << typeName(*basis_sys)
          << "\' supports the BasisSystemPerm interface.\n"
          << "Using DecompositionSystemVarReductPermStd to support basis permutations ...\n";
      *decomp_sys = Teuchos::rcp(
        new DecompositionSystemVarReductPermStd(
          decomp_sys_imp
          ,basis_sys_perm_
          ,false // basis_selected
          ,D_imp
          ) );
    }
    else {
#endif
      if(trase_out)
        *trase_out
          << "\nThe BasisSystem object with concreate type \'" << typeName(*basis_sys)
          << "\' does not support the BasisSystemPerm interface.\n"
          << "Using " << typeName(*decomp_sys_imp) << " with a fixed basis ...\n";
      decomp_sys_imp->initialize(
        nlp.space_x()
        ,nlp.space_c()
        ,basis_sys   // Must already be ready to go with a basis selection!
        );
      *decomp_sys = decomp_sys_imp;
    }
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
  }
#endif
}

void DecompositionSystemStateStepBuilderStd::add_iter_quantities(
  std::ostream                                           *trase_out
  ,NLP                                                   &nlp
  ,NLPFirstOrder                                         *nlp_foi
  ,NLPSecondOrder                                        *nlp_soi
  ,NLPDirect                                             *nlp_fod
  ,bool                                                  tailored_approach
  ,const Teuchos::RCP<DecompositionSystem>       &decomp_sys
  ,const Teuchos::RCP<NLPAlgoState>              &state
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  
  const size_type
    n = nlp.n(),
    m = nlp.m();

  if( tailored_approach ) {
    // NLPDirect
    TEUCHOS_TEST_FOR_EXCEPT( !(  nlp_fod->con_undecomp().size() == 0  ) );
    // ToDo: Add the necessary iteration quantities when con_undecomp().size() > 0 is supported!
  }
  else {
    // NLPFirstOrder
    if(m)
      state->set_iter_quant(
        Gc_name
        ,Teuchos::rcp(
          new IterQuantityAccessContiguous<MatrixOp>(
            1
            ,Gc_name
            ,nlp_foi->factory_Gc()
            )
          )
        );
    if(nlp_soi)
      state->set_iter_quant(
        HL_name
        ,Teuchos::rcp(
          new IterQuantityAccessContiguous<MatrixSymOp>(
            1
            ,HL_name
            ,nlp_soi->factory_HL()
            )
          )
        );
  }

  //
  // Set the algorithm specific matrix objects
  //
    
  // Add range/null decomposition matrices

  if( m ) {
    if(tailored_approach) {
      // Z
      state->set_iter_quant(
        Z_name
        ,Teuchos::rcp(
          new IterQuantityAccessContiguous<MatrixOp>(
            1
            ,Z_name
            ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOp,MatrixIdentConcatStd>)
            )
          )
        );
      // Y
      state->set_iter_quant(
        Y_name
        ,Teuchos::rcp(
          new IterQuantityAccessContiguous<MatrixOp>(
            1
            ,Y_name
            ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOp,MatrixIdentConcatStd>)
            )
          )
        );
      // ToDo: Add matrix iq object for Uz
      // ToDo: Add matrix iq object for Uy
    }
    else {
      if( n > m ) {
        // Z
        state->set_iter_quant(
          Z_name
          ,Teuchos::rcp(
            new IterQuantityAccessContiguous<MatrixOp>(
              1
              ,Z_name
              ,decomp_sys->factory_Z()
              )
            )
          );
        // Y
        state->set_iter_quant(
          Y_name
          ,Teuchos::rcp(
            new IterQuantityAccessContiguous<MatrixOp>(
              1
              ,Y_name
              ,decomp_sys->factory_Y()
              )
            )
          );
        // R
        state->set_iter_quant(
          R_name
          ,Teuchos::rcp(
            new IterQuantityAccessContiguous<MatrixOpNonsing>(
              1
              ,R_name
              ,decomp_sys->factory_R()
              )
            )
          );
        // Uz
        state->set_iter_quant(
          Uz_name
          ,Teuchos::rcp(
            new IterQuantityAccessContiguous<MatrixOp>(
              1
              ,Uz_name
              ,decomp_sys->factory_Uz()
              )
            )
          );
        // Uy
        state->set_iter_quant(
          Uy_name
          ,Teuchos::rcp(
            new IterQuantityAccessContiguous<MatrixOp>(
              1
              ,Uy_name
              ,decomp_sys->factory_Uy()
              )
            )
          );
      }
      else {
        // Y
        state->set_iter_quant(
          Y_name
          ,Teuchos::rcp(
            new IterQuantityAccessContiguous<MatrixOp>(
              1
              ,Y_name
              ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSymIdent>())
              )
            )
          );
        // R
        state->set_iter_quant(
          R_name
          ,Teuchos::rcp(
            new IterQuantityAccessContiguous<MatrixOpNonsing>(
              1
              ,R_name
              ,dyn_cast<NLPFirstOrder>(nlp).basis_sys()->factory_C()
              )
            )
          );
      }
    }
  }
  else {
    // Z
    state->set_iter_quant(
      Z_name
      ,Teuchos::rcp(
        new IterQuantityAccessContiguous<MatrixOp>(
          1
          ,Z_name
          ,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSymIdent>())
          )
        )
      );
  }
  
}

void DecompositionSystemStateStepBuilderStd::create_eval_new_point(
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
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::RCP;

  const size_type
    m  = nlp.m(),
    nb = nlp.num_bounded_x();

  typedef RCP<DecompositionSystemHandler_Strategy>           decomp_sys_handler_ptr_t;
  decomp_sys_handler_ptr_t             decomp_sys_handler               = Teuchos::null;
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
  typedef RCP<DecompositionSystemHandlerSelectNew_Strategy>  decomp_sys_handler_select_new_ptr_t;
  decomp_sys_handler_select_new_ptr_t  decomp_sys_handler_select_new   = Teuchos::null;
#endif

  // Create the variable bounds testing object.
  if(nb) { // has variable bounds?
    const value_type var_bounds_warning_tol = 1e-10;
    const value_type var_bounds_error_tol   = 1e-5;
    *bounds_tester = Teuchos::rcp(
      new VariableBoundsTester(
        var_bounds_warning_tol      // default warning tolerance
        ,var_bounds_error_tol       // default error tolerance
        ) );
    if(options_.get()) {
      ConstrainedOptPack::VariableBoundsTesterSetOptions
        options_setter( bounds_tester->get() );
      options_setter.set_options(*options_);
    }
  }

  // Create the finite difference class
  *calc_fd_prod = Teuchos::rcp(new CalcFiniteDiffProd());
  if(options_.get()) {
    ConstrainedOptPack::CalcFiniteDiffProdSetOptions
      options_setter( calc_fd_prod->get() );
    options_setter.set_options(*options_);
  }

  if( m ) {

    // Decomposition system handler
    if( nlp_foi ) {
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
      if( basis_sys_perm_.get() )
        decomp_sys_handler = decomp_sys_handler_select_new
          = Teuchos::rcp( new DecompositionSystemHandlerVarReductPerm_Strategy );
      else
#endif
        decomp_sys_handler = Teuchos::rcp( new DecompositionSystemHandlerStd_Strategy );
    }
  
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
    // NewDecompositionSelectionStd_Strategy
    if( decomp_sys_handler_select_new.get() ) {
      *new_decomp_selection_strategy = Teuchos::rcp(
        new NewDecompositionSelectionStd_Strategy(decomp_sys_handler_select_new)
        );
    }
#else
    *new_decomp_selection_strategy = Teuchos::null;
#endif

  }
  else {
    decomp_sys_handler             = Teuchos::null;
    *new_decomp_selection_strategy = Teuchos::null;
  }

  //
  // EvalNewPoint_Step
  //

  // Create the step object
  if( tailored_approach ) {
    // create and setup the derivative tester
    typedef Teuchos::RCP<NLPDirectTester>   deriv_tester_ptr_t;
    deriv_tester_ptr_t
      deriv_tester = Teuchos::rcp(
        new NLPDirectTester(
          *calc_fd_prod
          ,NLPDirectTester::FD_DIRECTIONAL    // Gf testing
          ,NLPDirectTester::FD_DIRECTIONAL    // -Inv(C)*N testing
          ) );
    if(options_.get()) {
      NLPInterfacePack::NLPDirectTesterSetOptions
        options_setter(deriv_tester.get());
      options_setter.set_options(*options_);
    }
    // create the step
    typedef Teuchos::RCP<EvalNewPointTailoredApproach_Step>  _eval_new_point_step_ptr_t;
    _eval_new_point_step_ptr_t
      _eval_new_point_step = Teuchos::null;
    switch( cov_.range_space_matrix_type_ ) {
      case RANGE_SPACE_MATRIX_COORDINATE:
        _eval_new_point_step
          = Teuchos::rcp(new EvalNewPointTailoredApproachCoordinate_Step(deriv_tester,*bounds_tester));
        break;
      case RANGE_SPACE_MATRIX_ORTHOGONAL:
        _eval_new_point_step
          = Teuchos::rcp(new EvalNewPointTailoredApproachOrthogonal_Step(deriv_tester,*bounds_tester) );
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);	// only a local error
    }
    if(options_.get()) {
      EvalNewPointTailoredApproach_StepSetOptions
        options_setter(_eval_new_point_step.get());
      options_setter.set_options(*options_);
    }
    *eval_new_point_step = _eval_new_point_step;
  }
  else {
    // create and setup the derivative tester
    typedef Teuchos::RCP<NLPFirstDerivTester>   deriv_tester_ptr_t;
    deriv_tester_ptr_t
      deriv_tester = Teuchos::rcp(
        new NLPFirstDerivTester(
          *calc_fd_prod
          ,NLPFirstDerivTester::FD_DIRECTIONAL
          ) );
    if(options_.get()) {
      NLPInterfacePack::NLPFirstDerivTesterSetOptions
        options_setter(deriv_tester.get());
        options_setter.set_options(*options_);
    }
      // create and setup the decomposition system tester
    typedef Teuchos::RCP<DecompositionSystemTester>   decomp_sys_tester_ptr_t;
    decomp_sys_tester_ptr_t
      decomp_sys_tester = Teuchos::rcp( new DecompositionSystemTester() );
    if(options_.get()) {
      DecompositionSystemTesterSetOptions
        options_setter(decomp_sys_tester.get());
      options_setter.set_options(*options_);
    }
    Teuchos::RCP<EvalNewPointStd_Step>
      _eval_new_point_step = Teuchos::rcp(
        new EvalNewPointStd_Step(
          decomp_sys_handler
          ,deriv_tester
          ,*bounds_tester
          ,decomp_sys_tester
          ) );
    if(options_.get()) {
      EvalNewPointStd_StepSetOptions
        options_setter(_eval_new_point_step.get());
      options_setter.set_options(*options_);
    }
    *eval_new_point_step = _eval_new_point_step;
  }
}

// static

void DecompositionSystemStateStepBuilderStd::readin_options(
  const OptionsFromStreamPack::OptionsFromStream& options
  ,SOptionValues *ov, std::ostream* trase_out
  )
{
  namespace	ofsp = OptionsFromStreamPack;
  using		ofsp::OptionsFromStream;
  typedef		OptionsFromStream::options_group_t		options_group_t;
  using		ofsp::StringToIntMap;
  using		ofsp::StringToBool;

  TEUCHOS_TEST_FOR_EXCEPT( !( ov ) );	// only a local class error

  const std::string opt_grp_name = "DecompositionSystemStateStepBuilderStd";
  const OptionsFromStream::options_group_t optgrp = options.options_group( opt_grp_name );
  if( OptionsFromStream::options_group_exists( optgrp ) ) {

    const int num_opts = 3;
    enum EBuilder {
      NULL_SPACE_MATRIX
      ,RANGE_SPACE_MATRIX
      ,MAX_DOF_QUASI_NEWTON_DENSE
    };
    const char* SBuilder[num_opts]	= {
      "null_space_matrix"
      ,"range_space_matrix"
      ,"max_dof_quasi_newton_dense"
    };
    StringToIntMap	map( opt_grp_name, num_opts, SBuilder );

    options_group_t::const_iterator itr = optgrp.begin();
    for( ; itr != optgrp.end(); ++itr ) {
      switch( (EBuilder)map( ofsp::option_name(itr) ) ) {
        case NULL_SPACE_MATRIX:
        {
          const std::string &opt_val = ofsp::option_value(itr);
          if( opt_val == "EXPLICIT" ) {
            ov->null_space_matrix_type_ = NULL_SPACE_MATRIX_EXPLICIT;
          } else if( opt_val == "IMPLICIT" ) {
            ov->null_space_matrix_type_ = NULL_SPACE_MATRIX_IMPLICIT;
          } else if( opt_val == "AUTO" ) {
            ov->null_space_matrix_type_ = NULL_SPACE_MATRIX_AUTO;
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::invalid_argument
              ,"NLPAlgoConfigMamaJama::readin_options(...) : "
              "Error, incorrect value for \"null_space_matrix\" "
              ", Only the options for Z of EXPLICIT, IMPLICIT"
              ", and AUTO are avalible."	);
          }
          break;
        }
        case RANGE_SPACE_MATRIX:
        {
          const std::string &opt_val = ofsp::option_value(itr);
          if( opt_val == "COORDINATE" )
            ov->range_space_matrix_type_ = RANGE_SPACE_MATRIX_COORDINATE;
          else if( opt_val == "ORTHOGONAL" )
            ov->range_space_matrix_type_ = RANGE_SPACE_MATRIX_ORTHOGONAL;
          else if( opt_val == "AUTO" )
            ov->range_space_matrix_type_ = RANGE_SPACE_MATRIX_AUTO;
          else
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::invalid_argument
              ,"NLPAlgoConfigMamaJama::readin_options(...) : "
              "Error, incorrect value for \"range_space_matrix\" "
              ", Only the options for Z of COORDINATE,"
              ", ORTHOGONAL and AUTO are avalible."	);
          break;
        }
        case MAX_DOF_QUASI_NEWTON_DENSE:
          ov->max_dof_quasi_newton_dense_ = std::atoi( ofsp::option_value(itr).c_str() );
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true);	// this would be a local programming error only.
      }
    }
  }
  else {
    if(trase_out)
      *trase_out
        << "\n\n*** Warning!  The options group \"DecompositionSystemStateStepBuilderStd\" was not found.\n"
        << "Using a default set of options instead ... \n";
  }
}

void DecompositionSystemStateStepBuilderStd::set_default_options(
  const SOptionValues& uov
  ,SOptionValues *cov
  ,std::ostream* trase_out
  )
{

  TEUCHOS_TEST_FOR_EXCEPT( !( cov ) );

  if(trase_out)
    *trase_out
      << "\n*** Setting option defaults for options not set by the user or determined some other way ...\n";

  if( cov->null_space_matrix_type_ == NULL_SPACE_MATRIX_AUTO && uov.null_space_matrix_type_ == NULL_SPACE_MATRIX_AUTO ) {
    if(trase_out)
      *trase_out
        << "\nnull_space_matrix_type == AUTO: Let the algorithm deside as it goes along\n";
  }
  else if(cov->null_space_matrix_type_ == NULL_SPACE_MATRIX_AUTO) {
    cov->null_space_matrix_type_ = uov.null_space_matrix_type_;
  }
  if( cov->range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO && uov.range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO ) {
    if(trase_out)
      *trase_out
        << "\nrange_space_matrix_type == AUTO: Let the algorithm deside as it goes along\n";
  }
  else if(cov->range_space_matrix_type_ == RANGE_SPACE_MATRIX_AUTO) {
    cov->range_space_matrix_type_ = uov.range_space_matrix_type_;
  }
  if( cov->max_dof_quasi_newton_dense_ < 0 && uov.max_dof_quasi_newton_dense_ < 0 ) {
    if(trase_out)
      *trase_out
        << "\nmax_dof_quasi_newton_dense < 0 : setting max_dof_quasi_newton_dense = "<<DEFAULT_MAX_DOF_QUASI_NEWTON_DENSE<<"\n";
    cov->max_dof_quasi_newton_dense_ = DEFAULT_MAX_DOF_QUASI_NEWTON_DENSE;
  }
  else if(cov->max_dof_quasi_newton_dense_ < 0) {
    cov->max_dof_quasi_newton_dense_ = uov.max_dof_quasi_newton_dense_;
  }
  if(trase_out)
    *trase_out
      << "\n*** End setting default options\n";
}

} // end class MoochoPack
