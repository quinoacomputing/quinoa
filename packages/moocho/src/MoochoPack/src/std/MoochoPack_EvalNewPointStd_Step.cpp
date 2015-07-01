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

#include <ostream>
#include <typeinfo>

#include "MoochoPack_EvalNewPointStd_Step.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "ConstrainedOptPack_DecompositionSystemVarReduct.hpp"
#include "AbstractLinAlgPack_MatrixSymIdent.hpp"
#include "AbstractLinAlgPack_PermutationOut.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

#ifdef TEUCHOS_DEBUG
#include "DenseLinAlgPack_PermVecMat.hpp"
#endif

namespace MoochoPack {

EvalNewPointStd_Step::EvalNewPointStd_Step(
  const decomp_sys_handler_ptr_t                              &decomp_sys_handler
  ,const deriv_tester_ptr_t                                   &deriv_tester
  ,const bounds_tester_ptr_t                                  &bounds_tester
  ,const decomp_sys_tester_ptr_t                              &decomp_sys_tester
  ,EFDDerivTesting                                            fd_deriv_testing
  ,DecompositionSystemHandler_Strategy::EDecompSysTesting     decomp_sys_testing
  ,DecompositionSystemHandler_Strategy::EDecompSysPrintLevel  decomp_sys_testing_print_level
  )
  :decomp_sys_handler_(decomp_sys_handler)
  ,deriv_tester_(deriv_tester)
  ,bounds_tester_(bounds_tester)
  ,decomp_sys_tester_(decomp_sys_tester)
  ,fd_deriv_testing_(fd_deriv_testing)
  ,decomp_sys_testing_(decomp_sys_testing)
  ,decomp_sys_testing_print_level_(decomp_sys_testing_print_level)
{}

bool EvalNewPointStd_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  using Teuchos::rcp;
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using IterationPack::print_algorithm_step;
  using NLPInterfacePack::NLPFirstOrder;

  NLPAlgo         &algo   = rsqp_algo(_algo);
  NLPAlgoState    &s      = algo.rsqp_state();
  NLPFirstOrder   &nlp    = dyn_cast<NLPFirstOrder>(algo.nlp());

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  EJournalOutputLevel ns_olevel = algo.algo_cntr().null_space_journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  if(!nlp.is_initialized())
    nlp.initialize(algo.algo_cntr().check_results());

  Teuchos::VerboseObjectTempState<NLP>
    nlpOutputTempState(rcp(&nlp,false),Teuchos::getFancyOStream(rcp(&out,false)),convertToVerbLevel(olevel));

  const size_type
    n  = nlp.n(),
    nb = nlp.num_bounded_x(),
    m  = nlp.m();
  size_type
    r  = 0;

  // Get the iteration quantity container objects
  IterQuantityAccess<value_type>
    &f_iq   = s.f();
  IterQuantityAccess<VectorMutable>
    &x_iq   = s.x(),
    *c_iq   = m > 0  ? &s.c() : NULL,
    &Gf_iq  = s.Gf();
  IterQuantityAccess<MatrixOp>
    *Gc_iq  = m  > 0 ? &s.Gc() : NULL,
    *Z_iq   = NULL,
    *Y_iq   = NULL,
    *Uz_iq  = NULL,
    *Uy_iq  = NULL;
  IterQuantityAccess<MatrixOpNonsing>
    *R_iq   = NULL;

  MatrixOp::EMatNormType mat_nrm_inf = MatrixOp::MAT_NORM_INF;
  const bool calc_matrix_norms = algo.algo_cntr().calc_matrix_norms();
  const bool calc_matrix_info_null_space_only = algo.algo_cntr().calc_matrix_info_null_space_only();
  
  if( x_iq.last_updated() == IterQuantity::NONE_UPDATED ) {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "\nx is not updated for any k so set x_k = nlp.xinit() ...\n";
    }
    x_iq.set_k(0) = nlp.xinit();
  }
  
  // Validate x
  if( nb && algo.algo_cntr().check_results() ) {
    assert_print_nan_inf(
      x_iq.get_k(0), "x_k", true
      , int(olevel) >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
      );
    if( nlp.num_bounded_x() > 0 ) {
      if(!bounds_tester().check_in_bounds(
           int(olevel)  >= int(PRINT_ALGORITHM_STEPS) ? &out : NULL
           ,int(olevel) >= int(PRINT_VECTORS)                // print_all_warnings
           ,int(olevel) >= int(PRINT_ITERATION_QUANTITIES)  // print_vectors
           ,nlp.xl(),        "xl"
           ,nlp.xu(),        "xu"
           ,x_iq.get_k(0),   "x_k"
           ))
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, TestFailed
          ,"EvalNewPointStd_Step::do_step(...) : Error, "
          "the variables bounds xl <= x_k <= xu where violated!" );
      }
    }
  }

  Vector &x = x_iq.get_k(0);

  Range1D  var_dep(Range1D::INVALID), var_indep(Range1D::INVALID);
  if( s.get_decomp_sys().get() ) {
    const ConstrainedOptPack::DecompositionSystemVarReduct
      *decomp_sys_vr = dynamic_cast<ConstrainedOptPack::DecompositionSystemVarReduct*>(&s.decomp_sys());
    if(decomp_sys_vr) {
      var_dep   = decomp_sys_vr->var_dep();
      var_indep = decomp_sys_vr->var_indep();
    }
    s.var_dep(var_dep);
    s.var_indep(var_indep);
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\n||x_k||inf            = " << x.norm_inf();
    if( var_dep.size() )
      out << "\n||x(var_dep)_k||inf   = " << x.sub_view(var_dep)->norm_inf();
    if( var_indep.size() )
      out << "\n||x(var_indep)_k||inf = " << x.sub_view(var_indep)->norm_inf();
    out << std::endl;
  }
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out << "\nx_k = \n" << x;
    if( var_dep.size() )
      out << "\nx(var_dep)_k = \n" << *x.sub_view(var_dep);
  }
  if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    if( var_indep.size() )
      out << "\nx(var_indep)_k = \n" << *x.sub_view(var_indep);
  }

  // Set the references to the current point's quantities to be updated
  const bool f_k_updated  = f_iq.updated_k(0);
  const bool Gf_k_updated = Gf_iq.updated_k(0);
  const bool c_k_updated  = m  > 0 ? c_iq->updated_k(0)  : false;
  const bool Gc_k_updated = m  > 0 ? Gc_iq->updated_k(0) : false;
  nlp.unset_quantities();
  if(!f_k_updated) nlp.set_f( &f_iq.set_k(0) );
  if(!Gf_k_updated) nlp.set_Gf( &Gf_iq.set_k(0) );
  if( m > 0 ) {
    if(!c_k_updated) nlp.set_c( &c_iq->set_k(0) );
    if(!Gc_k_updated) nlp.set_Gc( &Gc_iq->set_k(0) );
  }

  // Calculate Gc at x_k
  bool new_point = true;
  if(m > 0) {
    if(!Gc_k_updated) nlp.calc_Gc( x, new_point );
    new_point = false;
  }

  //
  // Update (or select a new) range/null decomposition
  //
  bool new_decomp_selected = false;
  if( m > 0 ) {
    
    // Update the range/null decomposition
    decomp_sys_handler().update_decomposition(
      algo, s, nlp, decomp_sys_testing(), decomp_sys_testing_print_level()
      ,&new_decomp_selected
      );

    r  = s.equ_decomp().size();

    Z_iq   = ( n > m && r > 0 )      ? &s.Z()  : NULL;
    Y_iq   = ( r > 0 )               ? &s.Y()  : NULL;
    Uz_iq  = ( m  > 0 && m  > r )    ? &s.Uz() : NULL;
    Uy_iq  = ( m  > 0 && m  > r )    ? &s.Uy() : NULL;
    R_iq   = ( m > 0 )               ? &s.R()  : NULL;

    // Determine if we will test the decomp_sys or not
    const DecompositionSystem::ERunTests
      ds_test_what = ( ( decomp_sys_testing() == DecompositionSystemHandler_Strategy::DST_TEST
                 || ( decomp_sys_testing() == DecompositionSystemHandler_Strategy::DST_DEFAULT
                  && algo.algo_cntr().check_results() ) )
               ? DecompositionSystem::RUN_TESTS
               : DecompositionSystem::NO_TESTS );
    
    // Determine the output level for decomp_sys				
    DecompositionSystem::EOutputLevel ds_olevel;
    switch(olevel) {
      case PRINT_NOTHING:
      case PRINT_BASIC_ALGORITHM_INFO:
        ds_olevel = DecompositionSystem::PRINT_NONE;
        break;
      case PRINT_ALGORITHM_STEPS:
      case PRINT_ACTIVE_SET:
        ds_olevel = DecompositionSystem::PRINT_BASIC_INFO;
        break;
      case PRINT_VECTORS:
        ds_olevel = DecompositionSystem::PRINT_VECTORS;
        break;
      case PRINT_ITERATION_QUANTITIES:
        ds_olevel = DecompositionSystem::PRINT_EVERY_THING;
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
    };

    // Test the decomposition system
    if( ds_test_what == DecompositionSystem::RUN_TESTS ) {
      // Set the output level
      if( decomp_sys_tester().print_tests() == DecompositionSystemTester::PRINT_NOT_SELECTED ) {
        DecompositionSystemTester::EPrintTestLevel  ds_olevel;
        switch(olevel) {
          case PRINT_NOTHING:
          case PRINT_BASIC_ALGORITHM_INFO:
            ds_olevel = DecompositionSystemTester::PRINT_NONE;
            break;
          case PRINT_ALGORITHM_STEPS:
          case PRINT_ACTIVE_SET:
            ds_olevel = DecompositionSystemTester::PRINT_BASIC;
            break;
          case PRINT_VECTORS:
            ds_olevel = DecompositionSystemTester::PRINT_MORE;
            break;
          case PRINT_ITERATION_QUANTITIES:
            ds_olevel = DecompositionSystemTester::PRINT_ALL;
            break;
          default:
            TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
        }
        decomp_sys_tester().print_tests(ds_olevel);
        decomp_sys_tester().dump_all( olevel == PRINT_ITERATION_QUANTITIES );
      }
      // Run the tests
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
        out << "\nTesting the range/null decompostion ...\n";
      }
      const bool
        decomp_sys_passed = decomp_sys_tester().test_decomp_system(
          s.decomp_sys()
          ,Gc_iq->get_k(0)                   // Gc
          ,Z_iq ? &Z_iq->get_k(0) : NULL     // Z
          ,&Y_iq->get_k(0)                   // Y
          ,&R_iq->get_k(0)                   // R
          ,m > r  ? &Uz_iq->get_k(0) : NULL  // Uz
          ,m > r  ? &Uy_iq->get_k(0) : NULL  // Uy
          ,( olevel >= PRINT_ALGORITHM_STEPS ) ? &out : NULL
          );
      TEUCHOS_TEST_FOR_EXCEPTION(
        !decomp_sys_passed, TestFailed
        ,"EvalNewPointStd_Step::do_step(...) : Error, "
        "the tests of the decomposition system failed!" );
    }
  }
  else {
    // Unconstrained problem
    Z_iq = &s.Z();
    dyn_cast<MatrixSymIdent>(Z_iq->set_k(0)).initialize( nlp.space_x() );
    s.equ_decomp(Range1D::Invalid);
    s.equ_undecomp(Range1D::Invalid);
  }

  // Calculate the rest of the quantities.  If decomp_sys is a variable
  // reduction decomposition system object, then nlp will be hip to the
  // basis selection and will permute these quantities to that basis.
  // Note that x will already be permuted to the current basis.
  if(!Gf_k_updated) { nlp.calc_Gf( x, new_point ); new_point = false; }
  if( m && (!c_k_updated || new_decomp_selected ) ) {
    if(c_k_updated) nlp.set_c( &c_iq->set_k(0) ); // This was not set earlier!
    nlp.calc_c( x, false);
  }
  if(!f_k_updated) {
    nlp.calc_f( x, false);
  }
  nlp.unset_quantities();
  
  // Check for NaN and Inf
  assert_print_nan_inf(f_iq.get_k(0),"f_k",true,&out); 
  if(m)
    assert_print_nan_inf(c_iq->get_k(0),"c_k",true,&out); 
  assert_print_nan_inf(Gf_iq.get_k(0),"Gf_k",true,&out); 
  
  // Print the iteration quantities before we test the derivatives for debugging

  // Update the selection of dependent and independent variables
  if( s.get_decomp_sys().get() ) {
    const ConstrainedOptPack::DecompositionSystemVarReduct
      *decomp_sys_vr = dynamic_cast<ConstrainedOptPack::DecompositionSystemVarReduct*>(&s.decomp_sys());
    if(decomp_sys_vr) {
      var_dep   = decomp_sys_vr->var_dep();
      var_indep = decomp_sys_vr->var_indep();
    }
  }
  
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\nPrinting the updated iteration quantities ...\n";
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out	<< "\nf_k                      = "     << f_iq.get_k(0);
    out << "\n||Gf_k||inf              = "     << Gf_iq.get_k(0).norm_inf();
    if( var_dep.size() )
      out << "\n||Gf_k(var_dep)_k||inf   = " << Gf_iq.get_k(0).sub_view(var_dep)->norm_inf();
    if( var_indep.size() )
      out << "\n||Gf_k(var_indep)_k||inf = " << Gf_iq.get_k(0).sub_view(var_indep)->norm_inf();
    if(m) {
      out << "\n||c_k||inf               = " << c_iq->get_k(0).norm_inf();
      if( calc_matrix_norms && !calc_matrix_info_null_space_only )
        out << "\n||Gc_k||inf              = " << Gc_iq->get_k(0).calc_norm(mat_nrm_inf).value;
      if( n > r && calc_matrix_norms && !calc_matrix_info_null_space_only )
        out << "\n||Z||inf                 = " << Z_iq->get_k(0).calc_norm(mat_nrm_inf).value;
      if( r && calc_matrix_norms && !calc_matrix_info_null_space_only )
        out << "\n||Y||inf                 = " << Y_iq->get_k(0).calc_norm(mat_nrm_inf).value;
      if( r && calc_matrix_norms && !calc_matrix_info_null_space_only  )
        out << "\n||R||inf                 = " << R_iq->get_k(0).calc_norm(mat_nrm_inf).value;
      if( algo.algo_cntr().calc_conditioning() && !calc_matrix_info_null_space_only ) {
        out << "\ncond_inf(R)              = " << R_iq->get_k(0).calc_cond_num(mat_nrm_inf).value;
      }
      if( m > r && calc_matrix_norms && !calc_matrix_info_null_space_only ) {
        out << "\n||Uz_k||inf              = " << Uz_iq->get_k(0).calc_norm(mat_nrm_inf).value;
        out << "\n||Uy_k||inf              = " << Uy_iq->get_k(0).calc_norm(mat_nrm_inf).value;
      }
    }
    out << std::endl;
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
    if(m)
      out << "\nGc_k =\n" << Gc_iq->get_k(0);
    if( n > r )
      out << "\nZ_k =\n" << Z_iq->get_k(0);
    if(r) {
      out << "\nY_k =\n" << Y_iq->get_k(0);
      out << "\nR_k =\n" << R_iq->get_k(0);
    }
    if( m > r ) {
      out << "\nUz_k =\n" << Uz_iq->get_k(0);
      out << "\nUy_k =\n" << Uy_iq->get_k(0);
    }
  }
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    out	<< "\nGf_k =\n" << Gf_iq.get_k(0);
    if( var_dep.size() )
      out << "\nGf(var_dep)_k =\n " << *Gf_iq.get_k(0).sub_view(var_dep);
  }
  if( static_cast<int>(ns_olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    if( var_indep.size() )
      out << "\nGf(var_indep)_k =\n" << *Gf_iq.get_k(0).sub_view(var_indep);
  }
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
    if(m)
      out	<< "\nc_k = \n" << c_iq->get_k(0);
  }
  
  // Check the derivatives if we are checking the results
  if(		fd_deriv_testing() == FD_TEST
    || ( fd_deriv_testing() == FD_DEFAULT && algo.algo_cntr().check_results() )  )
  {
    
    if( olevel >= PRINT_ALGORITHM_STEPS ) {
      out	<< "\n*** Checking derivatives by finite differences\n";
    }

    const bool
      nlp_passed = deriv_tester().finite_diff_check(
        &nlp
        ,x
        ,nb ? &nlp.xl() : NULL
        ,nb ? &nlp.xu() : NULL
        ,m  ? &Gc_iq->get_k(0) : NULL
        ,&Gf_iq.get_k(0)
        ,olevel >= PRINT_VECTORS
        ,( olevel >= PRINT_ALGORITHM_STEPS ) ? &out : NULL
        );
    TEUCHOS_TEST_FOR_EXCEPTION(
      !nlp_passed, TestFailed
      ,"EvalNewPointStd_Step::do_step(...) : Error, "
      "the tests of the nlp derivatives failed!" );
  }

  return true;
}

void EvalNewPointStd_Step::print_step(
   const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
{
  const NLPAlgo       &algo = rsqp_algo(_algo);
  const NLPAlgoState  &s    = algo.rsqp_state();
  const NLP           &nlp  = algo.nlp();
  const size_type
    m = nlp.m();
  out
    << L << "*** Evaluate the new point and update the range/null decomposition\n"
    << L << "if nlp is not initialized then initialize the nlp\n"
    << L << "if x is not updated for any k then set x_k = xinit\n";
  if(m) {
    out
      << L << "if Gc_k is not updated Gc_k = Gc(x_k) <: space_x|space_c\n"
      << L << "For Gc_k = [ Gc_k(:,equ_decomp), Gc_k(:,equ_undecomp) ] where:\n"
      << L << "  Gc_k(:,equ_decomp) <: space_x|space_c(equ_decomp) has full column rank r\n"
      << L << "Find:\n"
      << L << "  Z_k  <: space_x|space_null    s.t. Gc_k(:,equ_decomp)' * Z_k = 0\n"
      << L << "  Y_k  <: space_x|space_range   s.t. [Z_k Y_k] is nonsigular \n"
      << L << "  R_k  <: space_c(equ_decomp)|space_range\n"
      << L << "                                s.t. R_k = Gc_k(:,equ_decomp)' * Y_k\n"
      << L << "  if m > r : Uz_k <: space_c(equ_undecomp)|space_null\n"
      << L << "                                s.t. Uz_k = Gc_k(:,equ_undecomp)' * Z_k\n"
      << L << "  if m > r : Uy_k <: space_c(equ_undecomp)|space_range\n"
      << L << "                                s.t. Uy_k = Gc_k(:,equ_undecomp)' * Y_k\n"
      << L << "begin update decomposition (class \'" << typeName(decomp_sys_handler()) << "\')\n"
      ;
    decomp_sys_handler().print_update_decomposition( algo, s, out, L + "  " );
    out
      << L << "end update decomposition\n"
      << L << "if ( (decomp_sys_testing==DST_TEST)\n"
      << L << "  or (decomp_sys_testing==DST_DEFAULT and check_results==true)\n"
      << L << "  ) then\n"
      << L << "  check properties for Z_k, Y_k, R_k, Uz_k and Uy_k\n"
      << L << "end\n"
      ;
  }
  else {
    out 
      << L << "Z_k = eye(space_x)\n";
  }
  if(m) {
    out
      << L << "end\n";
  }
  out
    << L << "Gf_k = Gf(x_k) <: space_x\n"
    << L << "if m > 0 and c_k is not updated c_k = c(x_k) <: space_c\n"
    << L << "if f_k is not updated f_k = f(x_k) <: REAL\n"
    << L << "if ( (fd_deriv_testing==FD_TEST)\n"
    << L << "  or (fd_deriv_testing==FD_DEFAULT and check_results==true)\n"
    << L << "  ) then\n"
    << L << "  check Gc_k (if m > 0) and Gf_k by finite differences.\n"
    << L << "end\n"
    ;
}

} // end namespace MoochoPack
