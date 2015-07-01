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

#include <assert.h>

#include "NLPInterfacePack_test_nlp_direct.hpp"
#include "NLPInterfacePack_CalcFiniteDiffProd.hpp"
#include "NLPInterfacePack_CalcFiniteDiffProdSetOptions.hpp"
#include "NLPInterfacePack_NLPTester.hpp"
#include "NLPInterfacePack_NLPTesterSetOptions.hpp"
#include "NLPInterfacePack_NLPDirectTester.hpp"
#include "NLPInterfacePack_NLPDirectTesterSetOptions.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorSpaceTester.hpp"
#include "AbstractLinAlgPack_VectorSpaceTesterSetOptions.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "TestingHelperPack_update_success.hpp"

bool NLPInterfacePack::test_nlp_direct(
  NLPDirect                                    *nlp
  ,OptionsFromStreamPack::OptionsFromStream    *options
  ,std::ostream                                *out
  )
{
  using TestingHelperPack::update_success;

  bool result;
  bool success = true;
  
  Teuchos::VerboseObjectTempState<NLP>
    nlpOutputTempState(Teuchos::rcp(nlp,false),Teuchos::getFancyOStream(Teuchos::rcp(out,false)),Teuchos::VERB_LOW);

  if(out)
    *out
      << "\n****************************"
      << "\n*** test_nlp_direct(...) ***"
      << "\n*****************************\n";
  
  nlp->initialize(true);

  // Test the DVector spaces
  if(out)
    *out << "\nTesting the vector spaces ...\n";
  
  VectorSpaceTester vec_space_tester;
  if(options) {
    VectorSpaceTesterSetOptions
      opt_setter(&vec_space_tester);
    opt_setter.set_options(*options);
  }

  if(out)
    *out << "\nTesting nlp->space_x() ...\n";
  result = vec_space_tester.check_vector_space(*nlp->space_x(),out);
  if(out) {
    if(result)
      *out << "nlp->space_x() checks out!\n";
    else
      *out << "nlp->space_x() check failed!\n";
  }
  update_success( result, &success );
  if(out)
    *out << "\nTesting nlp->space_x()->sub_space(nlp->var_dep()) ...\n";
  result = vec_space_tester.check_vector_space(
    *nlp->space_x()->sub_space(nlp->var_dep()),out);
  if(out) {
    if(result)
      *out << "nlp->space_x()->sub_space(nlp->var_dep()) checks out!\n";
    else
      *out << "nlp->space_x()->sub_space(nlp->var_dep()) check failed!\n";
  }
  update_success( result, &success );
  if(out)
    *out << "\nTesting nlp->space_x()->sub_space(nlp->var_indep()) ...\n";
  result = vec_space_tester.check_vector_space(
    *nlp->space_x()->sub_space(nlp->var_indep()),out);
  if(out) {
    if(result)
      *out << "nlp->space_x()->sub_space(nlp->var_indep()) checks out!\n";
    else
      *out << "nlp->space_x()->sub_space(nlp->var_indep()) check failed!\n";
  }
  update_success( result, &success );
  if(out)
    *out << "\nTesting nlp->space_c() ...\n";
  result = vec_space_tester.check_vector_space(*nlp->space_c(),out);
  if(out) {
    if(result)
      *out << "nlp->space_c() checks out!\n";
    else
      *out << "nlp->space_c() check failed!\n";
  }
  update_success( result, &success );
  if(out)
    *out << "\nTesting nlp->space_c()->sub_space(nlp->con_decomp()) ...\n";
  result = vec_space_tester.check_vector_space(
    *nlp->space_c()->sub_space(nlp->con_decomp()),out);
  if(out) {
    if(result)
      *out << "nlp->space_c()->sub_space(nlp->con_decomp()) checks out!\n";
    else
      *out << "nlp->space_c()->sub_space(nlp->con_decomp()) check failed!\n";
  }
  update_success( result, &success );
  if( nlp->con_decomp().size() < nlp->m() ) {
    if(out)
      *out << "\nTesting nlp->space_c()->sub_space(nlp->con_undecomp()) ...\n";
    result = vec_space_tester.check_vector_space(
      *nlp->space_c()->sub_space(nlp->con_undecomp()),out);
    if(out) {
      if(result)
        *out << "nlp->space_c()->sub_space(nlp->con_undecomp()) checks out!\n";
      else
        *out << "nlp->space_c()->sub_space(nlp->con_undecomp()) check failed!\n";
    }
    update_success( result, &success );
  }

  // Test the NLP interface first!

  NLPTester nlp_tester;
  if(options) {
    NLPTesterSetOptions
      nlp_tester_opt_setter(&nlp_tester);
    nlp_tester_opt_setter.set_options(*options);
  }
  const bool print_all_warnings = nlp_tester.print_all();

  result = nlp_tester.test_interface(
    nlp, nlp->xinit(), print_all_warnings, out );
  update_success( result, &success );
  
  // Test the NLPDirect interface now!

  const bool supports_Gf = nlp->supports_Gf();

  if(out)
    *out << "\nCalling nlp->calc_point(...) at nlp->xinit() ...\n";
  const size_type
    n  = nlp->n(),
    m  = nlp->m();
  const Range1D
    var_dep      = nlp->var_dep(),
    var_indep    = nlp->var_indep(),
    con_decomp   = nlp->con_decomp(),
    con_undecomp = nlp->con_undecomp();
  VectorSpace::vec_mut_ptr_t
    c   =      nlp->space_c()->create_member(),
    Gf  =      nlp->space_x()->create_member(),
    py  =      nlp->space_x()->sub_space(var_dep)->create_member(),
    rGf =      nlp->space_x()->sub_space(var_indep)->create_member();
  NLPDirect::mat_fcty_ptr_t::element_type::obj_ptr_t
    GcU = con_decomp.size() < m ? nlp->factory_GcU()->create() : Teuchos::null,
    D   =                         nlp->factory_D()->create(),
    Uz   = con_decomp.size() < m ? nlp->factory_Uz()->create()   : Teuchos::null;
  nlp->calc_point(
    nlp->xinit(), NULL, c.get(), true, supports_Gf?Gf.get():NULL, py.get(), rGf.get()
    ,GcU.get(), D.get(), Uz.get() );
  if(out) {
    if(supports_Gf) {
      *out << "\n||Gf||inf  = " << Gf->norm_inf();
      if(nlp_tester.print_all())
        *out << "\nGf =\n" << *Gf;
    }
    *out << "\n||py||inf  = " << py->norm_inf();
    if(nlp_tester.print_all())
      *out << "\npy =\n" << *py;
    *out << "\n||rGf||inf  = " << rGf->norm_inf();
    if(nlp_tester.print_all())
      *out << "\nrGf =\n" << *rGf;
    if(nlp_tester.print_all())
      *out << "\nD =\n" << *D;
    if( con_decomp.size() < m ) {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Print GcU and Uz
    }
    *out << "\n";
  }

  CalcFiniteDiffProd
    calc_fd_prod;
  if(options) {
    CalcFiniteDiffProdSetOptions
      options_setter( &calc_fd_prod );
    options_setter.set_options(*options);
  }
  NLPDirectTester
    nlp_first_order_direct_tester(Teuchos::rcp(&calc_fd_prod,false));
  if(options) {
    NLPDirectTesterSetOptions
      nlp_tester_opt_setter(&nlp_first_order_direct_tester);
    nlp_tester_opt_setter.set_options(*options);
  }
  result = nlp_first_order_direct_tester.finite_diff_check(
    nlp, nlp->xinit()
    ,nlp->num_bounded_x() ? &nlp->xl() : NULL
    ,nlp->num_bounded_x() ? &nlp->xu() : NULL
    ,c.get()
    ,supports_Gf?Gf.get():NULL,py.get(),rGf.get(),GcU.get(),D.get(),Uz.get()
    ,print_all_warnings, out
    );
  update_success( result, &success );

  return success;
}
