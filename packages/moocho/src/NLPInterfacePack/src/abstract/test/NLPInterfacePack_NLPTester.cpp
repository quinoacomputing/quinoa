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

#include <iomanip>
#include <ostream>

#include "NLPInterfacePack_NLPTester.hpp"
#include "NLPInterfacePack_NLP.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "TestingHelperPack_update_success.hpp"

namespace NLPInterfacePack {

NLPTester::NLPTester(
  bool     print_all
  ,bool    throw_exception
  )
  :print_all_(print_all), throw_exception_(throw_exception)
{}

bool NLPTester::test_interface(
  NLP                     *nlp
  ,const Vector           &xo
  ,bool                   print_all_warnings
  ,std::ostream           *out
  )
{
  using TestingHelperPack::update_success;
  using AbstractLinAlgPack::assert_print_nan_inf;

  bool result;
  bool success = true;

  if(out) {
    *out << std::boolalpha
       << std::endl
       << "**************************************\n"
       << "*** NLPTester::test_interface(...) ***\n"
       << "**************************************\n";
  }

  try {

    // Initialize the NLP if it has not been already and force in bounds
    if(out)
      *out << "\nnlp->force_xinit_in_bounds(true)";
    nlp->force_xinit_in_bounds();
    if(out)
      *out << "\nnlp->initialize(true)\n";
    nlp->initialize(true);
    
    const size_type
      n = nlp->n(),
      m = nlp->m();
    if(out)
      *out << "\n*** Dimensions of the NLP ...\n"
         << "\nnlp->n()  = " << n
         << "\nnlp->m()  = " << m
         << std::endl;
    if( n < m ) {
      if(*out)
        *out << "Error! n = " << n << " < m = " << m << " is not allowed!\n";
      TEUCHOS_TEST_FOR_EXCEPTION(
        throw_exception_, std::logic_error
        ,"NLPTester::test_interface(...): Error! n = " << n << " < m = " << m << " is not allowed!"
        );
    }

    // Validate the vector spaces
    if(out)
      *out << "\n*** Validate the dimensions of the vector spaces ...\n";
    
    result = nlp->space_x()->dim() == nlp->n();
    update_success( result, &success );
    if(out)
      *out << "\ncheck: nlp->space_x()->dim() = " << nlp->space_x()->dim()
           << " == nlp->n() = " << nlp->n() << ": " << result << std::endl;

    if( nlp->m() ) {
      result = nlp->space_c()->dim() == nlp->m();
      update_success( result, &success );
      if(out)
        *out << "\ncheck: nlp->space_c()->dim() = " << nlp->space_c()->dim()
             << " == nlp->m() = " << nlp->m() << ": " << result << std::endl;
    }
    else {
      result = nlp->space_c().get() == NULL;
      update_success( result, &success );
      if(out)
        *out << "\ncheck: nlp->space_c().get() = " << nlp->space_c().get()
             << " == NULL: " << result << std::endl;
    }

    // Validate the initial guess the bounds on the unknowns.
    if(out)
      *out << "\n*** Validate that the initial starting point is in bounds ...\n";
    const Vector &xinit = nlp->xinit();
    if(out) *out << "\n||nlp->xinit()||inf = " << xinit.norm_inf() << std::endl;
    if(out && print_all()) *out << "\nnlp->xinit() =\n" << xinit;
    assert_print_nan_inf(xinit,"xinit",true,out); 
    const Vector
      &xl = nlp->xl(),
      &xu = nlp->xu();
    if(out && print_all())
      *out << "\nnlp->xl() =\n" << xl
         << "\nnlp->xu() =\n" << xu;
    assert_print_nan_inf(xl,"xl",true,out); 
    assert_print_nan_inf(xu,"xu",true,out); 

    // Validate that xl <= xinit <= xu.
    VectorSpace::vec_mut_ptr_t
      d = nlp->space_x()->create_member();
    *d = 1.0;
    std::pair<value_type,value_type>
      u = AbstractLinAlgPack::max_near_feas_step(
        xinit, *d, nlp->xl(), nlp->xu(), 0.0
        );
    result = u.first >= 0.0;
    update_success( result, &success );
    if(out) {
      *out << "\ncheck: xl <= x <= xu : " << result;
      if(result)
        *out << "\nxinit is in bounds with { max |u| | xl <= x + u <= xu } -> "
           << ( u.first > -u.second ? u.first : u.second  ) << std::endl;
    }
    size_type 
      num_bounded_x = AbstractLinAlgPack::num_bounded(
        nlp->xl(), nlp->xu(), NLP::infinite_bound()
        );
    result = (num_bounded_x == nlp->num_bounded_x());
    update_success( result, &success );
    if(out)
      *out << "\ncheck: num_bounded(nlp->xl(),nlp->xu()) = " << num_bounded_x
         << " == nlp->num_bounded_x() = " << nlp->num_bounded_x()
         << ": " << result << std::endl;

    // Get the initial Lagrange multipliers
    if(out)
      *out << "\nGetting the initial estimates for the Lagrange mutipliers ...\n";
    VectorSpace::vec_mut_ptr_t  lambda, nu;
    nlp->get_init_lagrange_mult(
      (  nlp->m()
         ? (lambda  = nlp->space_c()->create_member()).get() 
         : (VectorMutable*)NULL )
      ,( nlp->num_bounded_x()
         ? (nu = nlp->space_x()->create_member()).get()
         : (VectorMutable*)NULL )
      );

    if(out) {
      if(lambda.get())
        *out << "\n||lambda||inf  = " << lambda->norm_inf();
      if(nu.get())
        *out << "\n||nu||inf      = " << nu->norm_inf()
           << "\nnu.nz()        = " << nu->nz();
      *out << std::endl;
      if(print_all()) {
        if(lambda.get())
          *out << "\nlambda =\n" << *lambda;
        if(nu.get())
          *out << "\nnu =\n" << *nu;
      }
    }
    if(lambda.get())
      assert_print_nan_inf(*lambda,"lambda",true,out); 
    if(nu.get())
      assert_print_nan_inf(*nu,"nu",true,out); 

    // Save the current reference that are set to be set back at the end
    value_type      *f_saved = NULL;
    VectorMutable   *c_saved = NULL;
    f_saved = nlp->get_f();
    if( nlp->m() )  c_saved = nlp->get_c();

    // Create calcualtion quantities
    value_type                   f;
    VectorSpace::vec_mut_ptr_t   c;
    if( nlp->m() )
      c = nlp->space_c()->create_member();

    // Set the calculation quantities
    nlp->set_f(&f);
    if( nlp->m() )  nlp->set_c(c.get());

    // Calculate the quantities at xo

    if(out)
      *out << "\n*** Evaluate the point xo ...\n";

    if(out)	*out << "\n||xo||inf = " << xo.norm_inf() << std::endl;
    if(out && print_all()) *out << "\nxo =\n" << xo;
    assert_print_nan_inf(xo,"xo",true,out); 

    nlp->calc_f(xo,true);
    if(nlp->m())  nlp->calc_c(xo,false);

    if(out) {
      *out << "\nf(xo) = " << f;
      if(nlp->m())
        *out << "\n||c(xo)||inf = " << nlp->c().norm_inf();
      *out << std::endl;
      if(print_all()) {
        if(nlp->m())
          *out << "\nc(xo) =\n" << nlp->c();
      }
    }

    if(c.get())
      assert_print_nan_inf(*c,"c(xo)",true,out); 

    // Report the final solution!
    if(out)
      *out << "\n*** Report this point to the NLP as suboptimal ...\n";
    nlp->report_final_solution(	xo, lambda.get(), nu.get(), false );

    // Print the number of evaluations!
    if(out) {
      *out << "\n*** Print the number of evaluations ...\n";
      *out << "\nnlp->num_f_evals() = " << nlp->num_f_evals();
      if(nlp->m())
        *out << "\nnlp->num_c_evals() = " << nlp->num_c_evals();
      *out << std::endl;
    }

    // Set the original quantities back
    nlp->set_f(f_saved);
    if(nlp->m())  nlp->set_c(c_saved);

  }
  catch(const std::exception& except) {
    if(out)
      *out << "Caught a std::exception: " << except.what() << std::endl;
    success = false;
    if(throw_exception())
      throw;
  }
  catch(...) {
    if(out)
      *out << "Caught an unknown exception!\n";
    success = false;
    if(throw_exception())
      throw;
  }

  return success;
}

} // namespace NLPInterfacePack
