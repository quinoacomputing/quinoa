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
#include <math.h>

#include <iomanip>
#include <sstream>
#include <limits>

#include "NLPInterfacePack_NLPFirstDerivTester.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "TestingHelperPack_update_success.hpp"
#include "Teuchos_Assert.hpp"

namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
} // end namespace

namespace NLPInterfacePack {

NLPFirstDerivTester::NLPFirstDerivTester(
  const calc_fd_prod_ptr_t  &calc_fd_prod
  ,ETestingMethod           fd_testing_method
  ,size_type                num_fd_directions
  ,value_type               warning_tol
  ,value_type               error_tol
  )
  :calc_fd_prod_(calc_fd_prod)
  ,fd_testing_method_(fd_testing_method)
  ,num_fd_directions_(num_fd_directions)
  ,warning_tol_(warning_tol)
  ,error_tol_(error_tol)
{}

bool NLPFirstDerivTester::finite_diff_check(
  NLP               *nlp
  ,const Vector     &xo
  ,const Vector     *xl
  ,const Vector     *xu
  ,const MatrixOp   *Gc
  ,const Vector     *Gf
  ,bool             print_all_warnings
  ,std::ostream     *out
  ) const
{
  namespace rcp = MemMngPack;
  using AbstractLinAlgPack::assert_print_nan_inf;

  const size_type
    //n  = nlp->n(),
    m  = nlp->m();

  // ///////////////////////////////////
  // Validate the input

  TEUCHOS_TEST_FOR_EXCEPTION(
    !m && Gc, std::invalid_argument
    ,"NLPFirstDerivTester::finite_diff_check(...) : "
    "Error, Gc must be NULL if m == 0" );

  assert_print_nan_inf(xo, "xo",true,out);
  if(Gf)
    assert_print_nan_inf(*Gf, "Gf",true,out); 

  bool success = true;

  try {

  switch(fd_testing_method_) {
    case FD_COMPUTE_ALL:
      return fd_check_all(nlp,xo,xl,xu,Gc,Gf,print_all_warnings,out);
    case FD_DIRECTIONAL:
      return fd_directional_check(nlp,xo,xl,xu,Gc,Gf,print_all_warnings,out);
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  } // end try
  catch( const AbstractLinAlgPack::NaNInfException& except ) {
    if(out)
      *out
        << "Error: found a NaN or Inf.  Stoping tests!\n";
    success = false;
  }
  
  return success;	// will never be executed
}

// private

bool NLPFirstDerivTester::fd_check_all(
  NLP               *nlp
  ,const Vector     &xo
  ,const Vector     *xl
  ,const Vector     *xu
  ,const MatrixOp   *Gc
  ,const Vector     *Gf
  ,bool             print_all_warnings
  ,std::ostream     *out
  ) const
{
  using std::setw;
  using std::endl;
  using std::right;

  namespace rcp = MemMngPack;
  using AbstractLinAlgPack::sum;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using AbstractLinAlgPack::random_vector;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_MtV;

   using TestingHelperPack::update_success;

  bool success = true;

  const size_type
    n  = nlp->n();
  //m  = nlp->m();

  // //////////////////////////////////////////////
  // Validate the input

  NLP::vec_space_ptr_t
    space_x = nlp->space_x(),
    space_c = nlp->space_c();

  const CalcFiniteDiffProd
    &fd_deriv_prod = this->calc_fd_prod();

  const value_type
    //rand_y_l = -1.0, rand_y_u = 1.0,
    small_num = ::pow(std::numeric_limits<value_type>::epsilon(),0.25);

  if(out)
    *out
      << "\nComparing Gf and/or Gc with finite-difference values "
        "FDGf and/or FDGc one variable i at a time ...\n";

  value_type  max_Gf_warning_viol = 0.0;
  int         num_Gf_warning_viol = 0;
  value_type  max_Gc_warning_viol = 0.0;
  int         num_Gc_warning_viol = 0;

  VectorSpace::vec_mut_ptr_t
    e_i       = space_x->create_member(),
    Gc_i      = ( Gc ? space_c->create_member()  : Teuchos::null ),
    FDGc_i    = ( Gc ? space_c->create_member()  : Teuchos::null );
  *e_i = 0.0;

  for( size_type i = 1; i <= n; ++ i ) {
    EtaVector eta_i(i,n);
    e_i->set_ele(i,1.0);
    // Compute exact??? values
    value_type
      Gf_i = Gf ? Gf->get_ele(i) : 0.0;
    if(Gc)
      V_MtV( Gc_i.get(), *Gc, BLAS_Cpp::trans, eta_i() );
    // Compute finite difference values
    value_type
      FDGf_i;
    const bool preformed_fd = fd_deriv_prod.calc_deriv_product(
      xo,xl,xu
      ,*e_i
      ,NULL // fo
      ,NULL // co
      ,true // check_nan_inf
      ,nlp
      ,Gf ? &FDGf_i : NULL
      ,Gc ? FDGc_i.get() : NULL
      ,out
      );
    if( !preformed_fd ) {
      if(out)
        *out
          << "\nError, the finite difference computation was not preformed due to cramped bounds\n"
          << "Finite difference test failed!\n" << endl;
      return false;
    }
    
    // Compare the quantities
    // Gf
    assert_print_nan_inf(FDGf_i, "FDGf'*e_i",true,out);
    const value_type
      Gf_err = ::fabs( Gf_i - FDGf_i ) / ( ::fabs(Gf_i) + ::fabs(FDGf_i) + small_num );
    if(out)
      *out
        << "\nrel_err(Gf("<<i<<"),FDGf'*e("<<i<<")) = "
        << "rel_err(" << Gf_i << "," << FDGf_i << ") = "
        << Gf_err << endl;
    if( Gf_err >= warning_tol() ) {
      max_Gf_warning_viol = my_max( max_Gf_warning_viol, Gf_err );
      ++num_Gf_warning_viol;
    }
    if( Gf_err >= error_tol() ) {
      if(out)
        *out
          << "\nError, exceeded Gf_error_tol = " << error_tol() << endl
          << "Stoping the tests!\n";
      if(print_all_warnings)
        *out
          << "\ne_i =\n"     << *e_i
          << "\nGf_i =\n"    << Gf_i << std::endl
          << "\nFDGf_i =\n"  << FDGf_i << std::endl;
      update_success( false, &success );
      return false;
    }
    // Gc
    if(Gc) {
      const value_type
        sum_Gc_i   = sum(*Gc_i),
        sum_FDGc_i = sum(*FDGc_i);
      assert_print_nan_inf(sum_FDGc_i, "sum(FDGc'*e_i)",true,out);
      const value_type
        calc_err = ::fabs( ( sum_Gc_i - sum_FDGc_i )
                   /( ::fabs(sum_Gc_i) + ::fabs(sum_FDGc_i) + small_num ) );
      if(out)
        *out
          << "\nrel_err(sum(Gc'*e("<<i<<")),sum(FDGc'*e("<<i<<"))) = "
          << "rel_err(" << sum_Gc_i << "," << sum_FDGc_i << ") = "
          << calc_err << endl;
      if( calc_err >= warning_tol() ) {
        max_Gc_warning_viol = my_max( max_Gc_warning_viol, calc_err );
        ++num_Gc_warning_viol;
      }
      if( calc_err >= error_tol() ) {
        if(out)
          *out
            << "\nError, rel_err(sum(Gc'*e("<<i<<")),sum(FDGc'*e("<<i<<"))) = "
            << "rel_err(" << sum_Gc_i << "," << sum_FDGc_i << ") = "
            << calc_err << endl
            << "exceeded error_tol = " << error_tol() << endl
            << "Stoping the tests!\n";
        if(print_all_warnings)
          *out << "\ne_i =\n"     << *e_i
             << "\nGc_i =\n"    << *Gc_i
             << "\nFDGc_i =\n"  << *FDGc_i;
        update_success( false, &success );
        return false;
      }
      if( calc_err >= warning_tol() ) {
        if(out)
          *out
            << "\nWarning, rel_err(sum(Gc'*e("<<i<<")),sum(FDGc'*e("<<i<<"))) = "
            << "rel_err(" << sum_Gc_i << "," << sum_FDGc_i << ") = "
            << calc_err << endl
            << "exceeded warning_tol = " << warning_tol() << endl;
      }
    }
    e_i->set_ele(i,0.0);
  }
  if(out && num_Gf_warning_viol)
    *out
      << "\nFor Gf, there were " << num_Gf_warning_viol << " warning tolerance "
      << "violations out of num_fd_directions = " << num_fd_directions()
      << " computations of FDGf'*e(i)\n"
      << "and the maximum violation was " << max_Gf_warning_viol
      << " > Gf_waring_tol = " << warning_tol() << endl;
  if(out && num_Gc_warning_viol)
    *out
      << "\nFor Gc, there were " << num_Gc_warning_viol << " warning tolerance "
      << "violations out of num_fd_directions = " << num_fd_directions()
      << " computations of FDGc'*e(i)\n"
      << "and the maximum violation was " << max_Gc_warning_viol
      << " > Gc_waring_tol = " << warning_tol() << endl;
  if(out)
    *out
      << "\nCongradulations!  All of the computed errors were within the specified error tolerance!\n";

  return true;

}

bool NLPFirstDerivTester::fd_directional_check(
  NLP               *nlp
  ,const Vector     &xo
  ,const Vector     *xl
  ,const Vector     *xu
  ,const MatrixOp   *Gc
  ,const Vector     *Gf
  ,bool             print_all_warnings
  ,std::ostream     *out
  ) const
{
  using std::setw;
  using std::endl;
  using std::right;

  namespace rcp = MemMngPack;
  using AbstractLinAlgPack::sum;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using AbstractLinAlgPack::random_vector;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_MtV;

   using TestingHelperPack::update_success;

  bool success = true;

  //const size_type
  //n  = nlp->n(),
  //m  = nlp->m();

  // //////////////////////////////////////////////
  // Validate the input

  NLP::vec_space_ptr_t
    space_x = nlp->space_x(),
    space_c = nlp->space_c();

  const CalcFiniteDiffProd
    &fd_deriv_prod = this->calc_fd_prod();

  const value_type
    rand_y_l = -1.0, rand_y_u = 1.0,
    small_num = ::pow(std::numeric_limits<value_type>::epsilon(),0.25);

  if(out)
    *out
      << "\nComparing directional products Gf'*y and/or Gc'*y with finite difference values "
        " FDGf'*y and/or FDGc'*y for random y's ...\n";

  value_type  max_Gf_warning_viol = 0.0;
  int         num_Gf_warning_viol = 0;

  VectorSpace::vec_mut_ptr_t
    y         = space_x->create_member(),
    Gc_prod   = ( Gc ? space_c->create_member()  : Teuchos::null ),
    FDGc_prod = ( Gc ? space_c->create_member()  : Teuchos::null );

  const int num_fd_directions_used = ( num_fd_directions() > 0 ? num_fd_directions() : 1 );

  for( int direc_i = 1; direc_i <= num_fd_directions_used; ++direc_i ) {
    if( num_fd_directions() > 0 ) {
      random_vector( rand_y_l, rand_y_u, y.get() );
      if(out)
        *out
          << "\n****"
          << "\n**** Random directional vector " << direc_i << " ( ||y||_1 / n = "
          << (y->norm_1() / y->dim()) << " )"
          << "\n***\n";
    }
    else {
      *y = 1.0;
      if(out)
        *out
          << "\n****"
          << "\n**** Ones vector y ( ||y||_1 / n = "<<(y->norm_1()/y->dim())<<" )"
          << "\n***\n";
    }
    // Compute exact??? values
    value_type
      Gf_y = Gf ? dot( *Gf, *y ) : 0.0;
    if(Gc)
      V_MtV( Gc_prod.get(), *Gc, BLAS_Cpp::trans, *y );
    // Compute finite difference values
    value_type
      FDGf_y;
    const bool preformed_fd = fd_deriv_prod.calc_deriv_product(
      xo,xl,xu
      ,*y
      ,NULL // fo
      ,NULL // co
      ,true // check_nan_inf
      ,nlp
      ,Gf ? &FDGf_y : NULL
      ,Gc ? FDGc_prod.get() : NULL
      ,out
      );
    if( !preformed_fd ) {
      if(out)
        *out
          << "\nError, the finite difference computation was not preformed due to cramped bounds\n"
          << "Finite difference test failed!\n" << endl;
      return false;
    }
    
    // Compare the quantities
    // Gf
    assert_print_nan_inf(FDGf_y, "FDGf'*y",true,out);
    const value_type
      Gf_err = ::fabs( Gf_y - FDGf_y ) / ( ::fabs(Gf_y) + ::fabs(FDGf_y) + small_num );
    if(out)
      *out
        << "\nrel_err(Gf'*y,FDGf'*y) = "
        << "rel_err(" << Gf_y << "," << FDGf_y << ") = "
        << Gf_err << endl;
    if( Gf_err >= warning_tol() ) {
      max_Gf_warning_viol = my_max( max_Gf_warning_viol, Gf_err );
      ++num_Gf_warning_viol;
    }
    if( Gf_err >= error_tol() ) {
      if(out)
        *out
          << "\nError, exceeded Gf_error_tol = " << error_tol() << endl
          << "Stoping the tests!\n";
      return false;
    }
    // Gc
    if(Gc) {
      const value_type
        sum_Gc_prod   = sum(*Gc_prod),
        sum_FDGc_prod = sum(*FDGc_prod);
      assert_print_nan_inf(sum_FDGc_prod, "sum(FDGc'*y)",true,out);
      const value_type
        calc_err = ::fabs( ( sum_Gc_prod - sum_FDGc_prod )
                   /( ::fabs(sum_Gc_prod) + ::fabs(sum_FDGc_prod) + small_num ) );
      if(out)
        *out
          << "\nrel_err(sum(Gc'*y),sum(FDGc'*y)) = "
          << "rel_err(" << sum_Gc_prod << "," << sum_FDGc_prod << ") = "
          << calc_err << endl;
      if( calc_err >= error_tol() ) {
        if(out)
          *out
            << "\nError, rel_err(sum(Gc'*y),sum(FDGc'*y)) = "
            << "rel_err(" << sum_Gc_prod << "," << sum_FDGc_prod << ") = "
            << calc_err << endl
            << "exceeded error_tol = " << error_tol() << endl
            << "Stoping the tests!\n";
        if(print_all_warnings)
          *out << "\ny =\n"          << *y
             << "\nGc_prod =\n"    << *Gc_prod
             << "\nFDGc_prod =\n"  << *FDGc_prod;
        update_success( false, &success );
        return false;
      }
      if( calc_err >= warning_tol() ) {
        if(out)
          *out
            << "\nWarning, rel_err(sum(Gc'*y),sum(FDGc'*y)) = "
            << "rel_err(" << sum_Gc_prod << "," << sum_FDGc_prod << ") = "
            << calc_err << endl
            << "exceeded warning_tol = " << warning_tol() << endl;
      }
    }
  }
  if(out && num_Gf_warning_viol)
    *out
      << "\nFor Gf, there were " << num_Gf_warning_viol << " warning tolerance "
      << "violations out of num_fd_directions = " << num_fd_directions()
      << " computations of FDGf'*y\n"
      << "and the maximum violation was " << max_Gf_warning_viol
      << " > Gf_waring_tol = " << warning_tol() << endl;

  if(out)
    *out
      << "\nCongradulations!  All of the computed errors were within the specified error tolerance!\n";

  return true;
}

}	// end namespace NLPInterfacePack
