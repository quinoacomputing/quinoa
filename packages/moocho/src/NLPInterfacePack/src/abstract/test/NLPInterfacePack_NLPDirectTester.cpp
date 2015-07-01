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

#include <ostream>
#include <iomanip>
#include <sstream>
#include <limits>

#include "NLPInterfacePack_NLPDirectTester.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "TestingHelperPack_update_success.hpp"
#include "Teuchos_Assert.hpp"

namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
} // end namespace

namespace NLPInterfacePack {

NLPDirectTester::NLPDirectTester(
  const calc_fd_prod_ptr_t  &calc_fd_prod
  ,ETestingMethod           Gf_testing_method
  ,ETestingMethod           Gc_testing_method
  ,value_type               Gf_warning_tol
  ,value_type               Gf_error_tol
  ,value_type               Gc_warning_tol
  ,value_type               Gc_error_tol
  ,size_type                num_fd_directions
  ,bool                     dump_all
  )
  :calc_fd_prod_(calc_fd_prod)
  ,Gf_testing_method_(Gf_testing_method)
  ,Gc_testing_method_(Gc_testing_method)
  ,Gf_warning_tol_(Gf_warning_tol)
  ,Gf_error_tol_(Gf_error_tol)
  ,Gc_warning_tol_(Gc_warning_tol)
  ,Gc_error_tol_(Gc_error_tol)
  ,num_fd_directions_(num_fd_directions)
  ,dump_all_(dump_all)
{
  if(calc_fd_prod_ == Teuchos::null)
    calc_fd_prod_ = Teuchos::rcp(new CalcFiniteDiffProd());
}

bool NLPDirectTester::finite_diff_check(
  NLPDirect         *nlp
  ,const Vector     &xo
  ,const Vector     *xl
  ,const Vector     *xu
  ,const Vector     *c
  ,const Vector     *Gf
  ,const Vector     *py
  ,const Vector     *rGf
  ,const MatrixOp   *GcU
  ,const MatrixOp   *D
  ,const MatrixOp   *Uz
  ,bool             print_all_warnings
  ,std::ostream     *out
  ) const
{

  using std::setw;
  using std::endl;
  using std::right;

  using AbstractLinAlgPack::sum;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::random_vector;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::Vp_MtV;
  using LinAlgOpPack::M_StM;
  using LinAlgOpPack::M_StMtM;

  typedef VectorSpace::vec_mut_ptr_t  vec_mut_ptr_t;

//  using AbstractLinAlgPack::TestingPack::CompareDenseVectors;
//  using AbstractLinAlgPack::TestingPack::CompareDenseSparseMatrices;

  using TestingHelperPack::update_success;

  bool success = true, preformed_fd;
  if(out) {
    *out << std::boolalpha
       << std::endl
       << "*********************************************************\n"
       << "*** NLPDirectTester::finite_diff_check(...) ***\n"
       << "*********************************************************\n";
  }

  const Range1D
    var_dep      = nlp->var_dep(),
    var_indep    = nlp->var_indep(),
    con_decomp   = nlp->con_decomp(),
    con_undecomp = nlp->con_undecomp();
  NLP::vec_space_ptr_t
    space_x = nlp->space_x(),
    space_c = nlp->space_c();

  // //////////////////////////////////////////////
  // Validate the input

  TEUCHOS_TEST_FOR_EXCEPTION(
    py && !c, std::invalid_argument
    ,"NLPDirectTester::finite_diff_check(...) : "
    "Error, if py!=NULL then c!=NULL must also be true!" );

  const CalcFiniteDiffProd
    &fd_deriv_prod = this->calc_fd_prod();

  const value_type
    rand_y_l = -1.0, rand_y_u = 1.0,
    small_num = ::sqrt(std::numeric_limits<value_type>::epsilon());

  try {

  // ///////////////////////////////////////////////
  // (1) Check Gf

  if(Gf) {
    switch( Gf_testing_method() ) {
      case FD_COMPUTE_ALL: {
        // Compute FDGf outright
        TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: update above!
        break;
      }
      case FD_DIRECTIONAL: {
        // Compute FDGF'*y using random y's
        if(out)
          *out
            << "\nComparing products Gf'*y with finite difference values FDGf'*y for "
            << "random y's ...\n";
        vec_mut_ptr_t y = space_x->create_member();
        value_type max_warning_viol = 0.0;
        int num_warning_viol = 0;
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
          value_type
            Gf_y = dot( *Gf, *y ),
            FDGf_y;
          preformed_fd = fd_deriv_prod.calc_deriv_product(
            xo,xl,xu
            ,*y,NULL,NULL,true,nlp,&FDGf_y,NULL,out,dump_all(),dump_all()
            );
          if( !preformed_fd )
            goto FD_NOT_PREFORMED;
          assert_print_nan_inf(FDGf_y, "FDGf'*y",true,out);
          const value_type
            calc_err = ::fabs( ( Gf_y - FDGf_y )/( ::fabs(Gf_y) + ::fabs(FDGf_y) + small_num ) );
          if( calc_err >= Gf_warning_tol() ) {
            max_warning_viol = my_max( max_warning_viol, calc_err );
            ++num_warning_viol;
          }
          if(out)
            *out
              << "\nrel_err(Gf'*y,FDGf'*y) = "
              << "rel_err(" << Gf_y << "," << FDGf_y << ") = "
              << calc_err << endl;
          if( calc_err >= Gf_error_tol() ) {
            if(out) {
              *out
                << "Error, above relative error exceeded Gf_error_tol = " << Gf_error_tol() << endl;
              if(dump_all()) {
                *out << "\ny =\n" << *y;
              }
            }
          }
        }
        if(out && num_warning_viol)
          *out
            << "\nThere were " << num_warning_viol << " warning tolerance "
            << "violations out of num_fd_directions = " << num_fd_directions()
            << " computations of FDGf'*y\n"
            << "and the maximum violation was " << max_warning_viol
            << " > Gf_waring_tol = " << Gf_warning_tol() << endl;
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true); // Invalid value
    }
  }

  // /////////////////////////////////////////////
  // (2) Check py = -inv(C)*c
  //
  // We want to check; 
  // 
  //  FDC * (inv(C)*c) \approx c
  //       \_________/
  //         -py
  //
  // We can compute this as:
  //           
  // FDC * py = [ FDC, FDN ] * [ -py ; 0 ]
  //            \__________/
  //                FDA'
  // 
  // t1 =  [ -py ; 0 ]
  // 
  // t2 = FDA'*t1
  // 
  // Compare t2 \approx c
  // 
  if(py) {
    if(out)
      *out
        << "\nComparing c with finite difference product FDA'*[ -py; 0 ] = -FDC*py ...\n";
    // t1 =  [ -py ; 0 ]
    VectorSpace::vec_mut_ptr_t
      t1 = space_x->create_member();
    V_StV( t1->sub_view(var_dep).get(), -1.0, *py );
    *t1->sub_view(var_indep) = 0.0;
    // t2 = FDA'*t1
    VectorSpace::vec_mut_ptr_t
      t2 = nlp->space_c()->create_member();
    preformed_fd = fd_deriv_prod.calc_deriv_product(
      xo,xl,xu
      ,*t1,NULL,NULL,true,nlp,NULL,t2.get(),out,dump_all(),dump_all()
      );
    if( !preformed_fd )
      goto FD_NOT_PREFORMED;
    const value_type
      sum_c  = sum(*c),
      sum_t2 = sum(*t2);
    assert_print_nan_inf(sum_t2, "sum(-FDC*py)",true,out);
    const value_type
      calc_err = ::fabs( ( sum_c - sum_t2 )/( ::fabs(sum_c) + ::fabs(sum_t2) + small_num ) );
    if(out)
      *out
        << "\nrel_err(sum(c),sum(-FDC*py)) = "
        << "rel_err(" << sum_c << "," << sum_t2 << ") = "
        << calc_err << endl;
    if( calc_err >= Gc_error_tol() ) {
      if(out)
        *out
          << "Error, above relative error exceeded Gc_error_tol = " << Gc_error_tol() << endl;
      if(print_all_warnings)
        *out << "\nt1 = [ -py; 0 ] =\n" << *t1
           << "\nt2 = FDA'*t1 = -FDC*py =\n"   << *t2;
      update_success( false, &success );
    }
    if( calc_err >= Gc_warning_tol() ) {
      if(out)
        *out
          << "\nWarning, above relative error exceeded Gc_warning_tol = " << Gc_warning_tol() << endl;
    }
  }

  // /////////////////////////////////////////////
  // (3) Check D = -inv(C)*N

  if(D) {
    switch( Gc_testing_method() ) {
      case FD_COMPUTE_ALL: {
        //
        // Compute FDN outright and check
        // -FDC * D \aprox FDN
        // 
        // We want to compute:
        // 
        // FDC * -D = [ FDC, FDN ] * [ -D; 0 ]
        //            \__________/
        //                FDA'
        // 
        // To compute the above we perform:
        // 
        // T = FDA' * [ -D; 0 ] (one column at a time)
        // 
        // Compare T \approx FDN
        //
/*
        // FDN
        DMatrix FDN(m,n-m);
        fd_deriv_computer.calc_deriv( xo, xl, xu
          , Range1D(m+1,n), nlp, NULL
          , &FDN() ,BLAS_Cpp::trans, out );

        // T = FDA' * [ -D; 0 ] (one column at a time)
        DMatrix T(m,n-m);
        DVector t(n);
        t(m+1,n) = 0.0;
        for( int s = 1; s <= n-m; ++s ) {
          // t = [ -D(:,s); 0 ]
          V_StV( &t(1,m), -1.0, D->col(s) );
          // T(:,s) =  FDA' * t
          fd_deriv_prod.calc_deriv_product(
            xo,xl,xu,t(),NULL,NULL,nlp,NULL,&T.col(s),out);
        }        

        // Compare T \approx FDN
        if(out)
          *out
            << "\nChecking the computed D = -inv(C)*N\n"
            << "where D(i,j) = (-FDC*D)(i,j), dM(i,j) = FDN(i,j) ...\n";
        result = comp_M.comp(
          T(), FDN, BLAS_Cpp::no_trans
          , CompareDenseSparseMatrices::FULL_MATRIX
          , CompareDenseSparseMatrices::REL_ERR_BY_COL
          , Gc_warning_tol(), Gc_error_tol()
          , print_all_warnings, out );
        update_success( result, &success );
        if(!result) return false;
*/
        TEUCHOS_TEST_FOR_EXCEPT(true); // Todo: Implement above!
        break;
      }
      case FD_DIRECTIONAL: {
        //
        // Compute -FDC * D * v \aprox FDN * v
        // for random v's
        //
        // We will compute this as:
        // 
        // t1 = [ 0; y ] <: R^(n)
        // 
        // t2 = FDA' * t1  (  FDN * y ) <: R^(m)
        //
        // t1 = [ -D * y ; 0 ]  <: R^(n)
        // 
        // t3 = FDA' * t1  ( -FDC * D * y ) <: R^(m)
        // 
        // Compare t2 \approx t3
        //
        if(out)
          *out
            << "\nComparing finite difference products -FDC*D*y with FDN*y for "
              "random vectors y ...\n";
        VectorSpace::vec_mut_ptr_t
          y  = space_x->sub_space(var_indep)->create_member(),
          t1 = space_x->create_member(),
          t2 = space_c->create_member(),
          t3 = space_c->create_member();
        value_type max_warning_viol = 0.0;
        int num_warning_viol = 0;
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
          // t1 = [ 0; y ] <: R^(n)
          *t1->sub_view(var_dep)   = 0.0;
          *t1->sub_view(var_indep) = *y;
          // t2 = FDA' * t1  (  FDN * y ) <: R^(m)
          preformed_fd = fd_deriv_prod.calc_deriv_product(
            xo,xl,xu
            ,*t1,NULL,NULL,true,nlp,NULL,t2.get(),out,dump_all(),dump_all()
            );
          if( !preformed_fd )
            goto FD_NOT_PREFORMED;
          // t1 = [ -D * y ; 0 ]  <: R^(n)
          V_StMtV( t1->sub_view(var_dep).get(), -1.0, *D, BLAS_Cpp::no_trans, *y );
          *t1->sub_view(var_indep) = 0.0;
          // t3 = FDA' * t1  ( -FDC * D * y ) <: R^(m)
          preformed_fd = fd_deriv_prod.calc_deriv_product(
            xo,xl,xu
            ,*t1,NULL,NULL,true,nlp,NULL,t3.get(),out,dump_all(),dump_all()
            );
          // Compare t2 \approx t3
          const value_type
            sum_t2 = sum(*t2),
            sum_t3 = sum(*t3);
          const value_type
            calc_err = ::fabs( ( sum_t2 - sum_t3 )/( ::fabs(sum_t2) + ::fabs(sum_t3) + small_num ) );
          if(out)
            *out
              << "\nrel_err(sum(-FDC*D*y),sum(FDN*y)) = "
              << "rel_err(" << sum_t3 << "," << sum_t2 << ") = "
              << calc_err << endl;
          if( calc_err >= Gc_warning_tol() ) {
            max_warning_viol = my_max( max_warning_viol, calc_err );
            ++num_warning_viol;
          }
          if( calc_err >= Gc_error_tol() ) {
            if(out)
              *out
                << "Error, above relative error exceeded Gc_error_tol = " << Gc_error_tol() << endl
                << "Stoping the tests!\n";
            if(print_all_warnings)
              *out << "\ny =\n" << *y
                   << "\nt1 = [ -D*y; 0 ] =\n" << *t1
                   << "\nt2 =  FDA' * [ 0; y ] = FDN * y =\n" << *t2
                   << "\nt3 =  FDA' * t1 = -FDC * D * y =\n" << *t3;
            update_success( false, &success );
          }
        }
        if(out && num_warning_viol)
          *out
            << "\nThere were " << num_warning_viol << " warning tolerance "
            << "violations out of num_fd_directions = " << num_fd_directions()
            << " computations of sum(FDC*D*y) and sum(FDN*y)\n"
            << "and the maximum relative iolation was " << max_warning_viol
            << " > Gc_waring_tol = " << Gc_warning_tol() << endl;
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }

  // ///////////////////////////////////////////////
  // (4) Check rGf = Gf(var_indep) + D'*Gf(var_dep)

  if(rGf) {
    if( Gf && D ) {
      if(out)
        *out
          << "\nComparing rGf_tmp = Gf(var_indep) - D'*Gf(var_dep) with rGf ...\n";
      VectorSpace::vec_mut_ptr_t
        rGf_tmp = space_x->sub_space(var_indep)->create_member();
      *rGf_tmp = *Gf->sub_view(var_indep);
      Vp_MtV( rGf_tmp.get(), *D, BLAS_Cpp::trans, *Gf->sub_view(var_dep) );
      const value_type
        sum_rGf_tmp  = sum(*rGf_tmp),
        sum_rGf      = sum(*rGf);
      const value_type
        calc_err = ::fabs( ( sum_rGf_tmp - sum_rGf )/( ::fabs(sum_rGf_tmp) + ::fabs(sum_rGf) + small_num ) );
      if(out)
        *out
          << "\nrel_err(sum(rGf_tmp),sum(rGf)) = "
          << "rel_err(" << sum_rGf_tmp << "," << sum_rGf << ") = "
          << calc_err << endl;
      if( calc_err >= Gc_error_tol() ) {
        if(out)
          *out
            << "Error, above relative error exceeded Gc_error_tol = " << Gc_error_tol() << endl;
        if(print_all_warnings)
          *out << "\nrGf_tmp =\n" << *rGf_tmp
             << "\nrGf =\n"   << *rGf;
        update_success( false, &success );
      }
      if( calc_err >= Gc_warning_tol() ) {
        if(out)
          *out
            << "\nWarning, above relative error exceeded Gc_warning_tol = "
            << Gc_warning_tol() << endl;
      }
    }
    else if( D ) {
      if(out)
        *out
          << "\nComparing rGf'*y with the finite difference product"
          << " fd_prod(f,[D*y;y]) for random vectors y ...\n";
      VectorSpace::vec_mut_ptr_t
        y  = space_x->sub_space(var_indep)->create_member(),
        t  = space_x->create_member();
      value_type max_warning_viol = 0.0;
      int num_warning_viol = 0;
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
        // t = [ D*y; y ]
        LinAlgOpPack::V_MtV(&*t->sub_view(var_dep),*D,BLAS_Cpp::no_trans,*y);
        *t->sub_view(var_indep) = *y;
        value_type fd_rGf_y = 0.0;
        // fd_Gf_y
        preformed_fd = fd_deriv_prod.calc_deriv_product(
          xo,xl,xu
          ,*t,NULL,NULL,true,nlp,&fd_rGf_y,NULL,out,dump_all(),dump_all()
          );
        if( !preformed_fd )
          goto FD_NOT_PREFORMED;
        if(out) *out << "fd_prod(f,[D*y;y]) = " << fd_rGf_y << "\n";
        // rGf_y = rGf'*y
        const value_type rGf_y = dot(*rGf,*y);
        if(out) *out << "rGf'*y = " << rGf_y << "\n";
        // Compare fd_rGf_y and rGf*y
        const value_type
          calc_err = ::fabs( ( rGf_y - fd_rGf_y )/( ::fabs(rGf_y) + ::fabs(fd_rGf_y) + small_num ) );
        if( calc_err >= Gc_warning_tol() ) {
          max_warning_viol = my_max( max_warning_viol, calc_err );
          ++num_warning_viol;
        }
        if(out)
          *out
            << "\nrel_err(rGf'*y,fd_prod(f,[D*y;y])) = "
            << "rel_err(" << fd_rGf_y << "," << rGf_y << ") = "
            << calc_err << endl;
        if( calc_err >= Gf_error_tol() ) {
          if(out)
            *out << "Error, above relative error exceeded Gc_error_tol = " << Gc_error_tol() << endl;
          if(print_all_warnings)
            *out << "\ny =\n" << *y
                 << "\nt = [ D*y; y ] =\n" << *t;
          update_success( false, &success );
        }
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Test rGf without D? (This is not going to be easy!)
    }
  }
  
  // ///////////////////////////////////////////////////
  // (5) Check GcU, and/or Uz (for undecomposed equalities)

  if(GcU || Uz) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
  }
  
FD_NOT_PREFORMED:

  if(!preformed_fd) {
    if(out)
      *out
        << "\nError, the finite difference computation was not preformed due to cramped bounds\n"
        << "Finite difference test failed!\n" << endl;
    return false;
  }

  } // end try
  catch( const AbstractLinAlgPack::NaNInfException& except ) {
    if(out)
      *out
        << "Error, found a NaN or Inf.  Stoping tests\n";
    success = false;
  }
  
  if( out ) {
    if( success )
      *out
        << "\nCongradulations, all the finite difference errors where within the\n"
        "specified error tolerances!\n";
    else
      *out
        << "\nOh no, at least one of the above finite difference tests failed!\n";
  }

  return success;

}

}  // end namespace NLPInterfacePack
