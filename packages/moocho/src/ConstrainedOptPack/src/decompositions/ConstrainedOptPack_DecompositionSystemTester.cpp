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

#include <math.h>

#include <limits>
#include <ostream>

#include "ConstrainedOptPack_DecompositionSystemTester.hpp"
#include "ConstrainedOptPack_DecompositionSystem.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsingTester.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_MatrixComposite.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace ConstrainedOptPack {

DecompositionSystemTester::DecompositionSystemTester(
  EPrintTestLevel  print_tests
  ,bool            dump_all
  ,bool            throw_exception
  ,size_type       num_random_tests
  ,value_type      mult_warning_tol
  ,value_type      mult_error_tol
  ,value_type      solve_warning_tol
  ,value_type      solve_error_tol
  )
  :print_tests_(print_tests)
  ,dump_all_(dump_all)
  ,throw_exception_(throw_exception)
  ,num_random_tests_(num_random_tests)
  ,mult_warning_tol_(mult_warning_tol)
  ,mult_error_tol_(mult_error_tol)
  ,solve_warning_tol_(solve_warning_tol)
  ,solve_error_tol_(solve_error_tol)
{}
 
bool DecompositionSystemTester::test_decomp_system(
  const DecompositionSystem   &ds
  ,const MatrixOp             &Gc
  ,const MatrixOp             *Z
  ,const MatrixOp             *Y
  ,const MatrixOpNonsing      *R
  ,const MatrixOp             *Uz
  ,const MatrixOp             *Uy
  ,std::ostream               *out
  )
{
  namespace rcp = MemMngPack;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using AbstractLinAlgPack::sum;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::Vp_StMtV;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using AbstractLinAlgPack::random_vector;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_VpV;
  using LinAlgOpPack::Vp_V;

  bool success = true, result, lresult, llresult;
  const value_type
    rand_y_l  = -1.0,
    rand_y_u  = 1.0,
    small_num = ::pow(std::numeric_limits<value_type>::epsilon(),0.25),
    alpha     = 2.0,
    beta      = 3.0;

  EPrintTestLevel
    print_tests = ( this->print_tests() == PRINT_NOT_SELECTED ? PRINT_NONE : this->print_tests() );

  // Print the input?
  if( out && print_tests != PRINT_NONE ) {
    if( print_tests >= PRINT_BASIC )
      *out << "\n**********************************************************"
         << "\n*** DecompositionSystemTester::test_decomp_system(...) ***"
         << "\n**********************************************************\n";
  }

  const size_type
    n = ds.n(),
    m = ds.m(),
    r = ds.r();
  const Range1D
    equ_decomp       = ds.equ_decomp(),
    equ_undecomp     = ds.equ_undecomp();

  // print dimensions, ranges
  if( out && print_tests >= PRINT_MORE ) {
    *out
      << "\nds.n()                  = " << n
      << "\nds.m()                  = " << m
      << "\nds.r()                  = " << r
      << "\nds.equ_decomp()         = ["<<equ_decomp.lbound()<<","<<equ_decomp.ubound()<<"]"
      << "\nds.equ_undecomp()       = ["<<equ_undecomp.lbound()<<","<<equ_undecomp.ubound()<<"]"
      << "\nds.space_range()->dim() = " << ds.space_range()->dim()
      << "\nds.space_null()->dim()  = " << ds.space_null()->dim()
      << std::endl;
  }

  // validate input matrices
  TEUCHOS_TEST_FOR_EXCEPTION(
      Z==NULL&&Y==NULL&&R==NULL&&Uz==NULL&&Uy==NULL
    , std::invalid_argument
    ,"DecompositionSystemTester::test_decomp_system(...) : Error, "
    "at least one of Z, Y, R, Uz or Uy can not be NULL!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    m == r && Uz != NULL, std::invalid_argument
    ,"DecompositionSystemTester::test_decomp_system(...) : Error, "
    "Uz must be NULL if m==r is NULL!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    m == r && Uy != NULL, std::invalid_argument
    ,"DecompositionSystemTester::test_decomp_system(...) : Error, "
    "Uy must be NULL if m==r is NULL!" );

  // Print the input?
  if( out && print_tests != PRINT_NONE ) {
    if(dump_all()) {
      *out << "\nGc =\n"       << Gc;
      if(Z)
        *out << "\nZ =\n"    << *Z;
      if(Y)
        *out << "\nY =\n"    << *Y;
      if(R)
        *out << "\nR =\n"    << *R;
      if(Uz)
        *out << "\nUz =\n"   << *Uz;
      if(Uy)
        *out << "\nUy =\n"   << *Uy;
    }
  }

  //
  // Check the dimensions of everything
  //

  if( out && print_tests >= PRINT_BASIC )
    *out << "\n1) Check the partitioning ranges and vector space dimensions ...";
  lresult = true;

  if( out && print_tests >= PRINT_MORE )
    *out << "\n\n1.a) check: equ_decomp.size() + equ_undecomp.size() == ds.m() : ";
  result = equ_decomp.size() + equ_undecomp.size() == ds.m();
  if(out && print_tests >= PRINT_MORE)
    *out << ( result ? "passed" : "failed" );
  if(!result) lresult = false;

  if( out && print_tests >= PRINT_MORE )
    *out << "\n\n1.b) check: equ_decomp.size() == ds.r() : ";
  result = equ_decomp.size() == ds.r();
  if(out && print_tests >= PRINT_MORE)
    *out << ( result ? "passed" : "failed" );
  if(!result) lresult = false;

  if( out && print_tests >= PRINT_MORE )
    *out << "\n\n1.c) check: ds.space_range()->dim() == ds.r() : ";
  result = ds.space_range()->dim() == ds.r();
  if(out && print_tests >= PRINT_MORE)
    *out << ( result ? "passed" : "failed" );
  if(!result) lresult = false;

  if( out && print_tests >= PRINT_MORE )
    *out << "\n\n1.d) check: ds.space_null()->dim() == ds.n() - ds.r() : ";
  result = ds.space_null()->dim() == ds.n() - ds.r();
  if(out && print_tests >= PRINT_MORE)
    *out << ( result ? "passed" : "failed" );
  if(!result) lresult = false;

  if(out && print_tests >= PRINT_MORE)
    *out << std::endl;

  if(!lresult) success = false;
  if( out && print_tests == PRINT_BASIC )
    *out << " : " << ( lresult ? "passed" : "failed" );

  //
  // Perform the tests
  //

  if(out && print_tests >= PRINT_BASIC)
    *out
      << "\n2) Check the compatibility of the vector spaces for Gc, Z, Y, R, Uz and Uy  ...";
  lresult = true;
  
  if(Z) {
    if(out && print_tests >= PRINT_MORE)
      *out
        << "\n2.a) Check consistency of the vector spaces for:"
        << "\n    Z.space_cols() == Gc.space_cols() and Z.space_rows() == ds.space_null()";
    llresult = true;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.a.1) Z->space_cols().is_compatible(Gc.space_cols()) == true : ";
    result = Z->space_cols().is_compatible(Gc.space_cols());	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.a.2) Z->space_cols().is_compatible(*ds.space_null()) == true : ";
    result = Z->space_rows().is_compatible(*ds.space_null());	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" );
  }

  if(Y) {
    if(out && print_tests >= PRINT_MORE)
      *out
        << "\n2.b) Check consistency of the vector spaces for:"
        << "\n    Y.space_cols() == Gc.space_cols() and Y.space_rows() == ds.space_range()";
    llresult = true;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.b.1) Y->space_cols().is_compatible(Gc.space_cols()) == true : ";
    result = Y->space_cols().is_compatible(Gc.space_cols());	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.b.2) Y->space_cols().is_compatible(*ds.space_range()) == true : ";
    result = Y->space_rows().is_compatible(*ds.space_range());	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" );
  }

  if(R) {
    if(out && print_tests >= PRINT_MORE)
      *out
        << "\n2.c) Check consistency of the vector spaces for:"
        << "\n    R.space_cols() == Gc.space_cols()(equ_decomp) and R.space_rows() == ds.space_range()";
    llresult = true;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.c.1) R->space_cols().is_compatible(*Gc.space_cols().sub_space(equ_decomp)) == true : ";
    result = R->space_cols().is_compatible(*Gc.space_cols().sub_space(equ_decomp));	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.c.2) R->space_cols().is_compatible(*ds.space_range()) == true : ";
    result = R->space_rows().is_compatible(*ds.space_range());	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" );
  }

  if(Uz) {
    if(out && print_tests >= PRINT_MORE)
      *out
        << "\n2.d) Check consistency of the vector spaces for:"
        << "\n    Uz.space_cols() == Gc.space_cols()(equ_undecomp) and Uz.space_rows() == ds.space_null()";
    llresult = true;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.d.1) Uz->space_cols().is_compatible(*Gc.space_cols().sub_space(equ_undecomp)) == true : ";
    result = Uz->space_cols().is_compatible(*Gc.space_cols().sub_space(equ_undecomp));	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.d.2) Uz->space_cols().is_compatible(*ds.space_null()) == true : ";
    result = Uz->space_rows().is_compatible(*ds.space_null());	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" );
  }

  if(Uy) {
    if(out && print_tests >= PRINT_MORE)
      *out
        << "\n2.e) Check consistency of the vector spaces for:"
        << "\n    Uy.space_cols() == Gc.space_cols()(equ_undecomp) and Uy.space_rows() == ds.space_range()";
    llresult = true;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.e.1) Uy->space_cols().is_compatible(*Gc.space_cols().sub_space(equ_undecomp)) == true : ";
    result = Uy->space_cols().is_compatible(*Gc.space_cols().sub_space(equ_undecomp));	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(out && print_tests >= PRINT_ALL)
      *out << "\n\n2.e.2) Uy->space_cols().is_compatible(*ds.space_range()) == true : ";
    result = Uy->space_rows().is_compatible(*ds.space_range());	
    if(out && print_tests >= PRINT_ALL)
      *out << ( result ? "passed" : "failed" )
         << std::endl;
    if(!result) llresult = false;
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" );
  }

  if(!lresult) success = false;
  if( out && print_tests == PRINT_BASIC )
    *out << " : " << ( lresult ? "passed" : "failed" );

  if(out && print_tests >= PRINT_BASIC)
    *out
      << "\n3) Check the compatibility of the matrices Gc, Z, Y, R, Uz and Uy numerically ...";
  
  if(Z) {

    if(out && print_tests >= PRINT_MORE)
      *out
        << std::endl
        << "\n3.a) Check consistency of:"
        << "\n     op ( alpha* Gc(:,equ_decomp)' * beta*Z ) * v"
        << "\n         \\__________________________________/"
        << "\n                         A"
        << "\n    ==  op( alpha*beta*Uz * v"
        << "\n            \\___________/"
        << "\n                     B"
        << "\nfor random vectors v ...";

    VectorSpace::vec_mut_ptr_t
      v_c       = Gc.space_rows().create_member(),
      v_c_tmp   = v_c->space().create_member(),
      v_x       = Gc.space_cols().create_member(),
      v_x_tmp   = v_x->space().create_member(),
      v_z       = ds.space_null()->create_member(),
      v_z_tmp   = v_z->space().create_member();
    
    if(out && print_tests >= PRINT_MORE)
      *out << "\n\n3.a.1) Testing non-transposed A*v == B*v ...";
    if(out && print_tests > PRINT_MORE)
      *out << std::endl;
    llresult = true;
    {for( int k = 1; k <= num_random_tests(); ++k ) {
      random_vector( rand_y_l, rand_y_u, v_z.get() );
      if(out && print_tests >= PRINT_ALL) {
        *out
          << "\n3.a.1."<<k<<") random vector " << k << " ( ||v_z||_1 / n = " << (v_z->norm_1() / v_z->dim()) << " )\n";
        if(dump_all() && print_tests >= PRINT_ALL)
          *out << "\nv_z =\n" << *v_z;
      }
      V_StMtV( v_x.get(), beta, *Z, no_trans, *v_z );
      V_StMtV( v_c.get(), alpha, Gc, trans, *v_x );
      *v_c_tmp->sub_view(equ_decomp) = 0.0;
      if(equ_undecomp.size()) {
        if(Uz)
          V_StMtV( v_c_tmp->sub_view(equ_undecomp).get(), alpha*beta, *Uz, no_trans, *v_z );
        else
          *v_c_tmp->sub_view(equ_undecomp).get() = *v_c->sub_view(equ_undecomp);
      }
      const value_type
        sum_Bv  = sum(*v_c_tmp), // should be zero if equ_undecomp.size() == 0 so scale by 1.0
        sum_Av  = sum(*v_c);
      assert_print_nan_inf(sum_Bv, "sum(B*v_z)",true,out);
      assert_print_nan_inf(sum_Av, "sum(A*v_z)",true,out);
      const value_type
        calc_err = ::fabs( ( sum_Av - sum_Bv )
                   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + (equ_undecomp.size() ? small_num : 1.0) ) );
      if(out && print_tests >= PRINT_ALL)
        *out
          << "\nrel_err(sum(A*v_z),sum(B*v_z)) = "
          << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
          << calc_err << std::endl;
      if( calc_err >= mult_warning_tol() ) {
        if(out && print_tests >= PRINT_ALL)
          *out
            << std::endl
            << ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
            << ", rel_err(sum(A*v_z),sum(B*v_z)) = "
            << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
            << calc_err
            << " exceeded "
            << ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
            << " = "
            << ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
            << std::endl;
        if(calc_err >= mult_error_tol()) {
          if(dump_all() && print_tests >= PRINT_ALL) {
            *out << "\nalpha = " << alpha << std::endl;
            *out << "\nbeta  = " << beta  << std::endl;
            *out << "\nv_z =\n"           << *v_z;
            *out << "\nbeta*Z*v_z =\n"    << *v_x;
            *out << "\nA*v_z =\n"         << *v_c;
            *out << "\nB*v_z =\n"         << *v_c_tmp;
          }
          llresult = false;
        }
      }
    }}
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" )
         << std::endl;
    
    if(out && print_tests >= PRINT_MORE)
      *out << "\n\n3.a.2) Testing transposed A'*v == B'*v ...";
    if(out && print_tests > PRINT_MORE)
      *out << std::endl;
    llresult = true;
    {for( int k = 1; k <= num_random_tests(); ++k ) {
      random_vector( rand_y_l, rand_y_u, v_c.get() );
      if(out && print_tests >= PRINT_ALL) {
        *out
          << "\n3.a.2."<<k<<") random vector " << k << " ( ||v_c||_1 / n = " << (v_c->norm_1() / v_c->dim()) << " )\n";
        if(dump_all() && print_tests >= PRINT_ALL)
          *out << "\nv_c =\n" << *v_c;
      }
      V_StMtV( v_x.get(), alpha, Gc, no_trans, *v_c );
      V_StMtV( v_z.get(), beta,  *Z, trans,    *v_x );
      *v_z_tmp = 0.0;
      if(equ_undecomp.size()) {
        if(Uz)
          V_StMtV( v_z_tmp.get(), alpha*beta, *Uz, trans, *v_c->sub_view(equ_undecomp) );
        else
          *v_z_tmp = *v_z;
      }
      const value_type
        sum_Bv  = sum(*v_z_tmp), // should be zero so scale by 1.0
        sum_Av  = sum(*v_z);
      assert_print_nan_inf(sum_Bv, "sum(B'*v_c)",true,out);
      assert_print_nan_inf(sum_Av, "sum(A'*v_c)",true,out);
      const value_type
        calc_err = ::fabs( ( sum_Av - sum_Bv )
                   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + (equ_undecomp.size() ? small_num : 1.0) ) );
      if(out && print_tests >= PRINT_ALL)
        *out
          << "\nrel_err(sum(A'*v_c),sum(B'*v_c)) = "
          << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
          << calc_err << std::endl;
      if( calc_err >= mult_warning_tol() ) {
        if(out && print_tests >= PRINT_ALL)
          *out
            << std::endl
            << ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
            << ", rel_err(sum(A'*v_c),sum(B'*v_c)) = "
            << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
            << calc_err
            << " exceeded "
            << ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
            << " = "
            << ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
            << std::endl;
        if(calc_err >= mult_error_tol()) {
          if(dump_all() && print_tests >= PRINT_ALL) {
            *out << "\nalpha = " << alpha << std::endl;
            *out << "\nbeta  = " << beta  << std::endl;
            *out << "\nv_c =\n"           << *v_c;
            *out << "\nalpha*Gc*v_c =\n"  << *v_x;
            *out << "\nA'*v_c =\n"        << *v_z;
            *out << "\nB'*v_c =\n"        << *v_z_tmp;
          }
          llresult = false;
        }
      }
    }}
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" )
         << std::endl;

  }
  else {
    if(out && print_tests >= PRINT_MORE)
      *out
        << std::endl
        << "\n3.a) Warning! Z ==NULL; Z, and Uz are not checked numerically ...\n";
  }

  if(Y) {

    if(out && print_tests >= PRINT_MORE)
      *out
        << std::endl
        << "\n3.b) Check consistency of:"
        << "\n     op ( alpha*[ Gc(:,equ_decomp)'   ]"
        << "\n                [ Gc(:,equ_undecomp)' ] * beta*Y ) * v"
        << "\n         \\_____________________________________/"
        << "\n                         A"
        << "\n    ==  op( alpha*beta*[ R  ]"
        << "\n                       [ Uy ] ) * v"
        << "\n            \\_______________/"
        << "\n                     B"
        << "\nfor random vectors v ...";

    VectorSpace::vec_mut_ptr_t
      v_c       = Gc.space_rows().create_member(),
      v_c_tmp   = v_c->space().create_member(),
      v_x       = Gc.space_cols().create_member(),
      v_x_tmp   = v_x->space().create_member(),
      v_y       = ds.space_range()->create_member(),
      v_y_tmp   = v_y->space().create_member();
    
    if(out && print_tests >= PRINT_MORE)
      *out << "\n\n3.b.1) Testing non-transposed A*v == B*v ...";
    if(out && print_tests > PRINT_MORE)
      *out << std::endl;
    llresult = true;
    {for( int k = 1; k <= num_random_tests(); ++k ) {
      random_vector( rand_y_l, rand_y_u, v_y.get() );
      if(out && print_tests >= PRINT_ALL) {
        *out
          << "\n3.b.1."<<k<<") random vector " << k << " ( ||v_y||_1 / n = " << (v_y->norm_1() / v_y->dim()) << " )\n";
        if(dump_all() && print_tests >= PRINT_ALL)
          *out << "\nv_y =\n" << *v_y;
      }
      V_StMtV( v_x.get(), beta, *Y, no_trans, *v_y );
      V_StMtV( v_c.get(), alpha, Gc, trans, *v_x );
      V_StMtV( v_c_tmp->sub_view(equ_decomp).get(), alpha*beta, *R, no_trans, *v_y );
      if(equ_undecomp.size()) {
        if(Uy)
          V_StMtV( v_c_tmp->sub_view(equ_undecomp).get(), alpha*beta, *Uy, no_trans, *v_y );
        else
          *v_c_tmp->sub_view(equ_undecomp) = *v_c->sub_view(equ_undecomp);
      }
      const value_type
        sum_Bv  = sum(*v_c_tmp),
        sum_Av  = sum(*v_c);
      assert_print_nan_inf(sum_Bv, "sum(B*v_y)",true,out);
      assert_print_nan_inf(sum_Av, "sum(A*v_y)",true,out);
      const value_type
        calc_err = ::fabs( ( sum_Av - sum_Bv )
                   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + small_num ) );
      if(out && print_tests >= PRINT_ALL)
        *out
          << "\nrel_err(sum(A*v_y),sum(B*v_y)) = "
          << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
          << calc_err << std::endl;
      if( calc_err >= mult_warning_tol() ) {
        if(out && print_tests >= PRINT_ALL)
          *out
            << std::endl
            << ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
            << ", rel_err(sum(A*v_y),sum(B*v_y)) = "
            << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
            << calc_err
            << " exceeded "
            << ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
            << " = "
            << ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
            << std::endl;
        if(calc_err >= mult_error_tol()) {
          if(dump_all() && print_tests >= PRINT_ALL) {
            *out << "\nalpha = " << alpha << std::endl;
            *out << "\nbeta  = " << beta  << std::endl;
            *out << "\nv_y =\n"           << *v_y;
            *out << "\nbeta*Y*v_y =\n"    << *v_x;
            *out << "\nA*v_y =\n"         << *v_c;
            *out << "\nB*v_y =\n"         << *v_c_tmp;
          }
          llresult = false;
        }
      }
    }}
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" )
         << std::endl;
    
    if(out && print_tests >= PRINT_MORE)
      *out << "\n\n3.b.2) Testing transposed A'*v == B'*v ...";
    if(out && print_tests > PRINT_MORE)
      *out << std::endl;
    llresult = true;
    {for( int k = 1; k <= num_random_tests(); ++k ) {
      random_vector( rand_y_l, rand_y_u, v_c.get() );
      if(out && print_tests >= PRINT_ALL) {
        *out
          << "\n3.a.2."<<k<<") random vector " << k << " ( ||v_c||_1 / n = " << (v_c->norm_1() / v_c->dim()) << " )\n";
        if(dump_all() && print_tests >= PRINT_ALL)
          *out << "\nv_c =\n" << *v_c;
      }
      V_StMtV( v_x.get(), alpha, Gc, no_trans, *v_c );
      V_StMtV( v_y.get(), beta,  *Y, trans,    *v_x );
      V_StMtV( v_y_tmp.get(), alpha*beta, *R, trans, *v_c->sub_view(equ_decomp) );
      if(equ_undecomp.size()) {
        if(Uy)
          Vp_StMtV( v_y_tmp.get(), alpha*beta, *Uy, trans, *v_c->sub_view(equ_undecomp) );
        else
          Vp_V( v_y_tmp.get(), *v_y );
      }
      const value_type
        sum_Bv  = sum(*v_y_tmp), // should be zero so scale by 1.0
        sum_Av  = sum(*v_y);
      assert_print_nan_inf(sum_Bv, "sum(B'*v_c)",true,out);
      assert_print_nan_inf(sum_Av, "sum(A'*v_c)",true,out);
      const value_type
        calc_err = ::fabs( ( sum_Av - sum_Bv )
                   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + small_num ) );
      if(out && print_tests >= PRINT_ALL)
        *out
          << "\nrel_err(sum(A'*v_c),sum(B'*v_c)) = "
          << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
          << calc_err << std::endl;
      if( calc_err >= mult_warning_tol() ) {
        if(out && print_tests >= PRINT_ALL)
          *out
            << std::endl
            << ( calc_err >= mult_error_tol() ? "Error" : "Warning" )
            << ", rel_err(sum(A'*v_c),sum(B'*v_c)) = "
            << "rel_err(" << sum_Av << "," << sum_Bv << ") = "
            << calc_err
            << " exceeded "
            << ( calc_err >= mult_error_tol() ? "mult_error_tol" : "mult_warning_tol" )
            << " = "
            << ( calc_err >= mult_error_tol() ? mult_error_tol() : mult_warning_tol() )
            << std::endl;
        if(calc_err >= mult_error_tol()) {
          if(dump_all() && print_tests >= PRINT_ALL) {
            *out << "\nalpha = " << alpha << std::endl;
            *out << "\nbeta  = " << beta  << std::endl;
            *out << "\nv_c =\n"           << *v_c;
            *out << "\nalpha*Gc*v_c =\n"  << *v_x;
            *out << "\nA'*v_c =\n"        << *v_y;
            *out << "\nB'*v_c =\n"        << *v_y_tmp;
          }
          llresult = false;
        }
      }
    }}
    if(!llresult) lresult = false;
    if( out && print_tests == PRINT_MORE )
      *out << " : " << ( llresult ? "passed" : "failed" )
         << std::endl;

  }
  else {
    if(out && print_tests >= PRINT_MORE)
      *out
        << std::endl
        << "\n3.b) Warning! Y ==NULL; Y, R and Uy are not checked numerically ...\n";
  }

  if(R) {
    if(out && print_tests >= PRINT_MORE)
      *out
        << std::endl
        << "\n3.b) Check consistency of: op(op(inv(R))*op(R)) == I ...\n";
    typedef MatrixOpNonsingTester  MWONST_t;
    MWONST_t::EPrintTestLevel
      olevel;
    switch(print_tests) {
      case PRINT_NONE:
      case PRINT_BASIC:
        olevel = MWONST_t::PRINT_NONE;
        break;
      case PRINT_MORE:
        olevel = MWONST_t::PRINT_MORE;
        break;
      case PRINT_ALL:
        olevel = MWONST_t::PRINT_ALL;
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here
    }
    MWONST_t
      R_tester(
        MWONST_t::TEST_LEVEL_2_BLAS
        ,olevel
        ,dump_all()
        ,throw_exception()
        ,num_random_tests()
        ,solve_warning_tol()
        ,solve_error_tol()
        );
    lresult = R_tester.test_matrix(*R,"R",out);
  }

  if(!lresult) success = false;
  if( out && print_tests == PRINT_BASIC )
    *out << " : " << ( lresult ? "passed" : "failed" );
  
  if( out && print_tests != PRINT_NONE ) {
    if(success)
      *out << "\nCongradulations! The DecompositionSystem object and its associated matrix objects seem to check out!\n";
    else
      *out << "\nOops! At least one of the tests did not check out!\n";
    if( print_tests >= PRINT_BASIC )
      *out << "\nEnd DecompositionSystemTester::test_decomp_system(...)\n";
  }

  return success;
}

} // end namespace ConstrainedOptPack
