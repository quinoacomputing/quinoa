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

#include <ostream>

#include "AbstractLinAlgPack_VectorSpaceTester.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "TestingHelperPack_update_success.hpp"
#include "Teuchos_Assert.hpp"


//#define ALAP_VECTOR_SPACE_TESTER_DUMP

#ifdef ALAP_VECTOR_SPACE_TESTER_DUMP
#  include "RTOpPack_SPMD_apply_op.hpp"
#endif // ALAP_VECTOR_SPACE_TESTER_DUMP


namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
} // end namespace


namespace AbstractLinAlgPack {


VectorSpaceTester::VectorSpaceTester(
  bool         print_all_tests
  ,bool        print_vectors
  ,bool        throw_exception
  ,size_type   num_random_tests
  ,value_type  warning_tol
  ,value_type  error_tol
  )
  :print_all_tests_(print_all_tests)
  ,print_vectors_(print_vectors)
  ,throw_exception_(throw_exception)
  ,num_random_tests_(num_random_tests)
  ,warning_tol_(warning_tol)
  ,error_tol_(error_tol)
{}

bool VectorSpaceTester::check_vector_space(
  const VectorSpace &space,
  std::ostream *out
  ) const
{

  using Teuchos::as;

  bool success = true, result = false;

  try {

  // Adapted from test_my_vector_library(...)

  const value_type 
    rand_l = -10.0,
    rand_u = +10.0;

  typedef VectorMutable::vec_ptr_t       vec_ptr_t;
  typedef VectorMutable::vec_mut_ptr_t   vec_mut_ptr_t;

  // Create three random non-mutable vectors
  vec_ptr_t            v_array[3];
  const Vector*  v[3];
  {for(int k = 0; k < 3; ++k) {
    vec_mut_ptr_t  r = space.create_member();
    random_vector( rand_l, rand_u, r.get() );
    v_array[k] = r;
    v[k] = v_array[k].get();
  }}

  // Create six random mutable vectors
  vec_mut_ptr_t          z_array[6];
  VectorMutable*   z[6];
  {for(int k = 0; k < 6; ++k) {
    random_vector( rand_l, rand_u, (z_array[k]= space.create_member()).get() );
    z[k] = z_array[k].get();
  }}

  if(out && print_all_tests())
    *out << std::boolalpha
       << "\n**************************************************"
       << "\n*** VectorSpaceTester::check_vector_space(...) ***"
       << "\n**************************************************\n";

  index_type
    n = space.dim();
  RTOp_value_type
    err     = 0.0,
    max_err = 0.0,
    sum_err = 0.0;
  char z_name[20], v_name[20];

  // Print the vector space dimension
  if(out && print_all_tests())
    *out << "\nspace->dim() = " << n << std::endl;

  // Print the initial vectors
  if(out && print_vectors()) {
    *out << "\n*** Printing the immutable vectors\n";
    {for(int j = 0; j < 3; ++j) {
      sprintf( v_name, "v[%d]", j );
      *out << std::endl << v_name << " : " << typeName(*v[j]) << std::endl
         << *v[j];
    }}
    *out << "\n*** Printing the mutable vectors\n";
    {for(int k = 0; k < 6; ++k) {
      sprintf( z_name, "z[%d]", k );
      *out << std::endl << z_name << " : " << typeName(*z[k]) << std::endl
         << *z[k];
    }}
  }

  ///////////////////////////////////////////////////////////////////////////
  if(out && print_all_tests())
    *out << "\n*** Testing the obvious assertions\n";
  {
    {for(int k = 0; k < 6; ++k) {
      const Vector *vec = NULL;
      if( k < 3 ) {
        sprintf( v_name, "v[%d]", k );
        vec = v[k];
      }
      else {
        sprintf( v_name, "z[%d]", k-3 );
        vec = z[k-3];
      }
      
      result = vec->space().is_compatible(space);
      if(out && (print_all_tests() || !result))
        *out << "\ncheck: " << v_name << "->space().is_compatible(space) : " << result << std::endl;
      check_test( result ? 0.0 : -10.0 , out, &success );
      
      result = space.is_compatible(vec->space());
      if(out && (print_all_tests() || !result))
        *out << "check: space.is_compatible(" << v_name << "->space()) : " << result << std::endl;
      check_test( result ? 0.0 : -10.0 , out, &success );
      
      err = vec->dim() - space.dim();
      if(out && (print_all_tests() || !result))
        *out << "check: " << v_name << "->dim() - space.dim() = " << vec->dim() << " - "
           << space.dim() << " = " << err << std::endl;
      check_test( err , out, &success );

      result = vec->nz() <= space.dim();
      if(out && (print_all_tests() || !result))
        *out << "check: " << v_name << "->nz() <= space.dim() = " << vec->nz() << " <= " << space.dim()
           << " : " << result << std::endl;
      check_test( result ? 0.0 : -10.0 , out, &success );

    }}
  }

  //////////////////////////////////////////////////////////////////////////
  if(out && print_all_tests())
    *out << "\n*** Testing scalar assignment and element access methods\n";
  {
    const index_type k = 0;
    sprintf( z_name, "z[%d]", (int)k );
    if(out && print_all_tests())
      *out << "\n0.0 -> " << z_name << std::endl;
    *z[k] = 0.0;
    if(out && print_vectors())
      *out << std::endl << z_name << " =\n" << *z[k];
    {for(index_type r = 0; r < num_random_tests(); ++r) {
      std::srand( n / (1+r) + r ); // This is very important in a parallel program!
      const index_type
        i = my_min(
          n,
          my_max(
            as<index_type>(((double)std::rand() / RAND_MAX) * n + 1.0),
            as<index_type>(1)
            )
          );
      const RTOp_value_type
        val = 10.0;
      
      if(out && print_all_tests())
        *out << std::endl << z_name << ".set_ele("<<i<<","<<val<<")\n";
#ifdef ALAP_VECTOR_SPACE_TESTER_DUMP
      RTOpPack::show_spmd_apply_op_dump = true;
#endif
      z[k]->set_ele(i,val);
#ifdef ALAP_VECTOR_SPACE_TESTER_DUMP
      RTOpPack::show_spmd_apply_op_dump = false;
#endif
      if(out && print_vectors())
        *out << std::endl << z_name << " =\n" << *z[k];
#ifdef ALAP_VECTOR_SPACE_TESTER_DUMP
      RTOpPack::show_spmd_apply_op_dump = true;
#endif
      RTOp_value_type
        val_get = z[k]->get_ele(i);
#ifdef ALAP_VECTOR_SPACE_TESTER_DUMP
      RTOpPack::show_spmd_apply_op_dump = false;
#endif
      err = val - val_get;
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
        *out << "check: " << val << " - " << z_name << ".get_ele(" << i << ") = "
           << val << " - " << val_get << " = " << err << std::endl;
      check_test(err,out,&success);
      
      RTOp_value_type
        z_k_sum = sum(*z[k]);
      err = val - z_k_sum;
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
        *out << "check: " << val << " - sum(" << z_name << ") = "
           << val << " - " << z_k_sum << " = " << err << std::endl;
      check_test(err,out,&success);

      if(out && print_all_tests())
        *out << z_name << ".set_ele("<<i<<",0.0)\n";
      z[k]->set_ele(i,0.0);
      if(out && print_vectors())
        *out << std::endl << z_name << " =\n" << *z[k];
      z_k_sum = sum(*z[k]);
      err = z_k_sum;
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
        *out << "check: sum(" << z_name << ") = " << z_k_sum << std::endl;
      check_test(err,out,&success);
    }}
  }

  //////////////////////////////////////////////////////////////////
  if(out && print_all_tests())
    *out << "\n*** Testing vector assignment\n";
  {
    {for( int k = 0; k < 3; ++k ) {
      sprintf( z_name, "z[%d]", k );
      sprintf( v_name, "v[%d]", k );
      if(out && print_all_tests())
        *out << "\n" << v_name << " -> " << z_name << "\n";
      *z[k] = *v[k];
      const RTOp_value_type
        sum_z_k = sum(*z[k]),
        sum_v_k = sum(*v[k]);
      err = (sum_z_k - sum_v_k)/n;
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
        *out << "check: (sum(" << z_name << ") - sum(" << v_name << "))/n = ("
           << sum_z_k << " - " << sum_v_k << ")/" << n << " = " << err << std::endl;
      check_test(err,out,&success);
    }}
  }

  /*

  //////////////////////////////////////////////////////////////////
  if(out && print_all_tests())
    *out << "\n*** Testing sub-vector and sub-space access\n";
  {
    const index_type k = 0;
    sprintf( z_name, "z[%d]", k );
    if(out && print_all_tests())
      *out << "\n0.0 -> " << z_name << std::endl;
    *z[k] = 0.0;
    if(out && print_vectors())
      *out << std::endl << z_name << " =\n" << *z[k];
    {for(int r = 0; r < num_random_tests(); ++r) {
      index_type
        i1 = my_min( n, (index_type)(((double)rand() / RAND_MAX) * n + 1) ),
        i2 = my_min( n, (index_type)(((double)rand() / RAND_MAX) * n + 1) );
      if( i1 > i2 ) std::swap( i1, i2 );
      index_type
        sub_vec_dim = i2-i1+1;
      const RTOp_value_type
        val = 10.0;
      
      if(out && print_all_tests())
        *out << "\nsub_vec_mut = " << z_name
           << ".sub_view(" << i1 << "," << i2 << ")\n";
      VectorMutable::vec_mut_ptr_t
        sub_vec_mut = z[k]->sub_view(i1,i2);
      if(out && print_all_tests())
        *out << "sub_vec_mut = " << val << std::endl;
      *sub_vec_mut = val;
      if(out && print_vectors())
        *out << std::endl << z_name << " =\n" << *z[k];

      err = sub_vec_dim - sub_vec_mut->dim();
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()))
        *out << "check: ("<<i2<<"-"<<i1<<"+1) - sub_vec_mut.dim() = "
           << sub_vec_dim <<" - " << sub_vec_mut->dim() << " = " << err << std::endl;
      check_test(err,out,&success);

      if(out && print_all_tests())
        *out << "\nsub_space = space.sub_space(" << i1 << "," << i2 << ")\n";
      VectorSpace::space_ptr_t
        sub_space = space.sub_space(i1,i2);

      result = sub_vec_mut->space().is_compatible(*sub_space);
      if(out && (print_all_tests() || !result))
        *out << "check: sub_vec_mut->space().is_compatible(*sub_space) : " << result << std::endl;
      check_test( result ? 0.0 : -10.0 , out, &success );
      
      result = sub_space->is_compatible(sub_vec_mut->space());
      if(out && (print_all_tests() || !result))
        *out << "check: sub_space->is_compatible(*sub_vec_mut->space()) : " << result << std::endl;
      check_test( result ? 0.0 : -10.0 , out, &success );
        
      RTOp_value_type
        expected_sum = val*sub_vec_dim,
        z_k_sum = sum( *z[k] );
      err = (expected_sum - z_k_sum)/sub_vec_mut->dim();
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()))
        *out << "check: ("<<val<<"*("<<i2<<"-"<<i1<<"+1) - sum("<<z_name<<"))/"<<sub_vec_dim
           << " = ("<<expected_sum<<" - "<<z_k_sum<<")/"<<sub_vec_dim<<" = "<<err<<std::endl;
      check_test(err,out,&success);
         
      if(out && print_all_tests())
        *out << "sub_vec = "<<z_name<<"{const}.sub_view("<<i1<<","<<i2<<")\n";
      Vector::vec_ptr_t
        sub_vec = static_cast<const Vector*>(z[k])->sub_view(i1,i2);

      err = sub_vec_dim - sub_vec->dim();
      if(out && print_all_tests())
        *out << "check: ("<<i2<<"-"<<i1<<"+1) - sub_vec.dim() = "
           << sub_vec_dim << " - " << sub_vec->dim() << " = " << err << std::endl;
      check_test(err,out,&success);
        
      expected_sum = val*sub_vec_dim;
      z_k_sum = sum(*sub_vec);
      err = (expected_sum - z_k_sum)/sub_vec_mut->dim();
      if(out && print_all_tests())
        *out << "check: ("<<val<<"*("<<i2<<"-"<<i1<<"+1) - sum(sub_vec))/"<<sub_vec_dim
           << " = ("<<expected_sum<<" - "<<z_k_sum<<")/"<<sub_vec_dim<<" = "<<err<<std::endl;

      if(out && print_all_tests())
        *out << "sub_vec_mut = 0.0\n";
      *sub_vec_mut = 0.0;
      if(out && print_vectors())
        *out << std::endl << z_name << " =\n" << *z[k];
    }}
  }

  //////////////////////////////////////////////////////////////////
  if(out && print_all_tests())
    *out << "\n*** Testing explicit sub-vector access\n";
  {
    const index_type k = 0;
    const Vector
      &v_from = *v[k];
    sprintf( v_name, "v[%d]", k );
    VectorMutable
      &z_to = *z[k];
    sprintf( z_name, "z[%d]", k );
    
    if(out && print_all_tests())
      *out << "\n0.0 -> " << z_name << std::endl;
    *z[k] = 0.0;
    if(out && print_vectors())
      *out << std::endl << z_name << " =\n" << *z[k];

    {for(int r = 0; r < num_random_tests(); ++r) {
      const index_type // Get random small sub-vectors so parallel efficiency will be good
        i1 = my_min( n, (index_type)(((double)rand() / RAND_MAX) * n + 1) ),
        i2 = my_min( (index_type)(i1 + ((double)rand() / RAND_MAX) * 9), n );
      const index_type
        sub_vec_dim = i2-i1+1;

      if(out && print_all_tests())
        *out << std::endl << v_name << ".get_sub_vector(Rang1D("<<i1<<","<<i2<<"),SPARSE,&sub_vec)\n";
      RTOpPack::SubVector sub_vec;
      v_from.get_sub_vector(Range1D(i1,i2),&sub_vec);	

      err = sub_vec_dim - sub_vec.subDim();
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()))
        *out << "check: ("<<i2<<"-"<<i1<<"+1) - sub_vec.subDim() = "
           << sub_vec_dim << " - " << sub_vec.subDim() << "  = " << err << std::endl;
      check_test(err,out,&success);

      if(out && print_all_tests())
        *out << z_name << ".set_sub_vector(sub_vec)\n";
      RTOpPack::SparseSubVector spc_sub_vec( sub_vec );
      z_to.set_sub_vector(spc_sub_vec);
      if(out && print_vectors())
        *out << std::endl << z_name << " =\n" << z_to;
      
      const RTOp_value_type
        v_sub_vec_sum = sum(*v_from.sub_view(i1,i2)),
        z_sum         = sum(z_to);
      err = (v_sub_vec_sum - z_sum)/sub_vec_dim;
      if(out && (print_all_tests() || ::fabs(err) >= warning_tol()))
        *out << "check: (sum(*"<<v_name<<".sub_view("<<i1<<","<<i2<<"))-sum("<<z_name<<"))/sub_vec.subDim()"
          " = ("<<v_sub_vec_sum<<"-"<<z_sum<<")/"<<sub_vec_dim<<" = "<<err<<std::endl;

      if(out && print_all_tests())
        *out << v_name << ".free_sub_vector(&sub_vec)\n";
      v_from.free_sub_vector(&sub_vec);	
        
      if(out && print_all_tests())
        *out << "*" << z_name<<".sub_view("<<i1<<","<<i2<<") = 0.0\n";
      *z_to.sub_view(i1,i2) = 0.0;
      if(out && print_vectors())
        *out << std::endl << z_name << " =\n" << z_to;

    }}
  }

  */

  //////////////////////////////////////////////////////////////////
  if(out && print_all_tests())
    *out << "\n*** Testing norms\n";
  if(n > 1) {
    const index_type k = 0;
    sprintf( z_name, "z[%d]", (int)k );

    const value_type val1 = -2.0, val2 = 3.0;
    const index_type i_mid = n/2;
    
    if(out && print_all_tests())
      *out << std::endl << val1 << " -> *" << z_name << ".sub_view(1,"<<i_mid<<")\n";
    *z[k]->sub_view(1,i_mid) = val1;
    if(out && print_all_tests())
      *out << val2 << " -> *" << z_name << ".sub_view("<<i_mid+1<<","<<n<<")\n";
    *z[k]->sub_view(i_mid+1,n) = val2;
    if(out && print_vectors())
      *out << std::endl << z_name << " =\n" << *z[k] << std::endl;

    value_type
      norm_1        = z[k]->norm_1(),
      expect_norm_1 = (::fabs(val1)*(i_mid) + ::fabs(val2)*(n-i_mid));
    err = (norm_1 - expect_norm_1)/n;
    if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
      *out << "check: (" << z_name << "->norm_1() - |"<<val1<<"|*("<<i_mid<<")+"
         << "|"<<val2<<"|*("<<n<<"-"<<i_mid<<"))/"<<n
         <<" = ("<<norm_1<<" - "<<expect_norm_1<<")/"<<n<<" = "<<err<<std::endl;
    check_test(err,out,&success);

    value_type
      norm_2        = z[k]->norm_2(),
      expect_norm_2 = ::sqrt(val1*val1*(i_mid) + val2*val2*(n-i_mid));
    err = (norm_2 - expect_norm_2)/n;
    if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
      *out << "check: (" << z_name << "->norm_2() - ("<<val1<<")^2*("<<i_mid<<")+"
         << "("<<val2<<")^2*("<<n<<"-"<<i_mid<<"))/"<<n
         <<" = ("<<norm_2<<" - "<<expect_norm_2<<")/"<<n<<" = "<<err<<std::endl;
    check_test(err,out,&success);

    value_type
      norm_inf        = z[k]->norm_inf(),
      expect_norm_inf = my_max(::fabs(val1),::fabs(val2));
    err = (norm_inf - expect_norm_inf)/n;
    if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
      *out << "check: (" << z_name << "->norm_inf() - max(|"<<val1<<"|,"
         << "|"<<val2<<"|)/"<<n<<" = ("<<norm_inf<<" - "<<expect_norm_inf<<")/"<<n
         <<" = "<<err<<std::endl;
    check_test(err,out,&success);

  }
  else {
    if(out && print_all_tests())
      *out << "space.dim() <= 1, can't test the norms...\n";
  }
  
  //////////////////////////////////////////////////////////////////
  if(out && print_all_tests())
    *out << "\n*** Testing clone() method\n";
  {
    if(out && print_all_tests())
      *out << "\n(*vec = space.create_member()) = v[0]\n";
    VectorSpace::vec_mut_ptr_t
      vec = space.create_member();
    *vec = *v[0];
    if(out && print_all_tests())
      *out << "vec_clone = vec->clone()\n";
    VectorSpace::vec_mut_ptr_t
      vec_clone = vec->clone();
    if(out && print_all_tests())
      *out << "vec = NULL\n";
    vec = Teuchos::null;
    const value_type
      sum_vec       = sum(*v[0]),
      sum_vec_clone = sum(*vec_clone);
    err = (sum_vec - sum_vec_clone)/n;
    if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) )
      *out << "check: (sum(v[0]) - sum(vec_clone))/n = ("
         << sum_vec << " - " << sum_vec_clone << ")/" << n
         << " = " << err << std::endl;
    check_test(err,out,&success);
  }

  } // end try
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

void VectorSpaceTester::check_test(value_type err, std::ostream* out, bool* success) const
{
  if( ::fabs(err) >= error_tol() ) *success = false;
  if(out && (print_all_tests() || ::fabs(err) >= warning_tol()) ) {
    if( ::fabs(err) >= error_tol() )
      *out << "Error!  |" << err << "| = " << ::fabs(err) << " >= error_tol = "
         << error_tol() << std::endl;
    else if( ::fabs(err) >= warning_tol() )
      *out << "Warning!  |" << err << "| = " << ::fabs(err) << " >= warning_tol = "
             << warning_tol() << std::endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    !*success && this->throw_exception(), std::logic_error
    ,"VectorSpaceTester::check_test(...): Error!  |" << err << "| = " << ::fabs(err) << " >= error_tol = "
    << error_tol() << std::endl );
}

} // end namespace AbstractLinAlgPack
