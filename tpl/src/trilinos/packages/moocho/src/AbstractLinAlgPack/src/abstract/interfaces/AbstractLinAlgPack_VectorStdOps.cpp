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

#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_AssertOp.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorView.hpp"
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_ROp_max_abs_ele.h"
#include "RTOp_ROp_sum.h"
#include "RTOp_TOp_add_scalar.h"
#include "RTOp_TOp_axpy.h"
#include "RTOp_TOp_ele_wise_divide.h"
#include "RTOp_TOp_ele_wise_prod.h"
//#include "RTOp_TOp_random_vector.h"
#include "RTOpPack_TOpRandomize.hpp"
#include "RTOp_TOp_scale_vector.h"
#include "RTOp_TOp_sign.h"
#include "RTOpPack_RTOpC.hpp"
#include "Teuchos_Assert.hpp"

namespace {

// sum
static RTOpPack::RTOpC                               sum_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  sum_targ;
// dot prod
static RTOpPack::RTOpC                               dot_prod_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  dot_prod_targ;
// number of bounded elements
static RTOpPack::RTOpC                               num_bounded_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  num_bounded_targ;
// add scalar to vector
static RTOpPack::RTOpC                               add_scalar_op;
// scale vector
static RTOpPack::RTOpC                               scale_vector_op;
// axpy
static RTOpPack::RTOpC                               axpy_op;
// random vector
//static RTOpPack::RTOpC                               random_vector_op;
static RTOpPack::TOpRandomize<AbstractLinAlgPack::value_type>  random_vector_op;
// element-wise division
static RTOpPack::RTOpC                               ele_wise_divide_op;
// element-wise product
static RTOpPack::RTOpC                               ele_wise_prod_op;

// Simple class for an object that will initialize the RTOp_Server.
class init_rtop_server_t {
public:
  init_rtop_server_t() {
    // Operator and target obj for sum
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_sum_construct(&sum_op.op()));
    sum_targ = sum_op.reduct_obj_create();
    // Operator and target obj for dot_prod
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_dot_prod_construct(&dot_prod_op.op()));
    dot_prod_targ = dot_prod_op.reduct_obj_create();
    // Operator add_scalar
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_add_scalar_construct(0.0,&add_scalar_op.op()));
    // Operator scale_vector
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_scale_vector_construct(0.0,&scale_vector_op.op()));
    // Operator axpy
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_axpy_construct(0.0,&axpy_op.op()));
    // Operator random_vector
    //TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_random_vector_construct(0.0,0.0,&random_vector_op.op()));
    // Operator ele_wise_divide
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_ele_wise_divide_construct(0.0,&ele_wise_divide_op.op()));
    // Operator ele_wise_prod
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_ele_wise_prod_construct(0.0,&ele_wise_prod_op.op()));
  }
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

AbstractLinAlgPack::value_type
AbstractLinAlgPack::sum( const Vector& v_rhs )
{
  sum_op.reduct_obj_reinit(sum_targ.ptr());
  const Vector* vecs[1] = { &v_rhs };
  apply_op(sum_op,1,vecs,0,NULL,&*sum_targ);
  return RTOp_ROp_sum_val(sum_op(*sum_targ));
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::dot( const Vector& v_rhs1, const Vector& v_rhs2 )
{
  dot_prod_op.reduct_obj_reinit(dot_prod_targ.ptr());
  const Vector* vecs[2] = { &v_rhs1, &v_rhs2 };
  apply_op(dot_prod_op,2,vecs,0,NULL,&*dot_prod_targ);
  return RTOp_ROp_dot_prod_val(dot_prod_op(*dot_prod_targ));
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::dot( const Vector& v_rhs1, const SpVectorSlice& sv_rhs2 )
{
  VopV_assert_compatibility(v_rhs1,sv_rhs2 );
  if( sv_rhs2.nz() ) {
    VectorSpace::vec_mut_ptr_t
      v_rhs2 = v_rhs1.space().create_member();
    v_rhs2->set_sub_vector(sub_vec_view(sv_rhs2));
    return dot(v_rhs1,*v_rhs2);
  }
  return 0.0;
}

void AbstractLinAlgPack::max_abs_ele(
  const Vector& v, value_type* max_v_j, index_type* max_j
  )
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  max_v_j && max_j  ) );
  RTOpPack::RTOpC op;
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_max_abs_ele_construct(&op.op()));
  Teuchos::RCP<RTOpPack::ReductTarget> reduct_obj = op.reduct_obj_create();
  const Vector* vecs[1] = { &v };
  apply_op(op,1,vecs,0,NULL,&*reduct_obj);
  RTOp_value_index_type val = RTOp_ROp_max_abs_ele_val(op(*reduct_obj));
  *max_v_j = val.value;
  *max_j   = val.index;
}

void AbstractLinAlgPack::Vp_S( VectorMutable* v_lhs, const value_type& alpha )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vp_S(...), Error!");
#endif
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_add_scalar_set_alpha(alpha,&add_scalar_op.op()));
  VectorMutable* targ_vecs[1] = { v_lhs };
  apply_op(add_scalar_op,0,NULL,1,targ_vecs,NULL);
}

void AbstractLinAlgPack::Vt_S( VectorMutable* v_lhs, const value_type& alpha )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vt_S(...), Error!");
#endif
  if( alpha == 0.0 ) {
    *v_lhs = 0.0;
  }
  else if( alpha != 1.0 ) {
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_scale_vector_set_alpha( alpha, &scale_vector_op.op() ));
    VectorMutable* targ_vecs[1] = { v_lhs };
    apply_op(scale_vector_op,0,NULL,1,targ_vecs,NULL);
  }
}

void AbstractLinAlgPack::Vp_StV(
  VectorMutable* v_lhs, const value_type& alpha, const Vector& v_rhs)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"Vp_StV(...), Error!");
#endif
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_axpy_set_alpha( alpha, &axpy_op.op() ));
  const Vector*  vecs[1]      = { &v_rhs };
  VectorMutable* targ_vecs[1] = { v_lhs  };
  apply_op(axpy_op,1,vecs,1,targ_vecs,NULL);
}

void AbstractLinAlgPack::Vp_StV(
  VectorMutable* y, const value_type& a, const SpVectorSlice& sx )
{
  Vp_V_assert_compatibility(y,sx);
  if( sx.nz() ) {
    VectorSpace::vec_mut_ptr_t
        x = y->space().create_member();
    x->set_sub_vector(sub_vec_view(sx));
    Vp_StV( y, a, *x );
  }
}

void AbstractLinAlgPack::ele_wise_prod(
  const value_type& alpha, const Vector& v_rhs1, const Vector& v_rhs2
  , VectorMutable* v_lhs )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_prod(...), Error");
#endif
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_ele_wise_prod_set_alpha(alpha,&ele_wise_prod_op.op()));
  const Vector*   vecs[2]      = { &v_rhs1, &v_rhs2 };
  VectorMutable*  targ_vecs[1] = { v_lhs };
  apply_op(ele_wise_prod_op,2,vecs,1,targ_vecs,NULL);
}

void AbstractLinAlgPack::ele_wise_divide(
  const value_type& alpha, const Vector& v_rhs1, const Vector& v_rhs2
  , VectorMutable* v_lhs )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(v_lhs==NULL,std::logic_error,"ele_wise_divide(...), Error");
#endif
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_ele_wise_divide_set_alpha(alpha,&ele_wise_divide_op.op()));
  const int num_vecs = 2;
  const Vector*   vecs[2]      = { &v_rhs1, &v_rhs2 };
  VectorMutable*  targ_vecs[1] = { v_lhs };
  apply_op(ele_wise_divide_op,2,vecs,1,targ_vecs,NULL);
}

void AbstractLinAlgPack::seed_random_vector_generator( unsigned int s )
{
  random_vector_op.set_seed(s);
}

void AbstractLinAlgPack::random_vector( value_type l, value_type u, VectorMutable* v )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(v==NULL,std::logic_error,"Vt_S(...), Error");
#endif
  //TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_random_vector_set_bounds( l, u, &random_vector_op.op() ));
  random_vector_op.set_bounds(l,u);
  VectorMutable* targ_vecs[1] = { v };
  apply_op(random_vector_op,0,NULL,1,targ_vecs,NULL);
  
}

void AbstractLinAlgPack::sign(
  const Vector      &v
  ,VectorMutable    *z
  )
{
  RTOpPack::RTOpC op;
  TEUCHOS_TEST_FOR_EXCEPT( !( 0==RTOp_TOp_sign_construct(&op.op()) ) );
  const Vector*   vecs[1]      = { &v };
  VectorMutable*  targ_vecs[1] = { z  };
  apply_op(op,1,vecs,1,targ_vecs,NULL);
}
