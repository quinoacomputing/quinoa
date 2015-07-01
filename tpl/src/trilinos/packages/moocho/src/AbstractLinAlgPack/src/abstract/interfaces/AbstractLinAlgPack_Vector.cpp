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


#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSubView.hpp"
#include "RTOp_ROp_dot_prod.h"
#include "RTOp_ROp_sum.h"
#include "RTOp_ROp_norms.h"
#include "RTOp_ROp_num_nonzeros.h"
#include "RTOp_ROp_get_sub_vector.h"
#include "RTOpPack_RTOpC.hpp"
#include "RTOpPack_print_sub_vector.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

#include <limits>
#include <ostream>

#include <assert.h>


// Uncomment to ignore cache of reduction data
//#define ALAP_VECTOR_IGNORE_CACHE_DATA


namespace {


// get element operator
static RTOpPack::RTOpC                               sum_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  sum_targ;
// number nonzros
static RTOpPack::RTOpC                               num_nonzeros_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  num_nonzeros_targ;
// Norm 1
static RTOpPack::RTOpC                               norm_1_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  norm_1_targ;
// Norm 2
static RTOpPack::RTOpC                               norm_2_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  norm_2_targ;
// Norm inf
static RTOpPack::RTOpC                               norm_inf_op;
static Teuchos::RCP<RTOpPack::ReductTarget>  norm_inf_targ;
// get sub-vector operator
static RTOpPack::RTOpC                               get_sub_vector_op;


// Simple class for an object that will initialize the RTOp_Server and operators.
class init_rtop_server_t {
public:
  init_rtop_server_t() {
    // Operator and target for getting a vector element
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_sum_construct(&sum_op.op()));
    sum_targ = sum_op.reduct_obj_create();
    // Operator and target for norm 1
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_num_nonzeros_construct(&num_nonzeros_op.op()));
    num_nonzeros_targ = num_nonzeros_op.reduct_obj_create();
    // Operator and target for norm 1
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_norm_1_construct(&norm_1_op.op()));
    norm_1_targ = norm_1_op.reduct_obj_create();
    // Operator and target for norm 1
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_norm_2_construct(&norm_2_op.op()));
    norm_2_targ = norm_2_op.reduct_obj_create();
    // Operator and target for norm 1
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_norm_inf_construct(&norm_inf_op.op()));
    norm_inf_targ = norm_inf_op.reduct_obj_create();
    // Get sub-vector operator
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_get_sub_vector_construct(1,1,&get_sub_vector_op.op()));
  }
}; 


// When the program starts, this object will be created
init_rtop_server_t  init_rtop_server;


} // end namespace


namespace AbstractLinAlgPack {


Vector::Vector()
  :num_nonzeros_(-1), norm_1_(-1.0), norm_2_(-1.0), norm_inf_(-1.0) // uninitalized
{}


index_type Vector::dim() const
{
  return this->space().dim();
}


index_type Vector::nz() const
{
#ifdef ALAP_VECTOR_IGNORE_CACHE_DATA
  {
#else
  if( num_nonzeros_ < 0 ) {
#endif
    num_nonzeros_op.reduct_obj_reinit(num_nonzeros_targ.ptr());
    const Vector *vecs[1] = { this };
    AbstractLinAlgPack::apply_op(num_nonzeros_op,1,vecs,0,NULL,&*num_nonzeros_targ);
    num_nonzeros_ = RTOp_ROp_num_nonzeros_val(num_nonzeros_op(*num_nonzeros_targ));
  }
  return num_nonzeros_;
}


std::ostream& Vector::output(
  std::ostream& out_arg, bool print_dim, bool newline,
  index_type global_offset
  ) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  RTOpPack::SubVector sub_vec;
  this->get_sub_vector( Range1D(), &sub_vec );
  RTOpPack::SubVector sub_vec_print( sub_vec.globalOffset() + global_offset, sub_vec.subDim(), sub_vec.values(), sub_vec.stride() );
  RTOpPack::output(*out,sub_vec_print,print_dim,newline);
  this->free_sub_vector( &sub_vec );
  return out_arg;
}


VectorMutable::vec_mut_ptr_t Vector::clone() const
{
  vec_mut_ptr_t
    vec = this->space().create_member();
  *vec = *this;
  return vec;
}


value_type Vector::get_ele(index_type i) const
{
  sum_op.reduct_obj_reinit(sum_targ.ptr());
  const Vector *vecs[1] = { this };
  AbstractLinAlgPack::apply_op(
    sum_op,1,vecs,0,NULL,&*sum_targ
    ,i,1,0 // first_ele, sub_dim, global_offset
    );
   return RTOp_ROp_sum_val(sum_op(*sum_targ));
}


value_type Vector::norm_1() const
{
#ifdef ALAP_VECTOR_IGNORE_CACHE_DATA
  {
#else
  if( norm_1_ < 0.0 ) {
#endif
    norm_1_op.reduct_obj_reinit(norm_1_targ.ptr());
    const Vector *vecs[1] = { this };
    AbstractLinAlgPack::apply_op(norm_1_op,1,vecs,0,NULL,&*norm_1_targ);
    norm_1_ = RTOp_ROp_norm_1_val(norm_1_op(*norm_1_targ));
  }
  return norm_1_;
}


value_type Vector::norm_2() const
{
#ifdef ALAP_VECTOR_IGNORE_CACHE_DATA
  {
#else
  if( norm_2_ < 0.0 ) {
#endif
    norm_2_op.reduct_obj_reinit(norm_2_targ.ptr());
    const Vector *vecs[1] = { this };
    AbstractLinAlgPack::apply_op(norm_2_op,1,vecs,0,NULL,&*norm_2_targ);
    norm_2_ = RTOp_ROp_norm_2_val(norm_2_op(*norm_2_targ));
  }
  return norm_2_;
}


value_type Vector::norm_inf() const
{
#ifdef ALAP_VECTOR_IGNORE_CACHE_DATA
  {
#else
  if( norm_inf_ < 0.0 ) {
#endif
    norm_inf_op.reduct_obj_reinit(norm_inf_targ.ptr());
    const Vector *vecs[1] = { this };
    AbstractLinAlgPack::apply_op(norm_inf_op,1,vecs,0,NULL,&*norm_inf_targ);
    norm_inf_ = RTOp_ROp_norm_inf_val(norm_inf_op(*norm_inf_targ));
  }
  return norm_inf_;
}


value_type Vector::inner_product(  const Vector& v ) const
{
  return this->space().inner_prod()->inner_prod(*this,v);
}

Vector::vec_ptr_t
Vector::sub_view( const Range1D& rng_in ) const
{
  namespace rcp = MemMngPack;
  const index_type dim = this->dim();
  const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    rng.ubound() > dim, std::out_of_range
    ,"Vector::sub_view(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
    "is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
  if( rng.lbound() == 1 && rng.ubound() == dim )
    return vec_ptr_t( this, false );
  return Teuchos::rcp(
    new VectorSubView(
      vec_ptr_t( this, false )
      ,rng ) );
}


void Vector::get_sub_vector( const Range1D& rng_in, RTOpPack::SubVector* sub_vec_inout ) const
{
  const Range1D rng = rng_in.full_range() ? Range1D(1,this->space().dim()) : rng_in;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->space().dim() < rng.ubound(), std::out_of_range
    ,"Vector::get_sub_vector(rng,...): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()
    <<"] is not in range = [1,"<<this->space().dim()<<"]" );
#endif
  // Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
  if( sub_vec_inout->values() ) {
    std::free( (void*)sub_vec_inout->values()  );
  }
  // Initialize the operator
  RTOpPack::RTOpC get_sub_vector_op;
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_ROp_get_sub_vector_construct(rng.lbound(),rng.ubound(),&get_sub_vector_op.op()));
  // Create the reduction object (another sub_vec)
  Teuchos::RCP<RTOpPack::ReductTarget> reduct_obj = get_sub_vector_op.reduct_obj_create(); // This is really of type RTOpPack::ConstSubVectorView<Scalar>!
  // Perform the reduction (get the sub-vector requested)
  const size_t  num_vecs = 1;
  const Vector* sub_vecs[num_vecs] = { this };
  AbstractLinAlgPack::apply_op(
    get_sub_vector_op,num_vecs,sub_vecs,0,NULL,&*reduct_obj
    ,rng.lbound(),rng.size(),rng.lbound()-1 // first_ele, sub_dim, global_offset
    );
  // Set the sub-vector.  Note reduct_obj will go out of scope so the sub_vec parameter will
  // own the memory allocated within get_sub_vector_op.create_reduct_obj_raw(...).  This is okay
  // since the client is required to call release_sub_vector(...) so release memory!
  RTOp_ReductTarget reduct_obj_raw = get_sub_vector_op(*reduct_obj);
  RTOp_SubVector sub_vec = RTOp_ROp_get_sub_vector_val(reduct_obj_raw);
  sub_vec_inout->initialize(sub_vec.global_offset,sub_vec.sub_dim,sub_vec.values,sub_vec.values_stride);\
  reduct_obj.release();  // Do not allow delete to be called!
  std::free(reduct_obj_raw); // Now *sub_vec owns the values[] and indices[] arrays!
}


void Vector::free_sub_vector( RTOpPack::SubVector* sub_vec ) const
{
  // Free sub_vec if needed (note this is dependent on the implemenation of this operator class!)
  if( sub_vec->values() )
    std::free( (void*)sub_vec->values() );
  sub_vec->set_uninitialized();
}


void Vector::has_changed() const
{
  num_nonzeros_= -1;  // uninitalized;
  norm_1_ = norm_2_ = norm_inf_ = -1.0;
}


// protected


void Vector::finalize_apply_op(
  const size_t num_targ_vecs, VectorMutable** targ_vecs
  ) const
{
  for( size_t k = 0; k < num_targ_vecs; ++k )
    targ_vecs[k]->has_changed();
}


} // end namespace AbstractLinAlgPack


// nonmember functions


void AbstractLinAlgPack::apply_op(
  const RTOpPack::RTOp       &op
  ,const size_t              num_vecs
  ,const Vector*             vecs[]
  ,const size_t              num_targ_vecs
  ,VectorMutable*            targ_vecs[]
  ,RTOpPack::ReductTarget    *reduct_obj
  ,const index_type          first_ele
  ,const index_type          sub_dim
  ,const index_type          global_offset
  )
{
  if(num_vecs)
    vecs[0]->apply_op(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
  else if (num_targ_vecs)
    targ_vecs[0]->apply_op(op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele,sub_dim,global_offset);
}
