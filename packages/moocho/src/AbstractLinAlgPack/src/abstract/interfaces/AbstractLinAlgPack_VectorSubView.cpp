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

#include <stdexcept>

#include "AbstractLinAlgPack_VectorMutableSubView.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

VectorSubView::VectorSubView( const vec_ptr_t& vec, const Range1D& rng )
{
  initialize(vec,rng);
}

void VectorSubView::initialize( const vec_ptr_t& vec, const Range1D& rng )
{
  namespace rcp = MemMngPack;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    vec.get() == NULL, std::invalid_argument
    ,"VectorSubView::initialize(...) : Error!" );
#endif
  typedef VectorSpace::space_ptr_t   space_ptr_t;
  space_.initialize( Teuchos::rcp(&vec->space(),false), rng );
  full_vec_ = vec;
  this->has_changed();
}

void VectorSubView::set_uninitialized()
{
  space_.set_uninitialized();
  full_vec_ = Teuchos::null;
  this->has_changed();
}

// Overridden from Vector

const VectorSpace& VectorSubView::space() const
{
  return space_;
}

index_type VectorSubView::dim() const
{
  return space_.dim();
}

void VectorSubView::apply_op(
  const RTOpPack::RTOp& op
  ,const size_t num_vecs, const Vector* vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget* reduct_obj
  ,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  // If these are in-core vectors then just let a default implementation
  // take care of this.
  if( this->space().is_in_core() ) {
    this->apply_op_serial(
      op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
      ,first_ele_in,sub_dim_in,global_offset_in
      );
    return;
  }
  // These are not in-core vectors so ...
  int k;
  const index_type this_dim = this->dim();
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    sub_dim_in < 0
    || !(1 <= first_ele_in && first_ele_in <= this_dim)
    || ( sub_dim_in > 0 && (sub_dim_in - (first_ele_in - 1) > this_dim) )
    , std::logic_error
    ,"VectorSubView::apply_op(...): Error, first_ele_in = "
    << first_ele_in << ", global_offset_in = " << global_offset_in
    << ", sub_dim_in = " << sub_dim_in << " and this->dim() = this_dim  = "
    << this_dim << " are not compatible." );
#endif
  const index_type this_offset = space_impl().rng().lbound() - 1;
  const index_type
    this_sub_dim = ( sub_dim_in 
             ? sub_dim_in
             : space_impl().rng().size() - (first_ele_in - 1)
                 );
  Workspace<const Vector*>    vecs_full(wss,num_vecs);
  for( k = 0; k < num_vecs; ++k )
    vecs_full[k] = dyn_cast<const VectorSubView>(*vecs[k]).full_vec().get();
  Workspace<VectorMutable*>   targ_vecs_full(wss,num_targ_vecs);
  for( k = 0; k < num_targ_vecs; ++k )
    targ_vecs_full[k] = dyn_cast<VectorMutableSubView>(*targ_vecs[k]).full_vec().get();
  AbstractLinAlgPack::apply_op(
    op
    ,num_vecs,      num_vecs      ? &vecs_full[0]      : NULL
    ,num_targ_vecs, num_targ_vecs ? &targ_vecs_full[0] : NULL
    ,reduct_obj
    ,this_offset + first_ele_in     // first_ele
    ,this_sub_dim                   // sub_dim
    ,global_offset_in               // global_dim
    );
}

value_type VectorSubView::get_ele(index_type i) const
{
  space_.validate_range(Range1D(i,i));
  return full_vec_->get_ele( space_.rng().lbound() + i - 1 );
}

Vector::vec_ptr_t
VectorSubView::sub_view( const Range1D& rng_in ) const
{
  namespace rcp = MemMngPack;
  const index_type this_dim = this->dim();
  const Range1D rng = RangePack::full_range(rng_in,1,this_dim);
  space_.validate_range(rng);
  const index_type this_offset = space_.rng().lbound() - 1;
  return Teuchos::rcp(
    new VectorSubView(
      full_vec_
      ,Range1D( 
        this_offset  + rng.lbound()
        ,this_offset + rng.ubound() )
      ) );
}

void VectorSubView::get_sub_vector( const Range1D& rng_in, RTOpPack::SubVector* sub_vec ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !sub_vec, std::logic_error, "VectorSubView::get_sub_vector(...): Error!" );
#endif
  const index_type this_dim = this->dim();
  const Range1D rng = RangePack::full_range(rng_in,1,this_dim);
  space_.validate_range(rng);
  const index_type this_offset = space_.rng().lbound() - 1;
  full_vec_->get_sub_vector( rng + this_offset, sub_vec );
  sub_vec->setGlobalOffset( sub_vec->globalOffset() - this_offset );
}

void VectorSubView::free_sub_vector( RTOpPack::SubVector* sub_vec ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !sub_vec, std::logic_error, "VectorSubView::free_sub_vector(...): Error!" );
#endif
  const index_type this_offset = space_.rng().lbound() - 1;
  sub_vec->setGlobalOffset( sub_vec->globalOffset() + this_offset );
  full_vec_->free_sub_vector( sub_vec );
}

} // end namespace AbstractLinAlgPack
