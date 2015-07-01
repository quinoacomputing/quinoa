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

#include <typeinfo>
#include <stdexcept>

#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_apply_op_helper.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"

#ifdef TEUCHOS_DEBUG
#define CLASS_MEMBER_PTRS \
const VectorMutableDense   *_this = this; \
const DVectorSlice                *_v; \
const release_resource_ptr_t     *_v_release; \
const VectorSpaceSerial          *_space;
#else
#define CLASS_MEMBER_PTRS
#endif

namespace AbstractLinAlgPack {

VectorMutableDense::VectorMutableDense(
  const size_type                   dim
  )
  :space_(dim)
{
  CLASS_MEMBER_PTRS
  this->initialize(dim);
}

VectorMutableDense::VectorMutableDense(
  DVectorSlice                        v
  ,const release_resource_ptr_t&     v_release
  )
  :space_(v.dim())
{
  CLASS_MEMBER_PTRS
  this->initialize(v,v_release);
}

void VectorMutableDense::initialize(
  const size_type                   dim
  )
{
  CLASS_MEMBER_PTRS
  namespace rcp = MemMngPack;
  namespace rmp = MemMngPack;
  typedef Teuchos::RCP<DVector> vec_ptr_t;
  vec_ptr_t vec_ptr = Teuchos::rcp(new DVector(dim));
  this->initialize(
    (*vec_ptr)()
    ,Teuchos::rcp(
      new rmp::ReleaseResource_ref_count_ptr<DVector>(
        vec_ptr
        )
      )
    );
}

void VectorMutableDense::initialize(
  DVectorSlice                       v
  ,const release_resource_ptr_t&     v_release
  )
{
  CLASS_MEMBER_PTRS
  v_.bind(v);
  v_release_ = v_release;
  space_.initialize(v.dim());
  this->has_changed();
}

// Overridden from Vector

const VectorSpace& VectorMutableDense::space() const
{
  CLASS_MEMBER_PTRS
  return space_;
}

void VectorMutableDense::apply_op(
  const RTOpPack::RTOp& op
  ,const size_t num_vecs,      const Vector*  vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
  ) const
{
  CLASS_MEMBER_PTRS
#ifdef TEUCHOS_DEBUG
  AbstractLinAlgPack::apply_op_validate_input(
    "VectorMutableDense::apply_op(...)"
    ,op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj,first_ele_in,sub_dim_in,global_offset_in
    );
#endif
  this->apply_op_serial(
    op,num_vecs,vecs,num_targ_vecs,targ_vecs,reduct_obj
    ,first_ele_in,sub_dim_in,global_offset_in
    );
}

index_type VectorMutableDense::dim() const
{
  return v_.dim();
}

value_type VectorMutableDense::get_ele(index_type i) const
{
  return v_(i);
}

void VectorMutableDense::get_sub_vector(
  const Range1D& rng_in, RTOpPack::SubVector* sub_vec
  ) const
{
  CLASS_MEMBER_PTRS
  const size_type  this_dim = v_.dim();
  const Range1D    rng = RangePack::full_range(rng_in,1,this_dim);
  TEUCHOS_TEST_FOR_EXCEPTION(
    rng.ubound() > this_dim, std::out_of_range
    ,"VectorMutableDense::get_sub_vector(...) : Error, "
    "rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
    "is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
  // Just return the dense view regardless of spare_or_dense argument
  sub_vec->initialize(
    rng.lbound()-1                             // global_offset
    ,rng.size()                                // sub_dim
    ,v_.raw_ptr()+v_.stride()*(rng.lbound()-1) // values
    ,v_.stride()
    );
}

void VectorMutableDense::free_sub_vector( RTOpPack::SubVector* sub_vec ) const
{
  sub_vec->set_uninitialized(); // No memory to deallocate!
}

// Overridden from VectorMutable

VectorMutable&
VectorMutableDense::operator=(value_type alpha)
{
  CLASS_MEMBER_PTRS
  v_ = alpha;
  this->has_changed();
  return *this;
}

VectorMutable&
VectorMutableDense::operator=(const Vector& v)
{
  CLASS_MEMBER_PTRS
  if( const VectorMutableDense *vp = dynamic_cast<const VectorMutableDense*>(&v) )
    v_ = vp->v_;
  else
    return VectorMutable::operator=(v); // Try the default implementation?
  this->has_changed();
  return *this;
}

VectorMutable&
VectorMutableDense::operator=(const VectorMutable& v)
{
  CLASS_MEMBER_PTRS
  if( const VectorMutableDense *vp = dynamic_cast<const VectorMutableDense*>(&v) )
    v_ = vp->v_;
  else
    return VectorMutable::operator=(v); // Try the default implementation?
  this->has_changed();
  return *this;
}

void VectorMutableDense::set_ele( index_type i, value_type val )
{
  CLASS_MEMBER_PTRS
  v_(i) = val;
  this->has_changed();
}

VectorMutableDense::vec_mut_ptr_t
VectorMutableDense::sub_view( const Range1D& rng_in )
{
  CLASS_MEMBER_PTRS
  namespace rcp = MemMngPack;
  const size_type this_dim = this->dim();
  const Range1D rng = RangePack::full_range( rng_in, 1, this_dim );
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    rng.ubound() > this_dim, std::out_of_range
    ,"VectorMutableDense::sub_view(...) : Error, "
    "rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
    "is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
#endif
  if( rng == Range1D(1,this_dim) )
    return Teuchos::rcp( this, false );
  this->has_changed(); // This will result in a change in the vector
  return Teuchos::rcp( new VectorMutableDense( v_(rng), Teuchos::null ) ); 
}

void VectorMutableDense::get_sub_vector(
  const Range1D& rng_in, RTOpPack::MutableSubVector* sub_vec )
{
  CLASS_MEMBER_PTRS
  const size_type  this_dim = v_.dim();
  const Range1D    rng = RangePack::full_range(rng_in,1,this_dim);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    rng.ubound() > this_dim, std::out_of_range
    ,"VectorMutableDense::get_sub_vector(...) : Error, "
    "rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
    "is not in the range [1,this->dim()] = [1," << this_dim << "]!" );
#endif
  sub_vec->initialize(
    rng.lbound()-1                             // global_offset
    ,rng.size()                                // sub_dim
    ,v_.raw_ptr()+v_.stride()*(rng.lbound()-1) // values
    ,v_.stride()
    );
}

void VectorMutableDense::commit_sub_vector( RTOpPack::MutableSubVector* sub_vec )
{
  CLASS_MEMBER_PTRS
  sub_vec->set_uninitialized(); // No memory to deallocate!
  this->has_changed();          // Be aware of any final changes!
}

void VectorMutableDense::set_sub_vector( const RTOpPack::SparseSubVector& sub_vec )
{
  CLASS_MEMBER_PTRS
  VectorMutable::set_sub_vector(sub_vec); // ToDo: Provide specialized implementation?
}

void VectorMutableDense::Vp_StMtV(
  value_type                       alpha
  ,const GenPermMatrixSlice        &P
  ,BLAS_Cpp::Transp                P_trans
  ,const Vector              &x
  ,value_type                      beta
  )
{
  CLASS_MEMBER_PTRS
  VectorDenseEncap  x_de(x);
  AbstractLinAlgPack::Vp_StMtV( &v_, alpha, P, P_trans, x_de(), beta );
}

} // end namespace AbstractLinAlgPack
