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

#include "AbstractLinAlgPack_VectorSpaceSubSpace.hpp"
#include "AbstractLinAlgPack_VectorMutableSubView.hpp"
#include "Teuchos_Assert.hpp"

namespace AbstractLinAlgPack {

VectorSpaceSubSpace::VectorSpaceSubSpace( const space_ptr_t& full_space, const Range1D& rng )
{
  this->initialize(full_space,rng);
}

void VectorSpaceSubSpace::initialize( const space_ptr_t& full_space, const Range1D& rng )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    full_space.get() == NULL, std::invalid_argument
    ,"VectorSpaceSubSpace::initialize(...): Error!" );
#endif
  const index_type n = full_space->dim();
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !rng.full_range() && rng.ubound() > n, std::out_of_range
    ,"VectorSpaceSubSpace::initialize(...): Error, "
    "rng = [" << rng.lbound() << "," << rng.ubound() << "] is not in the range "
    "[1,vec->dim()] = [1," << n << "]" );
#endif
  full_space_ = full_space;
  rng_ = rng.full_range() ? Range1D(1,n) : rng;
}

void VectorSpaceSubSpace::set_uninitialized()
{
  full_space_ = Teuchos::null;
  rng_        = Range1D::Invalid;
}

#ifdef TEUCHOS_DEBUG
void VectorSpaceSubSpace::validate_range(const Range1D& rng) const
{
  const index_type n = this->dim();
  TEUCHOS_TEST_FOR_EXCEPTION(
    full_space_.get() == NULL, std::logic_error
    ,"VectorSpaceSubSpace::validate_range(rng): Error, Uninitialized" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    full_space_.get() && !rng.full_range() && rng.ubound() > n, std::logic_error
    ,"VectorSpaceSubSpace::validate_range(rng): Error, "
    "rng = [" << rng.lbound() << "," << rng.ubound() << "] is not in the range "
    "[1,this->dim] = [1," << n << "]" );
}
#endif

// Overridden from VectorSpace

bool VectorSpaceSubSpace::is_compatible(const VectorSpace& another_space) const
{
  if( this->dim() == another_space.dim() && this->is_in_core() && another_space.is_in_core() )
    return true;
  const VectorSpaceSubSpace
    *a_space = dynamic_cast<const VectorSpaceSubSpace*>(&another_space);
  if(!a_space)
    return false;
  return
    ( this->full_space_.get() == NULL && a_space->full_space_.get() == NULL )
    ||
    ( this->rng_ == a_space->rng_ && this->full_space_->is_compatible(*a_space->full_space_) );
}

bool VectorSpaceSubSpace::is_in_core() const
{
  return full_space_->is_in_core();
}

index_type VectorSpaceSubSpace::dim() const
{
  return full_space_.get() ? rng_.size() : 0;
}

VectorSpace::vec_mut_ptr_t VectorSpaceSubSpace::create_member() const
{
  namespace rcp = MemMngPack;
  if( full_space_.get() )
    return Teuchos::rcp(
      new VectorMutableSubView(
        full_space_->create_member(), rng_ 
        ) );
  return Teuchos::null;
}

VectorSpace::space_ptr_t VectorSpaceSubSpace::clone() const
{
  namespace rcp = MemMngPack;
  if( full_space_.get() )
    return Teuchos::rcp(new VectorSpaceSubSpace( full_space_->clone(), rng_ ));
  return Teuchos::rcp(new VectorSpaceSubSpace());
}

VectorSpace::space_ptr_t VectorSpaceSubSpace::sub_space(const Range1D& rng_in) const
{
  namespace rcp = MemMngPack;
  if( full_space_.get() == NULL && rng_in == Range1D::Invalid )
    return Teuchos::rcp(this,false);
  validate_range(rng_in);
  const index_type dim         = this->dim();
  const Range1D    rng         = rng_in.full_range() ? Range1D(1,dim) : rng_in;
  if( rng.lbound() == 1 && rng.ubound() == dim )
    return space_ptr_t( this, false );
  const index_type this_offset = rng_.lbound() - 1;
  return Teuchos::rcp(
    new VectorSpaceSubSpace(
      full_space_
      ,Range1D( 
        this_offset  + rng.lbound()
        ,this_offset + rng.ubound() )
      ) );
}

} // end namespace AbstractLinAlgPack
