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
//
// Definitions to template functions

#ifndef ITER_QUANITY_ACCESS_CONTINUOUS_DEF_H
#define ITER_QUANITY_ACCESS_CONTINUOUS_DEF_H

#include <typeinfo>
#include <algorithm>
#include <iterator>

#include "IterationPack_IterQuantityAccessContiguousDecl.hpp"
#include "Teuchos_Assert.hpp"

namespace IterationPack {

// Constructors/initializers

template<class T_info>
IterQuantityAccessContiguous<T_info>::IterQuantityAccessContiguous(
  int                              num_quantities
  ,const std::string&              name
  ,const abstract_factory_ptr_t&   abstract_factory
  )
  :num_quantities_(0)
  ,name_(name)
  ,abstract_factory_(abstract_factory)
  ,max_offset_(0)
{
  resize( num_quantities );
}

template<class T_info>
IterQuantityAccessContiguous<T_info>::~IterQuantityAccessContiguous() {
  release_mem();
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::set_factory(
  const abstract_factory_ptr_t& abstract_factory
  )
{
  release_mem();
  abstract_factory_ = abstract_factory;
  max_offset_ = std::numeric_limits<int>::min() + num_quantities_;  // uninitialized
  // ToDo: Don't wipe out storage, just reallocate new objects
  // as iteration quantities are updated.  This will take a little bit of
  // work and more overhead but it should be worth it in some types of
  // applications.
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::resize( int num_quantities ) {
  TEUCHOS_TEST_FOR_EXCEPTION(
    num_quantities < 1, std::length_error
    ,"IterQuantityAccessContiguous::resize(num_quantities): Error, "
    "name = "<<name_<<", num_quantities = "<<num_quantities<<" must be greater than zero" );
  if( num_quantities_ != num_quantities )
    release_mem();
  num_quantities_ = num_quantities;
  max_offset_ = std::numeric_limits<int>::min() + num_quantities_; // uninitialized
}

// Overridden from IterQuantity

template<class T_info>
IterQuantity* IterQuantityAccessContiguous<T_info>::clone() const {
  return 0;
  // ToDo: replace above with the following when the copy
  // constructor is implemented.
  // return new IterQuantityAccessContiguous<T_info>(*this);
}

template<class T_info>
const char* IterQuantityAccessContiguous<T_info>::name() const {
  return name_.c_str();
}

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::has_storage_k(int offset) const {
  return is_initialized()
    ? offset >= max_offset_ - num_quantities_ + 1
    : true;
}

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::updated_k(int offset) const {
  if( !is_initialized() )
    return false;
  if( ( offset > max_offset_ ) || ( offset < max_offset_ - num_quantities_ + 1 ) )
    return false;
  return updated_[max_offset_ - offset];
}

template<class T_info>
int IterQuantityAccessContiguous<T_info>::last_updated() const {
  if( !is_initialized() )
    return base_t::NONE_UPDATED;
  // Find the last still set as updated.
  for(	int offset = max_offset_;
      offset >= max_offset_ - num_quantities_ + 1;
      --offset										)
  {
    if( updated_[max_offset_ - offset] )
      return offset;
  }
  return base_t::NONE_UPDATED;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::set_not_updated_k(int offset) {
  this->assert_updated_k(offset);
  updated_[max_offset_ - offset] = false;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::set_all_not_updated() {
  if(!is_initialized()) return;
  std::fill( updated_.begin(), updated_.end(), false );
}

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::will_loose_mem(int offset, int set_offset) const {
  this->assert_updated_k(offset);
  return set_offset - max_offset_ > num_quantities_ - (max_offset_ - offset) - 1;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::next_iteration()
{
  if( !is_initialized() ) return;
  --max_offset_;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::print_concrete_type( std::ostream& out ) const
{
  const int last_updated = this->last_updated();
  if(last_updated != base_t::NONE_UPDATED)
    out << typeName(get_k(last_updated));
  else if( abstract_factory_.get() == NULL )
    out << "NULL";
  else
    out << typeName(*abstract_factory_->create());
}

// Overridden from IterQuantityAccess

template<class T_info>
T_info& IterQuantityAccessContiguous<T_info>::get_k(int offset) {
  this->assert_updated_k(offset);
  return *quantities_[max_offset_ - offset];
}
template<class T_info>
const T_info& IterQuantityAccessContiguous<T_info>::get_k(int offset) const {
  this->assert_updated_k(offset);
  return *quantities_[max_offset_ - offset];
}

template<class T_info>
T_info& IterQuantityAccessContiguous<T_info>::set_k(int offset) {

  lazy_initialization();

  this->assert_has_storage_k(offset);	// assert that we are not trying to iterate backwards

  if(offset > max_offset_ + num_quantities_ - 1) {
    // There will be no back memory so you don't need to adjust the pointers
    max_offset_ = offset;
    std::fill(updated_.begin(), updated_.end(), false);
  }
  else {
    // Pointers may have to be rearranged
    if(offset > max_offset_) {
      // We need to rearrange quantities_ and updated_.
      int shifted = offset - max_offset_;

      // ///////////////////////////////////
      // Set the updated flags

      // Example: We are shifting from:
      //		[1, 0, -1, -2] to [ 3, 2, 1, 0 ]

      // Shift the flags for the memory that may be saved
      updated_t::iterator
        itr_updated_from        = updated_.begin(),
        itr_updated_from_end    = itr_updated_from + num_quantities_ - shifted,
        itr_updated_to          = itr_updated_from + shifted;

      std::copy(itr_updated_from, itr_updated_from_end, itr_updated_to);

      // make updated[] for the new quantities false
      std::fill_n( updated_.begin(), shifted, false );

      // /////////////////////////////////////
      // rotate the quantitiy pointer vector

      // example: [1, 0, -1, -2] => [3, 2, 1, 0] so reverse rotate 0 to the back.

      // Rotate the (num_quantities_ - shifted) pointer to the back and the
      // remainder back arround again.

#if defined(_INTEL_CXX)
      typedef std::reverse_iterator<T_info**, T_info*, T_info*&
        , T_info**, ptrdiff_t>									rev_t;
#else
      typedef std::reverse_iterator<T_info**>						rev_t;
#endif

      std::rotate(
          rev_t(&quantities_[0] + num_quantities_)
        , rev_t(&quantities_[0] + num_quantities_ - shifted)
        , rev_t(&quantities_[0])                              );

      max_offset_ = offset;

    }
    // else, no pointers need to be rearranged since we are not advancing the range
  }

  updated_[max_offset_ - offset] = true;
  return *quantities_[max_offset_ - offset];
}

template<class T_info>
T_info& IterQuantityAccessContiguous<T_info>::set_k(int set_offset, int get_offset) {
  T_info& iq = this->get_k(get_offset);   // May throw exception
  return this->set_k(set_offset) = iq;    // ""
}

// Private member functions

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::is_initialized() const {
  return store_.size() > 0;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::lazy_initialization() {
  if( !is_initialized() ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      abstract_factory_.get() == NULL, std::logic_error
      ,"IterQuantityAccessContiguous::lazy_initialization(): Error, "
      "iq_name = "<<name_<<" the abstract factory can not be NULL" );
    // Allocate storage
    updated_.resize(num_quantities_,false);
    store_.resize(num_quantities_);
    quantities_.resize(num_quantities_,NULL);
    // Set initial points to locations
    typename updated_t::iterator       itr_updated         = updated_.begin();
    typename store_t::iterator         itr_store           = store_.begin();
    typename quantities_t::iterator    itr_quantities      = quantities_.begin();
    for( ; itr_store != store_.end(); ++itr_updated, ++itr_store, ++itr_quantities )
    {
      *itr_updated     = false;
      *itr_store       = abstract_factory_->create();
      *itr_quantities  = itr_store->get();
    }
    max_offset_ = std::numeric_limits<int>::min() + num_quantities_ + 1;
  }
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::release_mem() {
  updated_.resize(0);
  store_.resize(0);
  quantities_.resize(0);
}

}	// end namespace IterationPack

#endif	// ITER_QUANITY_ACCESS_CONTINUOUS_DEF_H
