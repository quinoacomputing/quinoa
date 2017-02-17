// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_DYNAMICSTORAGE_HPP
#define SACADO_FAD_DYNAMICSTORAGE_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_DynamicArrayTraits.hpp"

namespace Sacado {

  namespace Fad {

    //! Derivative array storage class using dynamic memory allocation
    template <typename T, typename S = T> 
    class DynamicStorage {

    public:

      //! Default constructor
      DynamicStorage(const T & x) : 
	val_(x), sz_(0), len_(0), dx_(NULL) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      DynamicStorage(const int sz, const T & x) : 
	val_(x), sz_(sz), len_(sz) {
	dx_ = ds_array<S>::get_and_fill(sz_);
      }

      //! Copy constructor
      DynamicStorage(const DynamicStorage& x) : 
	val_(x.val_), sz_(x.sz_), len_(x.sz_) {
	dx_ = ds_array<S>::get_and_fill(x.dx_, sz_);
      }
      
      //! Destructor
      ~DynamicStorage() {
	if (len_ != 0)
	  ds_array<S>::destroy_and_release(dx_, len_);
      }

      //! Assignment
      DynamicStorage& operator=(const DynamicStorage& x) { 
	val_ = x.val_;
	if (sz_ != x.sz_) {
	  sz_ = x.sz_;
	  if (x.sz_ > len_) {
	    if (len_ != 0)
	      ds_array<S>::destroy_and_release(dx_, len_);
	    len_ = x.sz_;
	    dx_ = ds_array<S>::get_and_fill(x.dx_, sz_);
	  }
	  else 
	    ds_array<S>::copy(x.dx_, dx_, sz_);
	}
	else 
	  ds_array<S>::copy(x.dx_, dx_, sz_);

	return *this; 
      } 

      //! Returns number of derivative components
      int size() const { return sz_;}

      //! Returns array length
      int length() const { return len_; }

      //! Resize the derivative array to sz
      /*!
       * Note:  This does not necessarily preserve derivative components.
       */
      void resize(int sz) { 
	if (sz > len_) {
	  if (len_ != 0)
	    ds_array<S>::destroy_and_release(dx_, len_);
	  dx_ = ds_array<S>::get_and_fill(sz);
	  len_ = sz;
	}
	sz_ = sz;
      }

      //! Expand derivative array to size sz
      /*!
       * This method preserves any existing derivative components and
       * sets any that are added to zero.
       */
      void expand(int sz) {
        if (sz > len_) {
          S* dx_new = ds_array<S>::get_and_fill(sz);
          ds_array<S>::copy(dx_, dx_new, sz_);
          if (len_ > 0)
            ds_array<S>::destroy_and_release(dx_, len_);
          dx_ = dx_new;
          len_ = sz;
        }
        else if (sz > sz_) 
          ds_array<S>::zero(dx_+sz_, sz-sz_);
        sz_ = sz;
      }

      //! Zero out derivative array
      void zero() { 
	ds_array<S>::zero(dx_, sz_);
      }

      //! Returns value
      const T& val() const { return val_; }

      //! Returns value
      T& val() { return val_; }

      //! Returns derivative array
      const S* dx() const { return dx_;}

      //! Returns derivative component \c i with bounds checking
      S dx(int i) const { return sz_ ? dx_[i] : T(0.); }
    
      //! Returns derivative component \c i without bounds checking
      S& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      const S& fastAccessDx(int i) const { return dx_[i];}

    private:

      //! Value
      T val_;

      //! Derivative array size
      int sz_;

      //! Derivative array length
      int len_;

      //! Derivative array
      S* dx_;

    }; // class DynamicStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_DYNAMICSTORAGE_HPP
