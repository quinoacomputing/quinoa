// $Id$
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_DYNAMIC_STORAGE_HPP
#define STOKHOS_DYNAMIC_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

namespace Stokhos {

  //! Dynamically allocated storage class
  template <typename ordinal_t, typename value_t, typename node_t>
  class DynamicStorage {
  public:

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef node_t node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::DynArrayTraits<value_type,node_type> ds;

    //! Turn DynamicStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef DynamicStorage<ord_t,val_t,node_type> type;
    };

    //! Constructor
    DynamicStorage(const ordinal_type& sz,
		  const value_type& x = value_type(0.0));

    //! Copy constructor
    DynamicStorage(const DynamicStorage& s);

    //! Destructor
    ~DynamicStorage();

    //! Assignment operator
    DynamicStorage& operator=(const DynamicStorage& s);

    //! Initialize values to a constant value
    void init(const_reference v);

    //! Initialize values to an array of values
    void init(const_pointer v, const ordinal_type& sz_ = 0);

    //! Load values to an array of values
    void load(pointer v);

    //! Resize to new size (values are preserved)
    void resize(const ordinal_type& sz);

    //! Return size
    static ordinal_type size();

    //! Coefficient access (avoid if possible)
    const_reference operator[] (const ordinal_type& i) const;

    //! Coefficient access (avoid if possible)
    reference operator[] (const ordinal_type& i);

    //! Get coefficients
    const_pointer coeff() const;

    //! Get coefficients
    pointer coeff();

  };

}

// Host specialization
#include "KokkosArray_Host.hpp"
#include "KokkosArray_Host_macros.hpp"
#include "Stokhos_DynamicStorage_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

// Cuda specialization
#include "KokkosArray_Cuda.hpp"
#include "KokkosArray_Cuda_macros.hpp"
#include "Stokhos_DynamicStorage_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

#endif // STOKHOS_DYNAMIC_STORAGE_HPP
