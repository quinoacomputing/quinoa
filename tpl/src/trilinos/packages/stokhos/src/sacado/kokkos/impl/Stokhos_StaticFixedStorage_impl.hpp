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

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Stokhos_StaticFixedStorage_impl.hpp> without macros defined"

#else

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t, int Num>
  class StaticFixedStorage<ordinal_t, value_t, Num, KOKKOS_MACRO_DEVICE> {
  public:

    static const bool is_static = true;
    static const int static_size = Num;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef KOKKOS_MACRO_DEVICE node_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::StaticArrayTraits<value_type,node_type> ss;

    //! Turn StaticFixedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef StaticFixedStorage<ord_t,val_t,Num,node_type> type;
    };

    //! Constructor
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    StaticFixedStorage(const ordinal_type& sz,
		       const value_type& x = value_type(0.0)) { 
      ss::fill(coeff_, Num, x); 
    }

    //! Copy constructor
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    StaticFixedStorage(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Destructor
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    ~StaticFixedStorage() {}

    //! Assignment operator
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    StaticFixedStorage& operator=(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
      return *this;
    }

    //! Initialize values to a constant value
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_reference v) { 
      ss::fill(coeff_, Num, v); 
    }

    //! Initialize values to an array of values
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
      	ss::copy(v, coeff_, Num);
      else
      	ss::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    void load(pointer v) { 
      ss::copy(coeff_, v, Num); 
    }

    //! Resize to new size (values are preserved)
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    void resize(const ordinal_type& sz) {}

    //! Reset storage to given array, size, and stride
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz, 
		      const ordinal_type& stride, bool owned) {}

    //! Return size
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    static ordinal_type size() { return Num; }

    //! Coefficient access (avoid if possible)
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference operator[] (const ordinal_type& i) const { 
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    reference operator[] (const ordinal_type& i) { return coeff_[i]; }

    template <int i>
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    reference getCoeff() { return coeff_[i]; }

    template <int i>
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    const_reference getCoeff() const { return coeff_[i]; }

    //! Get coefficients
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
    pointer coeff() { return coeff_; }

  private:

    //! Coefficient values
    value_type coeff_[Num];

  };

}

#endif
