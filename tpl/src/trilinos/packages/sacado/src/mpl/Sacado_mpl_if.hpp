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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_MPL_IF_HPP
#define SACADO_MPL_IF_HPP

#include "Sacado_mpl_type_wrap.hpp"

namespace Sacado {

  namespace mpl {

    template <bool cond, class T1, class T2> struct mpl_if_c {};
    template <class T1, class T2> struct mpl_if_c<true,T1,T2> :
      type_wrap<T1> {};
    template <class T1, class T2> struct mpl_if_c<false,T1,T2> :
      type_wrap<T2> {};

    template <class C, class T1, class T2> struct mpl_if :
      mpl_if_c<C::value,T1,T2> {};

  }

}

#endif // SACADO_MPL_IF_HPP
