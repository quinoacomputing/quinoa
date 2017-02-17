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

#ifndef SACADO_MPL_FOR_EACH_HPP
#define SACADO_MPL_FOR_EACH_HPP

#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_deref.hpp"

namespace Sacado {

  namespace mpl {

    template <class Seq, 
	      class Iter1 = typename mpl::begin<Seq>::type, 
	      class Iter2 = typename mpl::end<Seq>::type>
    struct for_each {
      template <typename Op>
      for_each(const Op& op) {
	op(typename mpl::deref<Iter1>::type());
	for_each<Seq, typename mpl::next<Iter1>::type, Iter2> f(op);
      }
    };

    template <class Seq, class Iter1>
    struct for_each<Seq, Iter1, Iter1> {
      template <typename Op>
      for_each(const Op& op) {}
    };

  }

}

#endif // SACADO_MPL_FOR_EACH_HPP
