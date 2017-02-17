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

#ifndef STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP
#define STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP

#include "Stokhos_StandardStorage.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"

namespace Stokhos {

  //! Strategy interface for computing PCE of a/b
  template <typename ordinal_type, typename value_type, typename node_type> 
  class DivisionExpansionStrategy {
  public:

    //! Constructor
    DivisionExpansionStrategy() {}

    //! Destructor
    virtual ~DivisionExpansionStrategy() {}
 
    // Division operation:  c = \alpha*(a/b) + beta*c
    virtual void divide(
      Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& c,
      const value_type& alpha,
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
      const Stokhos::OrthogPolyApprox<ordinal_type, value_type, node_type>& b,
      const value_type& beta) = 0;

  private:

    // Prohibit copying
    DivisionExpansionStrategy(const DivisionExpansionStrategy&);

    // Prohibit Assignment
    DivisionExpansionStrategy& operator=(const DivisionExpansionStrategy& b);
    
  }; // class DivisionExpansionStrategy

} // namespace Stokhos

#endif // STOKHOS_DIVISION_EXPANSION_STRATEGY_HPP
