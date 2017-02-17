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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_QUADRATURE
#define STOKHOS_QUADRATURE

#include <ostream>
#include "Teuchos_Array.hpp"

namespace Stokhos {

  //! Abstract base class for quadrature methods
  template <typename ordinal_type, typename value_type>
  class Quadrature {
  public:

    //! Constructor
    Quadrature() {}

    //! Destructor
    virtual ~Quadrature() {}

    //! Get number of quadrature points
    virtual ordinal_type size() const = 0;

    //! Get quadrature points
    /*!
     * Array is dimensioned Q-by-d where Q is the number of quadrature
     * points and d is the dimension of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> >& 
    getQuadPoints() const = 0;

    //! Get quadrature weights
    /*!
     * Array is of size Q where Q is the number of quadrature points.
     */
    virtual const Teuchos::Array<value_type>& 
    getQuadWeights() const = 0;

    //! Get values of basis at quadrature points
    /*!
     * Array is dimensioned Q-by-P where Q is the number of quadrature
     * points and P is the size of the basis.
     */
    virtual const Teuchos::Array< Teuchos::Array<value_type> > & 
    getBasisAtQuadPoints() const = 0;

    //! Print quadrature data
    virtual std::ostream& print(std::ostream& os) const = 0;

  private:

    // Prohibit copying
    Quadrature(const Quadrature&);

    // Prohibit Assignment
    Quadrature& operator=(const Quadrature& b);

  }; // class Quadrature

  //! Print quadrature object to stream
  template <typename ordinal_type, typename value_type>
  std::ostream& operator << (std::ostream& os,
			     const Quadrature<ordinal_type, value_type>& quad) {
    return quad.print(os);
  }

} // namespace Stokhos

#endif // STOKHOS_QUADRATURE
