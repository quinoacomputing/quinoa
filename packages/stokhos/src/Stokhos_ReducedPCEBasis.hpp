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

#ifndef STOKHOS_REDUCED_PCE_BASIS_HPP
#define STOKHOS_REDUCED_PCE_BASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"

namespace Stokhos {

  /*! 
   * \brief Abstract base class for reduced basis strategies built from 
   * polynomial chaos expansions in some other basis.
   */
  template <typename ordinal_type, typename value_type>
  class ReducedPCEBasis : 
    public virtual OrthogPolyBasis<ordinal_type,value_type> {
  public:

    //! Default constructor
    ReducedPCEBasis() {}

    //! Destructor
    virtual ~ReducedPCEBasis() {}

    //! \name ReducedBasis virtual methods
    //@{

    //! Transform coefficients to original basis from this basis
    virtual void 
    transformToOriginalBasis(const value_type *in, 
			     value_type *out,
			     ordinal_type ncol = 1, 
			     bool transpose = false) const = 0;

    //! Transform coefficients from original basis to this basis
    virtual void 
    transformFromOriginalBasis(const value_type *in, 
			       value_type *out,
			       ordinal_type ncol = 1, 
			       bool transpose = false) const = 0;

    //! Get reduced quadrature object
    virtual Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
    getReducedQuadrature() const = 0;

    //@}

  private:

    // Prohibit copying
    ReducedPCEBasis(const ReducedPCEBasis&);

    // Prohibit Assignment
    ReducedPCEBasis& operator=(const ReducedPCEBasis&);

  }; // class ReducedPCEBasis

} // Namespace Stokhos

#endif
