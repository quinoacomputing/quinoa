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

#ifndef STOKHOS_PRODUCT_LANCZOS_PCE_BASIS_HPP
#define STOKHOS_PRODUCT_LANCZOS_PCE_BASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_ReducedPCEBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief Generate a basis from a given set of PCE expansions that is 
   * orthogonal with respect to the product measure induced by these expansions.
   */
  /*!
   * Given the PCE expansions, first build a an orthogonal basis for each
   * compnent, then form the multivariate basis by a total-order tensor 
   * product.  The resulting basis is not necessarily orthogonal with respect
   * to the full measure.
   */
  template <typename ordinal_type, typename value_type>
  class ProductLanczosPCEBasis : 
    public ProductBasis<ordinal_type,value_type>,
    public ReducedPCEBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param pce polynomial chaos expansions defining new measure
     * \param Cijk sparse triple product tensor for basis defining pce
     */
    ProductLanczosPCEBasis(
     ordinal_type p,
     const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
     const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
     const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
     const Teuchos::ParameterList& params = Teuchos::ParameterList());

    //! Destructor
    virtual ~ProductLanczosPCEBasis();

    //! \name Implementation of Stokhos::OrthogPolyBasis methods
    //@{

    //! Return order of basis
    ordinal_type order() const;

    //! Return dimension of basis
    ordinal_type dimension() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\Psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is size()-1.
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j,k=0,\dots,P\f$ where
     * \f$P\f$ is size()-1.
     */
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor() const;

    //! Compute linear triple product tensor where k = 0,1,..,d
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeLinearTripleProductTensor() const;

    //! Evaluate basis polynomial \c i at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to size size()-1.
     */
    virtual void evaluateBases(
      const Teuchos::ArrayView<const value_type>& point,
      Teuchos::Array<value_type>& basis_vals) const;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const;

    //! Return string name of basis
    virtual const std::string& getName() const;

    //@}

    //! \name Implementation of Stokhos::ProductBasis methods
    //@{

   //! Get orders of each coordinate polynomial given an index \c i
    /*!
     * The returned array is of size \f$d\f$, where \f$d\f$ is the dimension of
     * the basis, and entry \f$l\f$ is given by \f$i_l\f$ where
     * \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual const MultiIndex<ordinal_type>& term(ordinal_type i) const;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type index(const MultiIndex<ordinal_type>& term) const;

    //! Return coordinate bases
    /*!
     * Array is of size dimension().
     */
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, 
							   value_type> > > 
    getCoordinateBases() const;

    //! Return maximum order allowable for each coordinate basis
    virtual MultiIndex<ordinal_type> getMaxOrders() const;

    //@}

    //! \name Implementation of Stokhos::ReducedPCEBasis methods
    //@{

    //! Transform coefficients to original basis from this basis
    virtual void 
    transformToOriginalBasis(const value_type *in, 
			     value_type *out,
			     ordinal_type ncol = 1, 
			     bool transpose = false) const;

    //! Transform coefficients from original basis to this basis
    virtual void 
    transformFromOriginalBasis(const value_type *in, 
			       value_type *out,
			       ordinal_type ncol = 1, 
			       bool transpose = false) const;

    //! Get reduced quadrature object
    virtual Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
    getReducedQuadrature() const;

    //@}

  protected:

    // Determine if a pce is linear, in that it has a total degree of at
    // most 1.  If the pce is nonlinear, return -2, or if it is constant, 
    // return -1, otherwise return the index of the variable the pce is
    // linear in, ie, if return value is i, the pce = a_{i+1}*\xi_i
    ordinal_type
    isInvariant(const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce) const;

  private:

    // Prohibit copying
    ProductLanczosPCEBasis(const ProductLanczosPCEBasis&);

    // Prohibit Assignment
    ProductLanczosPCEBasis& operator=(const ProductLanczosPCEBasis&);
    
  protected:

    typedef Teuchos::SerialDenseVector<ordinal_type,value_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> SDM;

    //! Name of basis
    std::string name;

    //! Algorithm parameters
    Teuchos::ParameterList params;

    //! Product Lanczos basis
    Teuchos::RCP< Stokhos::CompletePolynomialBasis<ordinal_type,value_type> > tensor_lanczos_basis;

    //! Values of transformed basis at quadrature points
    SDM Phi;

    //! Transition matrix from reduced basis to original
    SDM A;

    //! Projection matrix from original matrix to reduced
    SDM Ainv;

    //! Reduced quadrature object
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > reduced_quad;

    //! Temporary pce used in invariant subspace calculations
    mutable Stokhos::OrthogPolyApprox<ordinal_type, value_type> tmp_pce;

  }; // class ProductLanczosPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_ProductLanczosPCEBasisImp.hpp"

#endif
