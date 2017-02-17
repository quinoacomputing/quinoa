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

#ifndef STOKHOS_ONEDORTHOGPOLYBASIS_HPP
#define STOKHOS_ONEDORTHOGPOLYBASIS_HPP

#include <ostream>
#include <string>
#include "Stokhos_Dense3Tensor.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_DAKOTA
#include "sandia_rules.hpp"
#endif

//! Top-level namespace for Stokhos classes and functions.
namespace Stokhos {

  //! Abstract base class for 1-D orthogonal polynomials.
  /*!
   * This class provides an abstract interface for univariate orthogonal
   * polynomials.  Orthogonality is defined by the inner product
   * \f[
   *      (f,g) = \langle fg \rangle =
   *              \int_{-\infty}^{\infty} f(x)g(x) \rho(x) dx
   * \f]
   * where \f$\rho\f$ is the density function of the measure associated with
   * the orthogonal polynomials.
   * See Stokhos::RecurrenceBasis for a general implementation
   * of this interface based on the three-term recurrence satisfied by
   * these polynomials.  Multivariate polynomials can be formed from
   * a collection of univariate polynomials through tensor products (see
   * Stokhos::CompletePolynomialBasis).
   *
   * Like most classes in Stokhos, the class is templated on the ordinal
   * and value types.  Typically \c ordinal_type = \c int and \c value_type
   * = \c double.
   */
  template <typename ordinal_type, typename value_type>
  class OneDOrthogPolyBasis {
  public:

    //! Default constructor
    OneDOrthogPolyBasis() {};

    //! Destructor
    virtual ~OneDOrthogPolyBasis() {};

    //! Return order of basis (largest monomial degree \f$P\f$).
    virtual ordinal_type order() const = 0;

    //! Return total size of basis (given by order() + 1).
    virtual ordinal_type size() const = 0;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is given by order().
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const = 0;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const = 0;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is size()-1 and \f$k=0,\dots,p\f$ where \f$p\f$
     * is the supplied \c order.
     */
    virtual
    Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> >
    computeTripleProductTensor() const = 0;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is size()-1 and \f$k=0,\dots,p\f$ where \f$p\f$
     * is the supplied \c order.
     */
    virtual
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
    computeSparseTripleProductTensor(ordinal_type order) const = 0;

    //! Compute derivative double product tensor
    /*!
     * The \f$(i,j)\f$ entry of the tensor \f$B_{ij}\f$ is given by
     * \f$B_{ij} = \langle\psi_i'\psi_j\rangle\f$ where \f$\psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is the order of the basis.
     */
    virtual
    Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
    computeDerivDoubleProductTensor() const = 0;

    //! Evaluate each basis polynomial at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to order order().
     */
    virtual void evaluateBases(const value_type& point,
                               Teuchos::Array<value_type>& basis_pts) const = 0;

    /*!
     * \brief Evaluate basis polynomial given by order \c order at given
     * point \c point.
     */
    virtual value_type evaluate(const value_type& point,
                                ordinal_type order) const = 0;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const {};

    //! Return string name of basis
    virtual const std::string& getName() const = 0;

    /*!
     * \brief Compute quadrature points, weights, and values of
     * basis polynomials at given set of points \c points.
     */
    /*!
     * \c quad_order specifies the order to which the quadrature should be
     * accurate, not the number of quadrature points.  The number of points
     * is given by (\c quad_order + 1) / 2.   Note however the passed arrays
     * do NOT need to be sized correctly on input as they will be resized
     * appropriately.
     */
    virtual void
    getQuadPoints(ordinal_type quad_order,
                  Teuchos::Array<value_type>& points,
                  Teuchos::Array<value_type>& weights,
                  Teuchos::Array< Teuchos::Array<value_type> >& values) const = 0;

    /*!
     * Return polynomial degree of exactness for a given number of quadrature
     * points.
     */
    virtual ordinal_type quadDegreeOfExactness(ordinal_type n) const = 0;

    /*!
     * \brief Clone this object with the option of building a higher order
     * basis.
     */
    /*!
     * This method is following the Prototype pattern (see Design Pattern's textbook).
     * The slight variation is that it allows the order of the polynomial to be modified,
     * otherwise an exact copy is formed. The use case for this is creating basis functions
     * for column indices in a spatially varying adaptive refinement context.
     */
    virtual Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > cloneWithOrder(ordinal_type p) const = 0;

    //! Evaluate coefficient growth rule for Smolyak-type bases
    virtual ordinal_type coefficientGrowth(ordinal_type n) const = 0;

    //! Evaluate point growth rule for Smolyak-type bases
    virtual ordinal_type pointGrowth(ordinal_type n) const = 0;

    //! Function pointer needed for level_to_order mappings
    typedef int ( *LevelToOrderFnPtr ) ( int level, int growth );

    //! Get sparse grid level_to_order mapping function
    /*!
     * Predefined functions are:
     *  webbur::level_to_order_linear_wn Symmetric Gaussian linear growth
     *  webbur::level_to_order_linear_nn Asymmetric Gaussian linear growth
     *  webbur::level_to_order_exp_cc    Clenshaw-Curtis exponential growth
     *  webbur::level_to_order_exp_gp    Gauss-Patterson exponential growth
     *  webbur::level_to_order_exp_hgk   Genz-Keister exponential growth
     *  webbur::level_to_order_exp_f2    Fejer-2 exponential growth
     */
    virtual LevelToOrderFnPtr getSparseGridGrowthRule() const = 0;

    //! Set sparse grid rule
    virtual void setSparseGridGrowthRule(LevelToOrderFnPtr ptr) = 0;

  private:

    // Prohibit copying
    OneDOrthogPolyBasis(const OneDOrthogPolyBasis&);

    // Prohibit Assignment
    OneDOrthogPolyBasis& operator=(const OneDOrthogPolyBasis& b);


  }; // class OrthogPolyBasis

  //! Print basis to stream \c os.
  template <typename ordinal_type, typename value_type>
  std::ostream&
  operator << (std::ostream& os,
               const OneDOrthogPolyBasis<ordinal_type, value_type>& b) {
    b.print(os);
    return os;
  }

} // Namespace Stokhos

#endif
