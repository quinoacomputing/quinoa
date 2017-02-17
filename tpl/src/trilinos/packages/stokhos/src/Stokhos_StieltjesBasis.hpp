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

#ifndef STOKHOS_STIELTJESBASIS_HPP
#define STOKHOS_STIELTJESBASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Stokhos_RecurrenceBasis.hpp"
#include "Stokhos_Quadrature.hpp"

namespace Stokhos {

  /*! 
   * \brief Generates three-term recurrence using the Discretized Stieltjes 
   * procedure applied to a functional mapping another basis.
   */
  template <typename ordinal_type, typename value_type, typename func_type>
  class StieltjesBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param func mapping defining new density function
     * \param quad quadrature data for basis of PC expansion
     * \param use_pce_quad_points whether to use quad to define quadrature
     *        points for the new basis, or whether to use the Golub-Welsch
     *        system.
     * \param normalize whether polynomials should be given unit norm
     */
    StieltjesBasis(
      ordinal_type p,
      const Teuchos::RCP<const func_type>& func,
      const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
      bool use_pce_quad_points,
      bool normalize = false);

    //! Destructor
    ~StieltjesBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

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
    virtual Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > cloneWithOrder(ordinal_type p) const;

    //@}

  protected:

    //! \name Implementation of Stokhos::RecurrenceBasis methods
    //@{ 

    //! Compute recurrence coefficients
    virtual bool
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta,
				  Teuchos::Array<value_type>& gamma) const;

    //! Setup basis after computing recurrence coefficients
    virtual void setup();

    //@}

    //! Compute 3-term recurrence using Stieljtes procedure
    void stieltjes(ordinal_type nstart,
		   ordinal_type nfinish,
		   const Teuchos::Array<value_type>& weights,
		   const Teuchos::Array<value_type>& points,
		   Teuchos::Array<value_type>& a,
		   Teuchos::Array<value_type>& b,
		   Teuchos::Array<value_type>& nrm,
		   Teuchos::Array< Teuchos::Array<value_type> >& phi_vals) const;

    /*! 
     * \brief Compute \f$\int\psi^2_k(t) d\lambda(t)\f$ and 
     * \f$\int t\psi^2_k(t) d\lambda(t)\f$
     */
    void integrateBasisSquared(
      ordinal_type k, 
      const Teuchos::Array<value_type>& a,
      const Teuchos::Array<value_type>& b,
      const Teuchos::Array<value_type>& weights,
      const Teuchos::Array<value_type>& points,
      Teuchos::Array< Teuchos::Array<value_type> >& phi_vals,
      value_type& val1, value_type& val2) const;
    
    //! Evaluate polynomials via 3-term recurrence
    void evaluateRecurrence(ordinal_type k,
			    const Teuchos::Array<value_type>& a,
			    const Teuchos::Array<value_type>& b,
			    const Teuchos::Array<value_type>& points,
			    Teuchos::Array< Teuchos::Array<value_type> >& values) const;

    //! Copy constructor with specified order
    StieltjesBasis(ordinal_type p, const StieltjesBasis& basis);

  private:

    // Prohibit copying
    StieltjesBasis(const StieltjesBasis&);

    // Prohibit Assignment
    StieltjesBasis& operator=(const StieltjesBasis& b);

  protected:

    //! PC expansion
    Teuchos::RCP<const func_type> func;

    //! Quadrature object
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > quad;

    //! PCE quadrature weights
    const Teuchos::Array<value_type>& pce_weights;

    //! Values of PCE basis functions at quadrature points
    const Teuchos::Array< Teuchos::Array<value_type> >& basis_values;

    //! Values of func at quadrature points
    mutable Teuchos::Array<value_type> func_vals;

    //! Values of generated polynomials at PCE quadrature points
    mutable Teuchos::Array< Teuchos::Array<value_type> > phi_vals;

    //! Use underlying pce's quadrature data
    bool use_pce_quad_points;

  }; // class StieltjesBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_StieltjesBasisImp.hpp"

#endif
