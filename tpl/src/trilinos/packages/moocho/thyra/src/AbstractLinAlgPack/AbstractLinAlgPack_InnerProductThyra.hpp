// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef ALAP_INNER_PRODUCT_Thyra_H
#define ALAP_INNER_PRODUCT_Thyra_H

#include "AbstractLinAlgPack_InnerProduct.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace AbstractLinAlgPack {


/** \brief Implements the inner product using
 * <tt>Thyra::VectorSpaceBase::scalarProd()</tt>.
 */
class InnerProductThyra : public InnerProduct {
public:

  /** @name Constructors / Initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  InnerProductThyra();

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  InnerProductThyra( const RCP<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc );

  /** \brief Initalize given a smart pointer to a <tt>Thyra::VetorSpace</tt> object.
   *
   * \param thyra_vec_spc [in] Smart pointer to Thyra vector
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_vec_spc.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == thyra_vec_spc.get()</tt>
   * </ul>
   */
  void initialize( const RCP<const Thyra::VectorSpaceBase<value_type> > &thyra_vec_spc );

  /** \brief Set to uninitialized and return smart pointer to the internal
   * <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == NULL</tt>
   * </ul>
   */
  RCP<const Thyra::VectorSpaceBase<value_type> > set_uninitialized();

  /** \brief Return a (converted) smart pointer to the internal smart pointer
   * to the <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * If <tt>this->thyra_vec_spc().count() == 1</tt>, then <tt>this</tt>
   * has sole ownership of the <tt>*this->thyra_vec_spc()</tt> object.
   */
  const RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc() const;

  //@}

  /** @name Overridden from InnerProduct */
  //@{

  /** \brief . */
  value_type inner_prod(const Vector& v1, const Vector& v2) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  const Thyra::VectorSpaceBase<value_type> *thyra_vector_space;
#else
  RCP<const Thyra::VectorSpaceBase<value_type> > thyra_vec_spc_;
#endif

}; // end class InnerProductThyra


// ///////////////////////////////
// Inline functions


inline
const RCP<const Thyra::VectorSpaceBase<value_type> >&
InnerProductThyra::thyra_vec_spc() const
{
  return thyra_vec_spc_;
}


} // end namespace AbstractLinAlgPack


#endif  // ALAP_INNER_PRODUCT_Thyra_H
