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

#ifndef ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP
#define ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP

#include "AbstractLinAlgPack_VectorSpaceFactory.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"

namespace AbstractLinAlgPack {

/** \brief <tt>VectorSpaceFactory</tt> adapter subclass for <tt>Thyra::VectorSpaceBase</tt>.
 */
class VectorSpaceFactoryThyra : public VectorSpaceFactory {
public:

  /** @name Constructors / Initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  VectorSpaceFactoryThyra();
  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceFactoryThyra( const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> > &thyra_vec_spc_fcty );
  /** \brief Initalize given a smart pointer to a <tt>Thyra::VetorSpaceFactory</tt> object.
   *
   * @param  thyra_vec_spc_fcty  [in] Smart pointer to Thyra vector
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_vec_spc_fcty.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc_fcty().get() == thyra_vec_spc_fcty.get()</tt>
   * </ul>
   */
  void initialize( const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> > &thyra_vec_spc_fcty );
  /** \brief Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorSpaceBase</tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc_fcty().get() == NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> > set_uninitialized();
  /** \brief Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorSpaceBase</tt> object.
   *
   * If <tt>this->thyra_vec_spc_fcty().count() == 1</tt>, then <tt>this</tt>
   * has sole ownership of the <tt>*this->thyra_vec_spc_fcty()</tt> object.
   */
  const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> >& thyra_vec_spc_fcty() const;

  //@}

  /** @name Overridden from VectorSpaceFactory */
  //@{

  /** \brief . */
  space_ptr_t create_vec_spc(index_type dim) const;

  //@}
  
private:

#ifdef DOXYGEN_COMPILE
  const Thyra::VectorSpaceFactoryBase<value_type>                             *thyra_vector_space_factory;
#else
  Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> >  thyra_vec_spc_fcty_;
#endif

}; // end class VectorSpaceFactoryThyra

// ///////////////////////////////
// Inline functions

inline
const Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<value_type> >&
VectorSpaceFactoryThyra::thyra_vec_spc_fcty() const
{
  return thyra_vec_spc_fcty_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_SPACE_FACTORY_Thyra_HPP
