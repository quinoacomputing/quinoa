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

#ifndef ALAP_VECTOR_SPACE_Thyra_HPP
#define ALAP_VECTOR_SPACE_Thyra_HPP

#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Thyra_VectorSpaceBase.hpp"

namespace AbstractLinAlgPack {

/** \brief <tt>VectorSpace</tt> adapter subclass for <tt>Thyra::VectorSpaceBase<value_type> </tt>.
 *
 * Note that the default copy constructor and assignment operators are
 * allowed which yield in shallow copy, not deep copy.
 */
class VectorSpaceThyra : public VectorSpace {
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  VectorSpaceThyra();
  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceThyra(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
    ,const inner_prod_ptr_t                                                  &inner_prod    = Teuchos::null
    );
  /** \brief Initalize given a smart pointer to a <tt>Thyra::VetorSpace</tt> object.
   *
   * @param  thyra_vec_spc  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
   * @param  inner_prod    [in] Smart pointer to an inner product.  If <tt>inner_prod.get()==NULL</tt>
   *                       then a <tt>InnerProductThyra</tt> object will be used which will
   *                       point to this.
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_vec_spc.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == thyra_vec_spc.get()</tt>
   * <li>[<tt>inner_prod.get()!=NULL</tt>]
   *     <tt>this->inner_prod().get()==inner_prod.get()</tt>
   * <li>[<tt>inner_prod.get()==NULL</tt>]
   *     <tt>dynamic_cast<const InnerProductThyra*>(this->inner_prod().get()).thyra_vec_spc().get()==thyra_vec_spc.get()</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >    &thyra_vec_spc
    ,const inner_prod_ptr_t                                                  &inner_prod    = Teuchos::null
    );
  /** \brief Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec_spc().get() == NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> > set_uninitialized();
  /** \brief Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorSpaceBase<value_type> </tt> object.
   *
   * If <tt>this->thyra_vec_spc().count() == 1</tt>, then <tt>this</tt>
   * has sole ownership of the <tt>*this->thyra_vec_spc()</tt> object.
   */
  const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >& thyra_vec_spc() const;

  //@}

  /** @name Overridden from VectorSpace */
  //@{

  /** \brief . */
  space_ptr_t clone() const;
  /** \brief . */
  bool is_compatible(const VectorSpace& vec_spc ) const;
  /** \brief . */
  bool is_in_core() const;
  /** \brief . */
  index_type dim() const;
  /** \brief . */
  vec_mut_ptr_t create_member() const;
  /** \brief . */
  space_fcty_ptr_t small_vec_spc_fcty() const;
  /** \brief . */
  multi_vec_mut_ptr_t create_members(size_type num_vecs) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  const Thyra::VectorSpaceBase<value_type>                              *thyra_vector_space;
#else
  Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >  thyra_vec_spc_;
#endif

}; // end class VectorSpaceThyra

// ///////////////////////////////
// Inline functions

inline
const Teuchos::RCP<const Thyra::VectorSpaceBase<value_type> >&
VectorSpaceThyra::thyra_vec_spc() const
{
  return thyra_vec_spc_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_SPACE_Thyra_HPP
