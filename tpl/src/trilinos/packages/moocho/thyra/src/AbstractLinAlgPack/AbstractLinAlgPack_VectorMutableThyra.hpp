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

#ifndef ALAP_VECTOR_MUTABLE_Thyra_HPP
#define ALAP_VECTOR_MUTABLE_Thyra_HPP

#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorApplyOpSerialBase.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "Thyra_VectorBase.hpp"

namespace AbstractLinAlgPack {

/** \brief <tt>VectorMutable</tt> adapter subclass for <tt>Thyra::VectorBase</tt>.
 */
class VectorMutableThyra : public VectorMutable, private VectorApplyOpSerialBase {
public:

  /** @name Constructors / Initializers */
  //@{

  /** \brief Construct to uninitialized.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  VectorMutableThyra();
  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorMutableThyra( const Teuchos::RCP<Thyra::VectorBase<value_type> >& thyra_vec );
  /** \brief Initalize given a smart pointer to a <tt>Thyra::Vetor</tt> object.
   *
   * @param  thyra_vec  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
   *
   * Preconditioins:<ul>
   * <li><tt>thyra_vec.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == thyra_vec.get()</tt>
   * </ul>
   */
  void initialize( const Teuchos::RCP<Thyra::VectorBase<value_type> >& thyra_vec );
  /** \brief Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorBase</tt> object.
   *
   * Postconditioins:<ul>
   * <li><tt>this->thyra_vec().get() == NULL</tt>
   * </ul>
   */
  Teuchos::RCP<Thyra::VectorBase<value_type> > set_uninitialized();
  /** \brief Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorBase</tt> object.
   *
   * If <tt>this->thyra_vec().count() == 2</tt>, then <tt>this</tt>
   * has so ownership of the <tt>*this->thyra_vec()</tt> object.
   */
  Teuchos::RCP<const Thyra::VectorBase<value_type> > thyra_vec() const;

  //@}

  /** @name Methods overridden from Vector */
  //@{

  /** \brief . */
  const VectorSpace& space() const;
  /** \brief . */
  void apply_op(
    const RTOpPack::RTOp       &op
    ,const size_t              num_vecs
    ,const Vector*             vecs[]
    ,const size_t              num_targ_vecs
    ,VectorMutable*            targ_vecs[]
    ,RTOpPack::ReductTarget    *reduct_obj
    ,const index_type          first_ele
    ,const index_type          sub_dim
    ,const index_type          global_offset
    ) const;
  /** \brief . */
  index_type dim() const;
  /** \brief . */
  void get_sub_vector( const Range1D& rng, RTOpPack::SubVector* sub_vec ) const;
  /** \brief . */
  void free_sub_vector( RTOpPack::SubVector* sub_vec ) const;

  //@}

  /** @name Methods overridden from VectorMutable */
  //@{

  /** \brief . */
  void get_sub_vector( const Range1D& rng, RTOpPack::MutableSubVector* sub_vec );
  /** \brief . */
  void commit_sub_vector( RTOpPack::MutableSubVector* sub_vec );
  /** \brief . */
  void set_sub_vector( const RTOpPack::SparseSubVector& sub_vec );
  /** \brief . */
  void Vp_StMtV(
    value_type                       alpha
    ,const GenPermMatrixSlice        &P
    ,BLAS_Cpp::Transp                P_trans
    ,const Vector                    &x
    ,value_type                      beta
    );


  //@}

private:

#ifdef DOXYGEN_COMPILE
  Thyra::VectorBase<value_type>                          *thyra_vector;
#else
  Teuchos::RCP<Thyra::VectorBase<value_type> >   thyra_vec_;
#endif
  VectorSpaceThyra                                       space_;

  // Not defined and not to be called!
  //VectorMutableThyra(const VectorMutableThyra&);
  VectorMutableThyra& operator=(const VectorMutableThyra&);

}; // class VectorMutableThyra 

// ////////////////////////////////
// Inline functions

inline
Teuchos::RCP<const Thyra::VectorBase<value_type> >
VectorMutableThyra::thyra_vec() const
{
  return thyra_vec_;
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_MUTABLE_Thyra_HPP
