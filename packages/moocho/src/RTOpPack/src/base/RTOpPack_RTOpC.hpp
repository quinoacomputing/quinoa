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

// ///////////////////////////////
// RTOpPack_RTOpC.hpp

#ifndef RTOPPACK_RTOP_NEW_C_HPP
#define RTOPPACK_RTOP_NEW_C_HPP

#include "RTOpPack_OldTypes.hpp"
#include "RTOpPack_RTOpT.hpp"
#include "RTOp.h"
#include "Teuchos_dyn_cast.hpp"

namespace RTOpPack {

/** \brief Adapter subclass that uses a <tt>RTOp_RTOp</tt> object.
 *
 * ToDo: Finish documentation!
 */
class RTOpC : public RTOpT<RTOp_value_type> {
public:

  /** \brief . */
  typedef RTOp_value_type Scalar;
  /** \brief . */
  RTOpC();
  /** \brief . */
  ~RTOpC();
  /** \brief . */
  RTOp_RTOp& op();
  /** \brief . */
  const RTOp_RTOp& op() const;
  /** \brief . */
  RTOp_ReductTarget& operator()(ReductTarget& reduct_obj) const;
  /** \brief . */
  const RTOp_ReductTarget& operator()(const ReductTarget& reduct_obj) const;

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void get_reduct_type_num_entries_impl(
    const Teuchos::Ptr<int> &num_values,
    const Teuchos::Ptr<int> &num_indexes,
    const Teuchos::Ptr<int> &num_chars
    ) const;
  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create_impl() const;
  /** \brief . */
  void reduce_reduct_objs_impl(
    const ReductTarget &in_reduct_obj,
    const Teuchos::Ptr<ReductTarget> &inout_reduct_obj
    ) const;
  /** \brief . */
  void reduct_obj_reinit_impl(
    const Teuchos::Ptr<ReductTarget> &reduct_obj
    ) const;
  /** \brief . */
  void extract_reduct_obj_state_impl(
    const ReductTarget &reduct_obj,
    const Teuchos::ArrayView<primitive_value_type> &value_data,
    const Teuchos::ArrayView<index_type> &index_data,
    const Teuchos::ArrayView<char_type> &char_data
    ) const;
  /** \brief . */
  void load_reduct_obj_state_impl(
    const Teuchos::ArrayView<const primitive_value_type> &value_data,
    const Teuchos::ArrayView<const index_type> &index_data,
    const Teuchos::ArrayView<const char_type> &char_data,
    const Teuchos::Ptr<ReductTarget> &reduct_obj
    ) const;
  /** \brief . */
  bool coord_invariant_impl() const;
  /** \brief . */
  std::string op_name_impl() const;
  /** \brief . */
  void apply_op_impl(
    const Teuchos::ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const Teuchos::ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Teuchos::Ptr<ReductTarget> &reduct_obj
    ) const;

  //@}

private:

  RTOp_RTOp op_;

}; // class RTOpC

/** \brief Adapter subclass for <tt>RTOp_ReductTarget</tt>
 */
class ReductTargetC : public ReductTarget {
public:
  inline ReductTargetC( const RTOp_RTOp& op, RTOp_ReductTarget obj );
  inline ~ReductTargetC();
  inline RTOp_ReductTarget& obj();
  inline const RTOp_ReductTarget& obj() const;
private:
  const RTOp_RTOp      &op_;
  RTOp_ReductTarget    obj_;
  ReductTargetC(); // Not defined and not to be called
};

// ///////////////////////////////
// Inline member functions

// RTOpC

inline
RTOp_RTOp& RTOpC::op()
{
  return op_;
}

inline
const RTOp_RTOp& RTOpC::op() const
{
  return op_;
}

inline
RTOp_ReductTarget&
RTOpC::operator()(ReductTarget& reduct_obj) const
{
  return Teuchos::dyn_cast<ReductTargetC>(reduct_obj).obj();
}

inline
const RTOp_ReductTarget&
RTOpC::operator()(const ReductTarget& reduct_obj) const
{
  return Teuchos::dyn_cast<const ReductTargetC>(reduct_obj).obj();
}

// ReductTargetC

inline
ReductTargetC::ReductTargetC( const RTOp_RTOp& op, RTOp_ReductTarget obj )
  : op_(op), obj_(obj)
{} 

inline
ReductTargetC::~ReductTargetC()
{
  if( obj() != RTOp_REDUCT_OBJ_NULL ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      0!=RTOp_reduct_obj_free(&op_,&obj_)
      ,UnknownError
      ,"RTOpC::reduct_obj_free(...): Error, "
      "RTOp_reduct_obj_free(...) returned != 0"
      );
  }
} 

inline
RTOp_ReductTarget& ReductTargetC::obj()
{
  return obj_;
}

inline
const RTOp_ReductTarget& ReductTargetC::obj() const
{
  return obj_;
}

} // namespace RTOpPack

#endif // RTOPPACK_RTOP_NEW_C_HPP
