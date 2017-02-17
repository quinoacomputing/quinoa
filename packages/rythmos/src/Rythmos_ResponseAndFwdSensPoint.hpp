//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_RESPONSE_AND_FWD_SEND_POINT_HPP
#define RYTHMOS_RESPONSE_AND_FWD_SEND_POINT_HPP


#include "Rythmos_Types.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_AssertOp.hpp"


namespace Rythmos {


/** \brief Simple class to combine a response and it's forward sensitivity at
 * a time point.
 *
 * NOTE: Compiler-generated copy constructor and assignment operator functions
 * are allowed and implement shallow copy.
 */
template<class Scalar>
class ResponseAndFwdSensPoint {
public:

  /** \brief . */
  ResponseAndFwdSensPoint()
    : t_(ScalarTraits<Scalar>::zero())
    {}

  // 2007/11/19: rabartl: ToDo: Add constructors that take only t and g or t
  // and DgDp if needed!
  
  /** \brief . */
  ResponseAndFwdSensPoint(
    const Scalar &t,
    const RCP<const Thyra::VectorBase<Scalar> > &g,
    const RCP<const Thyra::MultiVectorBase<Scalar> > &DgDp
    )
    :t_(t), g_(g), DgDp_(DgDp)
    {
#ifdef HAVE_RYTHMOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(is_null(g));
      TEUCHOS_TEST_FOR_EXCEPT(is_null(DgDp));
      THYRA_ASSERT_VEC_SPACES("Rythmos::ResponseAndFwdSensPoint()",
        *g->space(), * DgDp->range() );
#endif      
    }

  /** \brief . */
  Scalar t() const
    { return t_; }

  /** \brief . */
  const RCP<const Thyra::VectorBase<Scalar> > g() const
    { return g_; }

  /** \brief . */
  const RCP<const Thyra::MultiVectorBase<Scalar> > DgDp() const
    { return DgDp_; }

private:
  
  Scalar t_;
  RCP<const Thyra::VectorBase<Scalar> > g_;
  RCP<const Thyra::MultiVectorBase<Scalar> > DgDp_;

};


} // namespace Rythmos


#endif //RYTHMOS_RESPONSE_AND_FWD_SEND_POINT_HPP
