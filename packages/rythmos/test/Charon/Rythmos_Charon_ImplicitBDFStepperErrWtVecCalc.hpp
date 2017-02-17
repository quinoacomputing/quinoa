//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#ifndef RYTHMOS_CHARON_IMPLICIT_BDF_ERR_WT_VEC_CALC_BASE_HPP
#define RYTHMOS_CHARON_IMPLICIT_BDF_ERR_WT_VEC_CALC_BASE_HPP

#include "Rythmos_ErrWtVecCalcBase.hpp"
#include "Teuchos_RCP.hpp"

namespace RythmosCharon {

class CharonImplicitBDFStepperErrWtVecCalc
  : virtual public Rythmos::ErrWtVecCalcBase<double>
{
  public:
  CharonImplicitBDFStepperErrWtVecCalc();
  virtual ~CharonImplicitBDFStepperErrWtVecCalc();
  void errWtVecSet(
      Thyra::VectorBase<double>* weight, 
      const Thyra::VectorBase<double>& vector, 
      double relTol, 
      double absTol
      ) const;
  // Overridden from Teuchos::ParameterListAcceptor
  void setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& paramList );
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  private:
    Teuchos::RCP<Teuchos::ParameterList> paramList_;
};


} // namespace RythmosCharon

#endif // RYTHMOS_CHARON_IMPLICIT_BDF_ERR_WT_VEC_CALC_BASE_HPP

