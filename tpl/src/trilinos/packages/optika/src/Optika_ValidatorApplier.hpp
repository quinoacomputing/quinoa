// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#ifndef OPTIKA_VALIDATORAPPLIER_HPP_
#define OPTIKA_VALIDATORAPPLIER_HPP_
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QIntValidator>
#include <QDoubleValidator>
#include "Optika_ConfigDefs.hpp"
#include <QLineEdit>

/*! \file Optika_ValidatorApplier.hpp
    \brief A collection of objects which
    apply the restriction of a EnhancedNumberValidator
    to a SpinBoxes and LineEdits.
*/

namespace Optika{

/**
 * \brief A templated class that applies the minimum,
 * maximum and step specified in an EnhancedNumberValidator
 * onto a QSpingBox.
 */
template<class S> class ValidatorApplier{
public:
  //! @name Applier Functions
  //@{

  /**
   * \brief Applied attributes of the validator to the spin box.
   *
   * @param validator Validator whose attributes are to be applied to
   * the spinbox
   * @param spinBox A point to the spinbox upon which the validator's
   * attributes should be applied.
   */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<S> > validator, QSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<S>::min());
			spinBox->setMaximum(EnhancedNumberTraits<S>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<S>::defaultStep());
		}
	}

	static void applyToLineEdit(RCP<const EnhancedNumberValidator<S> > validator, QLineEdit *lineEdit){
    QIntValidator* qvalidator = new QIntValidator(lineEdit);
		if(!is_null(validator)){
      qvalidator->setRange(validator->getMin(), validator->getMax());
		}
		else{
      qvalidator->setRange(EnhancedNumberTraits<S>::min(), EnhancedNumberTraits<S>::max());
		}
    lineEdit->setValidator(qvalidator);
  }

  //@}

};

/**
 * \brief A template spcialization of the ValidatorApplier
 * class on the type double.
 */
template<>
class ValidatorApplier<double>{
public:

  //! @name Applier Functions
  //@{

	/**
	 * \brief Applies an EnhancedNumberValidator of type double to a QDoubleSpinBox
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<double> > validator, QDoubleSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<double>::min());
			spinBox->setMaximum(EnhancedNumberTraits<double>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<double>::defaultStep());
			spinBox->setDecimals(EnhancedNumberTraits<double>::defaultPrecision());
		}
	}

	static void applyToLineEdit(RCP<const EnhancedNumberValidator<double> > validator, QLineEdit *lineEdit){
    QDoubleValidator* qvalidator = new QDoubleValidator(lineEdit);
		if(!is_null(validator)){
      qvalidator->setRange(validator->getMin(), validator->getMax());
      qvalidator->setDecimals(validator->getPrecision());
		}
		else{
      qvalidator->setRange(EnhancedNumberTraits<double>::min(), EnhancedNumberTraits<double>::max());
      qvalidator->setDecimals(EnhancedNumberTraits<double>::defaultPrecision());
		}
    lineEdit->setValidator(qvalidator);
  }

  //@}
};

/**
 * \brief A template specialzation of the ValidatorApplier
 * class on the type float.
 */
template<>
class ValidatorApplier<float>{
public:

  //! @name Applier Functions
  //@{

	/**
	 * \brief Applies an EnhancedNumberValidator of type float to a QDoubleSpinBox.
	 *
	 * @param validator The validator to be useed.
	 * @param spinBox The SpinBox on which to apply the validator.
	 */
	static void applyToSpinBox(RCP<const EnhancedNumberValidator<float> > validator, QDoubleSpinBox *spinBox){
		if(!is_null(validator)){
			spinBox->setMinimum(validator->getMin());
			spinBox->setMaximum(validator->getMax());
			spinBox->setSingleStep(validator->getStep());
			spinBox->setDecimals(validator->getPrecision());
		}
		else{
			spinBox->setMinimum(EnhancedNumberTraits<float>::min());
			spinBox->setMaximum(EnhancedNumberTraits<float>::max());
			spinBox->setSingleStep(EnhancedNumberTraits<float>::defaultStep());
			spinBox->setDecimals(EnhancedNumberTraits<float>::defaultPrecision());
		}
	}

	static void applyToLineEdit(RCP<const EnhancedNumberValidator<float> > validator, QLineEdit *lineEdit){
    QDoubleValidator* qvalidator = new QDoubleValidator(lineEdit);
		if(!is_null(validator)){
      qvalidator->setRange(validator->getMin(), validator->getMax());
      qvalidator->setDecimals(validator->getPrecision());
		}
		else{
      qvalidator->setRange(EnhancedNumberTraits<float>::min(), EnhancedNumberTraits<float>::max());
      qvalidator->setDecimals(EnhancedNumberTraits<float>::defaultPrecision());
		}
    lineEdit->setValidator(qvalidator);
  }

  //@}
}; 

}
#endif // OPTIKA_VALIDATORAPPLIER_HPP_
