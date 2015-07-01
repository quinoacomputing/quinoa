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
#ifndef OPTIKA_ARRAYHELPERFUNCTIONS_HPP_
#define OPTIKA_ARRAYHELPERFUNCTIONS_HPP_
#include <QStringList>
#include <QVariant>
#include "Optika_Types.hpp"
#include "Optika_ConfigDefs.hpp"
/*! \file Opitka_ArrayHelperFunctions.hpp
    \brief Helper functions used for manipulating
    and querying Arrays
*/
namespace Optika{

/**
 * \brief Determines the type of array stored in a parameter.
 *
 * @param parameter The parameter whose array type is in question.
 * @return A QString containing the type of array in the parameter.
 */
QString determineArrayType(RCP<const ParameterEntry> parameter, bool twoD=false);

/**
 * \brief Creates a QVariant containing the array that is in
 * arrayEntry.
 *
 * @param arrayEntry The parameter entry containing the array.
 * @param type The array's template type.
 * @return A QVariant containing an Array that is equal to the Array in array entry.
 */
QVariant arrayEntryToVariant(
  RCP<const ParameterEntry> arrayEntry, QString type, bool twoD=false);

/**
 * \brief Given a type string, determines the template type of the Array.
 *
 * @param itemType The type string describing the array.
 * @return The template type of the Array.
 */
QString getArrayType(QString itemType);

/**
 * \brief Determines wether or no the array inside a ParameterEntry is empty.
 *
 * @param arrayEntry The parameter entry containging the array.
 * @param type The template type of the array.
 * @return True if the array is empty, false otherwise.
 */
bool isArrayEmpty(RCP<const ParameterEntry> arrayEntry, QString type);

}
#endif /* OPTIKA_ARRAYHELPERFUNCTIONS_HPP_ */
