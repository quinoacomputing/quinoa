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

#include "Teuchos_UnitTestHarness.hpp"
#include "Optika_treeitem.hpp"
#include "Optika_ArrayHelperFunctions.hpp"
#include "Optika_ValidatorApplier.hpp"

namespace Optika{

TEUCHOS_UNIT_TEST(Array_Helper_Functions, ArrayToQVariant){
  Array<std::string> testArray = Teuchos::tuple<std::string>("blah1", "blah2");
  ParameterEntry entry(testArray);
  QVariant testVariant = arrayEntryToVariant(rcpFromRef(entry), stringId);
  Array<std::string> extractedArray = testVariant.value<Array<std::string> >();
  TEST_EQUALITY(testArray,extractedArray);
  out << "Test Array: " << testArray << std::endl;
  out << "Extracted Array: " << extractedArray << std::endl;

  TwoDArray<float> twoDtestArray(3,3,4.0);
  ParameterEntry twoDEntry(twoDtestArray);
  QVariant twoDVariant = arrayEntryToVariant(rcpFromRef(twoDEntry), floatId, true);
  TwoDArray<float> twoDExtracted = twoDVariant.value<TwoDArray<float> >();
  TEST_EQUALITY(twoDExtracted, twoDtestArray)
  out << "TwoD Test: " << twoDtestArray << std::endl;
  out << "TwoD Extracted: " << twoDExtracted << std::endl;
  
}

TEUCHOS_UNIT_TEST(Array_Helper_Functions, DetermineArrayType){
  Array<int> intArray  = Teuchos::tuple<int>(5,6,7,8);
  ParameterEntry intParam(intArray);
  TEST_EQUALITY(
    determineArrayType(rcpFromRef(intParam)).toStdString(), 
    intId.toStdString());

  Array<double> doubleArray = 
    Teuchos::tuple<double>(5.0,6.4,7.5,8.9);
  ParameterEntry doubleParam(doubleArray);
  TEST_EQUALITY(
    determineArrayType(rcpFromRef(doubleParam)).toStdString(), 
    doubleId.toStdString());

  Array<std::string> stringArray = 
    Teuchos::tuple<std::string>("hello", "one");
  ParameterEntry stringParam(stringArray);
  TEST_EQUALITY(
    determineArrayType(rcpFromRef(stringParam)).toStdString(), 
    stringId.toStdString());
  TwoDArray<float> floatTwoDArray(4,3,5.0);
  ParameterEntry floatParam(floatTwoDArray);
  TEST_EQUALITY(
    determineArrayType(rcpFromRef(floatParam), true).toStdString(), 
    floatId.toStdString());

  TwoDArray<double> dDoubleArray(2,2,4.0);
  ParameterEntry ddoubleParam(dDoubleArray);
  TEST_EQUALITY(
    determineArrayType(rcpFromRef(ddoubleParam), true).toStdString(), 
    doubleId.toStdString());

  TwoDArray<int> dIntArray(2,2,4);
  ParameterEntry dintParam(dIntArray);
  TEST_EQUALITY(
    determineArrayType(rcpFromRef(dintParam), true).toStdString(), 
    intId.toStdString());
}

TEUCHOS_UNIT_TEST(Array_Helper_Functions, GetArrayType){
  TEST_EQUALITY(
    intId.toStdString(), 
    getArrayType(arrayId + " " + intId).toStdString());

  TEST_EQUALITY(
    doubleId.toStdString(), 
    getArrayType(arrayId + " " + doubleId).toStdString());

  TEST_EQUALITY(
    stringId.toStdString(), 
    getArrayType(arrayId + " " + stringId).toStdString());

}

TEUCHOS_UNIT_TEST(Array_Helper_Functions, IsArrayEmpty){
  Array<int> intArray;
  ParameterEntry intArrayParam(intArray);
  TEST_ASSERT(isArrayEmpty(rcpFromRef(intArrayParam), intId));

  Array<std::string> stringArray(4, "blah");
  ParameterEntry stringArrayParam(stringArray);
  TEST_ASSERT(!isArrayEmpty(rcpFromRef(stringArrayParam), stringId));
}

TEUCHOS_UNIT_TEST(TreeItem_Functions, ArrayTypes){
  RCP<ParameterEntry> twoDArrayParameter = 
    rcp(new ParameterEntry(Teuchos::TwoDArray<int>(2,2,4)));
  QString type = TreeItem::getTypeId(twoDArrayParameter);
  TEST_EQUALITY(
    type.toStdString(), 
    twoDArrayId.toStdString() + " " + intId.toStdString());
}


} //end namespace
