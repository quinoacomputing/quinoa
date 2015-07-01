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
#include "Optika_ArrayHelperFunctions.hpp"

namespace Optika{

QString determineArrayType(RCP<const ParameterEntry> parameter, bool twoD){
	any anyArray = parameter->getAny();
	if(anyArray.type() == (twoD ? typeid(TwoDArray<int>) : typeid(Array<int>))){
		return intId;
	}
	else if(anyArray.type() == (twoD ? typeid(TwoDArray<short>) : typeid(Array<short>))){
		return shortId;
	}
	else if(anyArray.type() == (twoD ? typeid(TwoDArray<double>) : typeid(Array<double>))){
		return doubleId;
	}
	else if(anyArray.type() == (twoD ? typeid(TwoDArray<float>) : typeid(Array<float>))){
		return floatId;
	}
	else if(anyArray.type() == (twoD ? typeid(TwoDArray<std::string>) : typeid(Array<std::string>))){
		return stringId;
	}
	else{
		return unrecognizedId;		
	}
}

QString determineArrayType(RCP<const ParameterEntry> parameter){
	any anyArray = parameter->getAny();
	if(anyArray.type() == typeid(Array<int>)){
		return intId;
	}
	else if(anyArray.type() == typeid(Array<short>)){
		return shortId;
	}
	else if(anyArray.type() == typeid(Array<double>)){
		return doubleId;
	}
	else if(anyArray.type() == typeid(Array<float>)){
		return floatId;
	}
	else if(anyArray.type() == typeid(Array<std::string>)){
		return stringId;
	}
	else{
		return unrecognizedId;		
	}
}

QVariant arrayEntryToVariant(
  RCP<const ParameterEntry> arrayEntry, QString type, bool twoD){
	if(type == intId){
    return (twoD ? 
      QVariant::fromValue<TwoDArray<int> >(
      getValue<TwoDArray<int> >(*arrayEntry))
      :
      QVariant::fromValue<Array<int> >(
      getValue<Array<int> >(*arrayEntry)));
	}
	else if(type == shortId){
    return (twoD ? 
      QVariant::fromValue<TwoDArray<short> >(
      getValue<TwoDArray<short> >(*arrayEntry))
      :
      QVariant::fromValue<Array<short> >(
      getValue<Array<short> >(*arrayEntry)));
	}
	else if(type == doubleId){
    return (twoD ? 
      QVariant::fromValue<TwoDArray<double> >(
      getValue<TwoDArray<double> >(*arrayEntry))
      :
      QVariant::fromValue<Array<double> >(
      getValue<Array<double> >(*arrayEntry)));
  }
	else if(type == floatId){
    return (twoD ? 
      QVariant::fromValue<TwoDArray<float> >(
      getValue<TwoDArray<float> >(*arrayEntry))
      :
      QVariant::fromValue<Array<float> >(
      getValue<Array<float> >(*arrayEntry)));
  }
	else if(type == stringId){
    return (twoD ? 
      QVariant::fromValue<TwoDArray<std::string> >(
      getValue<TwoDArray<std::string> >(*arrayEntry))
      :
      QVariant::fromValue<Array<std::string> >(
      getValue<Array<std::string> >(*arrayEntry)));
	}
  return QVariant();
}

QString getArrayType(QString itemType){
  return itemType.section(" ",-1);  
}

bool isArrayEmpty(RCP<const ParameterEntry> arrayEntry, QString type){
	if(type == intId){
    return getValue<Array<int> >(*arrayEntry).size() == 0;
	}
	else if(type == shortId){
    return getValue<Array<short> >(*arrayEntry).size() == 0;
	}
	else if(type == doubleId){
    return getValue<Array<double> >(*arrayEntry).size() == 0;
  }
	else if(type == floatId){
    return getValue<Array<float> >(*arrayEntry).size() == 0;
  }
	else if(type == stringId){
    return getValue<Array<std::string> >(*arrayEntry).size() == 0;
	}
  return true;
}


}
