/*
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
*/

#ifndef OPTIKA_CONFIG_DEFS_HPP
#define OPTIKA_CONFIG_DEFS_HPP

/*! \file Optika_ConfigDefs.hpp
 * \brief A Header file that includes some of the commonly used includes
 * throughtout Optika as well as some overall information used throughout
 * Optika.
 */
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_TwoDArray.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"


namespace Optika{

  //Common Teuchos classes that are used.
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::null;
  using Teuchos::is_null;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_static_cast;
  using Teuchos::ParameterEntryValidator;
  using Teuchos::EnhancedNumberValidator;
  using Teuchos::EnhancedNumberTraits;
  using Teuchos::FileNameValidator;
  using Teuchos::ArrayValidator;
  using Teuchos::DependencySheet;
  using Teuchos::Dependency;
  using Teuchos::VisualDependency;
  using Teuchos::any;
  using Teuchos::any_cast;
  using Teuchos::XMLParameterListWriter;
  using Teuchos::XMLObject;
  using Teuchos::getValue;
  using Teuchos::getParametersFromXmlFile;
  using Teuchos::TwoDArray;
  using Teuchos::TwoDArrayValidator;
  using Teuchos::rcp;


  
} //namespace Optika

//Declarations of the supported Array types
//as QMEATTYPES so that we can encapsulate them
//in QVariant objects.
Q_DECLARE_METATYPE(Teuchos::Array<int>)
Q_DECLARE_METATYPE(Teuchos::Array<short>)
Q_DECLARE_METATYPE(Teuchos::Array<float>)
Q_DECLARE_METATYPE(Teuchos::Array<double>)
Q_DECLARE_METATYPE(Teuchos::Array<std::string>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<int>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<short>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<float>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<double>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<std::string>)


#endif //OPTIKA_CONFIG_DEFS_HPP
