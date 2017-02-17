// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#ifndef TEUCHOS_FUNTIONOBJECTXMLCONVERTERDB_HPP
#define TEUCHOS_FUNTIONOBJECTXMLCONVERTERDB_HPP

/*! \file Teuchos_FunctionObjectXMLConverterDB.hpp
 * \brief A database for FunctionObjectXMLConverters.
*/

// All includes needed for convience macros below
#include "Teuchos_StandardFunctionObjectXMLConverters.hpp"
#include "Teuchos_StandardFunctionObjects.hpp"
#include "Teuchos_DummyObjectGetter.hpp"


namespace Teuchos {

/**
 * \brief Provides ability to lookup FunctionObjectXMLConverters
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT FunctionObjectXMLConverterDB {

public:

  /** \name Modifier Functions */
  //@{

  /** \brief Add a converter to the database.
   *
   * \param function A dummy FunctionObject representing the type of function
   * this converter will convert.
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(
    RCP<const FunctionObject> function,
    RCP<FunctionObjectXMLConverter> converterToAdd);

  //@}

  /** \name Converter Functions */
  //@{

  /** \brief Get an appropriate FunctionObjectXMLConverter given a
   *  FunctionObject.
   *
   * \param function The FunctionObject for which a converter is
   * desired.
   *
   * \return A converter for the function.
   */
  static RCP<const FunctionObjectXMLConverter>
    getConverter(const FunctionObject& function);

  /** \brief Get an appropriate FunctionObjectXMLConverter given a XMLObject.
   *
   * @param xmlObject The XMLObject for which a converter is desired.
   *
   * @return A converter for the XMLObject.
   */
  static RCP<const FunctionObjectXMLConverter>
    getConverter(const XMLObject& xmlObject);

  /**
   * \brief Given a FunctionObject, converts the FunctionObject to XML.
   *
   * \param function The FunctionObject to convert.
   *
   * \return XML representation of the function.
   */
  static XMLObject convertFunctionObject(RCP<const FunctionObject> function);

  /**
   * \brief Given an XMLObject, converts the XMLObject
   * to a FunctionObject.
   *
   * \param xmlObject The XMLObject to be converted.
   *
   * \return A FunctionObject that was represented by the XML.
   */
  static RCP<FunctionObject> convertXML(const XMLObject& xmlObject);
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief prints the xml tags associated with all known converters
   *
   * \param out Stream to which tags should be printed.
   */
  static void printKnownConverters(std::ostream& out){
    out << "Known ConditionXMLConverters: " << std::endl;
    for(
      ConverterMap::const_iterator it = getConverterMap().begin();
      it != getConverterMap().end();
      ++it)
    {
      out << "\t" << it->first <<std::endl;
    }
  }

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief convience class. */
  typedef std::map<std::string, RCP<FunctionObjectXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<FunctionObjectXMLConverter> > ConverterPair;

  /** \brief Gets the default converter to be used to convert
   * FunctionObjects.
   */
  static ConverterMap& getConverterMap();

  //@}

};


} // end namespace Teuchos

//
// Helper Macros
//

/** \brief Adds a SubtractionFunction, AdditionFunction, MultiplicationFunction,
 * and DivisionFunction of type T to the Database
 */
#define TEUCHOS_ADD_SIMPLEFUNCTIONCONVERTERS(T) \
  Teuchos::FunctionObjectXMLConverterDB::addConverter( \
    Teuchos::rcp(new Teuchos::SubtractionFunction< T >), \
    Teuchos::DummyObjectGetter<Teuchos::SubtractionFunctionXMLConverter< T > >:: \
      getDummyObject()); \
      \
  Teuchos::FunctionObjectXMLConverterDB::addConverter( \
    Teuchos::rcp(new Teuchos::AdditionFunction< T >), \
    Teuchos::DummyObjectGetter<Teuchos::AdditionFunctionXMLConverter< T > >:: \
      getDummyObject()); \
      \
  Teuchos::FunctionObjectXMLConverterDB::addConverter( \
    Teuchos::rcp(new Teuchos::MultiplicationFunction< T >), \
    Teuchos::DummyObjectGetter<Teuchos::MultiplicationFunctionXMLConverter< T > >:: \
      getDummyObject()); \
      \
  Teuchos::FunctionObjectXMLConverterDB::addConverter( \
    Teuchos::rcp(new Teuchos::DivisionFunction< T >), \
    Teuchos::DummyObjectGetter<Teuchos::DivisionFunctionXMLConverter< T > >:: \
      getDummyObject()); \
      \

#endif // TEUCHOS_FUNTIONOBJECTXMLCONVERTERDB_HPP
