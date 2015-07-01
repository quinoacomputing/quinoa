#ifndef CTEUCHOS_PARAMETERLIST_H
#define CTEUCHOS_PARAMETERLIST_H

/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"


/*! @file CTeuchos_ParameterList.h
 * @brief Wrappers for Teuchos::ParameterList */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Teuchos_ParameterList_Generalize ( 
  CT_Teuchos_ParameterList_ID_t id );

/*@}*/

/*! @name ParameterList constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::ParameterList::ParameterList()
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  );

/*! @brief Wrapper for 
   Teuchos::ParameterList::ParameterList(const std::string &name)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( 
  const char name[] );

/*! @brief Wrapper for 
   Teuchos::ParameterList::ParameterList(const ParameterList& source)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*@}*/

/*! @name ParameterList destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Teuchos::ParameterList::~ParameterList()
*/
void Teuchos_ParameterList_Destroy ( 
  CT_Teuchos_ParameterList_ID_t * selfID );

/*@}*/

/*! @name ParameterList member wrappers */
/*@{*/

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setName( const std::string &name )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setParameters(const ParameterList& source)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setParametersNotAlreadySet(const ParameterList& source)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::disableRecursiveValidation()
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/*! @brief Wrapper for 
   template<typename T> ParameterList& Teuchos::ParameterList::set( std::string const& name, T const& value, std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = Teuchos::null )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  double value, char const docString[] );

/*! @brief Wrapper for 
   template<typename T> ParameterList& Teuchos::ParameterList::set( std::string const& name, T const& value, std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = Teuchos::null )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  int value, char const docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::set( std::string const& name, char value[], std::string const& docString = "" ,RCP<const ParameterEntryValidator> const& validator = Teuchos::null )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set_str ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  char value[], char const docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::set( std::string const& name, ParameterList const& value, std::string const& docString = "" )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  CT_Teuchos_ParameterList_ID_t valueID, char const docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::setEntry(const std::string& name, const ParameterEntry& entry)
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  CT_Teuchos_ParameterEntry_ID_t entryID );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name, T def_value)
*/
double Teuchos_ParameterList_get_def_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double def_value );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name, T def_value)
*/
int Teuchos_ParameterList_get_def_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int def_value );

/*! @brief Wrapper for 
   std::string& Teuchos::ParameterList::get(const std::string& name, char def_value[])
*/
const char * Teuchos_ParameterList_get_def_char ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  char def_value[] );

/*! @brief Wrapper for 
   std::string& Teuchos::ParameterList::get(const std::string& name, const char def_value[])
*/
const char * Teuchos_ParameterList_get_def_const_char ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  const char def_value[] );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name)
*/
double Teuchos_ParameterList_get_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> T& Teuchos::ParameterList::get(const std::string& name)
*/
int Teuchos_ParameterList_get_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> const T& Teuchos::ParameterList::get(const std::string& name) const
*/
double Teuchos_ParameterList_get_const_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> const T& Teuchos::ParameterList::get(const std::string& name) const
*/
int Teuchos_ParameterList_get_const_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline T* Teuchos::ParameterList::getPtr(const std::string& name)
*/
double * Teuchos_ParameterList_getPtr_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline T* Teuchos::ParameterList::getPtr(const std::string& name)
*/
int * Teuchos_ParameterList_getPtr_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline const T* Teuchos::ParameterList::getPtr(const std::string& name) const
*/
const double * Teuchos_ParameterList_getPtr_const_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> inline const T* Teuchos::ParameterList::getPtr(const std::string& name) const
*/
const int * Teuchos_ParameterList_getPtr_const_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   ParameterEntry& Teuchos::ParameterList::getEntry(const std::string& name)
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   inline const ParameterEntry& Teuchos::ParameterList::getEntry(const std::string& name) const
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   inline ParameterEntry* Teuchos::ParameterList::getEntryPtr(const std::string& name)
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   inline const ParameterEntry* Teuchos::ParameterList::getEntryPtr(const std::string& name) const
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   bool Teuchos::ParameterList::remove( std::string const& name, bool throwIfNotExists = true )
*/
boolean Teuchos_ParameterList_remove ( 
  CT_Teuchos_ParameterList_ID_t selfID, char const name[], 
  boolean throwIfNotExists );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::sublist( const std::string& name, bool mustAlreadyExist = false ,const std::string& docString = "" )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  boolean mustAlreadyExist, const char docString[] );

/*! @brief Wrapper for 
   const ParameterList& Teuchos::ParameterList::sublist(const std::string& name) const
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   const std::string& Teuchos::ParameterList::name() const
*/
const char * Teuchos_ParameterList_name_it ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/*! @brief Wrapper for 
   bool Teuchos::ParameterList::isParameter(const std::string& name) const
*/
boolean Teuchos_ParameterList_isParameter ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   bool Teuchos::ParameterList::isSublist(const std::string& name) const
*/
boolean Teuchos_ParameterList_isSublist ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name) const
*/
boolean Teuchos_ParameterList_isType_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name) const
*/
boolean Teuchos_ParameterList_isType_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name, T* ptr) const
*/
boolean Teuchos_ParameterList_isType_type_double ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  double * ptr );

/*! @brief Wrapper for 
   template<typename T> bool Teuchos::ParameterList::isType(const std::string& name, T* ptr) const
*/
boolean Teuchos_ParameterList_isType_type_int ( 
  CT_Teuchos_ParameterList_ID_t selfID, const char name[], 
  int * ptr );

/*! @brief Wrapper for 
   std::string Teuchos::ParameterList::currentParametersString() const
*/
const char * Teuchos_ParameterList_currentParametersString ( 
  CT_Teuchos_ParameterList_ID_t selfID );

/*! @brief Wrapper for 
   void Teuchos::ParameterList::validateParameters( ParameterList const& validParamList, int const depth = 1000, EValidateUsed const validateUsed = VALIDATE_USED_ENABLED, EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED ) const
*/
void Teuchos_ParameterList_validateParameters ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, 
  const CT_EValidateUsed_E_t validateUsed, 
  const CT_EValidateDefaults_E_t validateDefaults );

/*! @brief Wrapper for 
   void Teuchos::ParameterList::validateParametersAndSetDefaults( ParameterList const& validParamList, int const depth = 1000 )
*/
void Teuchos_ParameterList_validateParametersAndSetDefaults ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t validParamListID, int const depth );

/*@}*/

/*! @name ParameterList operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterList::operator=(const ParameterList& source)
*/
void Teuchos_ParameterList_Assign ( 
  CT_Teuchos_ParameterList_ID_t selfID, 
  CT_Teuchos_ParameterList_ID_t sourceID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_PARAMETERLIST_H */

