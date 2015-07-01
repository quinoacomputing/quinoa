#ifndef CTEUCHOS_PARAMETERENTRY_H
#define CTEUCHOS_PARAMETERENTRY_H

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


/*! @file CTeuchos_ParameterEntry.h
 * @brief Wrappers for Teuchos::ParameterEntry */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ParameterEntry constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::ParameterEntry::ParameterEntry()
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Create (  );

/*! @brief Wrapper for 
   Teuchos::ParameterEntry::ParameterEntry(const ParameterEntry& source)
*/
CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Duplicate ( 
  CT_Teuchos_ParameterEntry_ID_t sourceID );

/*@}*/

/*! @name ParameterEntry destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::ParameterEntry::~ParameterEntry()
*/
void Teuchos_ParameterEntry_Destroy ( 
  CT_Teuchos_ParameterEntry_ID_t * selfID );

/*@}*/

/*! @name ParameterEntry member wrappers */
/*@{*/

/*! @brief Wrapper for 
   void Teuchos::ParameterEntry::setAnyValue( const any &value, bool isDefault = false )
*/
void Teuchos_ParameterEntry_setAnyValue ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_any_ID_t valueID, boolean isDefault );

/*! @brief Wrapper for 
   void Teuchos::ParameterEntry::setDocString(const std::string &docString)
*/
void Teuchos_ParameterEntry_setDocString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, const char docString[] );

/*! @brief Wrapper for 
   ParameterList& Teuchos::ParameterEntry::setList( bool isDefault = false, const std::string &docString = "" )
*/
CT_Teuchos_ParameterList_ID_t Teuchos_ParameterEntry_setList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean isDefault, 
  const char docString[] );

/*! @brief Wrapper for 
   template<typename T> inline T& Teuchos::ParameterEntry::getValue(T *ptr) const
*/
double Teuchos_ParameterEntry_getValue_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, double * ptr );

/*! @brief Wrapper for 
   template<typename T> inline T& Teuchos::ParameterEntry::getValue(T *ptr) const
*/
int Teuchos_ParameterEntry_getValue_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, int * ptr );

/*! @brief Wrapper for 
   inline any& Teuchos::ParameterEntry::getAny(bool activeQry = true)
*/
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

/*! @brief Wrapper for 
   inline const any& Teuchos::ParameterEntry::getAny(bool activeQry = true) const
*/
CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny_const ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

/*! @brief Wrapper for 
   inline bool Teuchos::ParameterEntry::isUsed() const
*/
boolean Teuchos_ParameterEntry_isUsed ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   bool Teuchos::ParameterEntry::isList() const
*/
boolean Teuchos_ParameterEntry_isList ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   template <typename T> inline bool Teuchos::ParameterEntry::isType() const
*/
boolean Teuchos_ParameterEntry_isType_double ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   template <typename T> inline bool Teuchos::ParameterEntry::isType() const
*/
boolean Teuchos_ParameterEntry_isType_int ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   inline bool Teuchos::ParameterEntry::isDefault() const
*/
boolean Teuchos_ParameterEntry_isDefault ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*! @brief Wrapper for 
   inline std::string Teuchos::ParameterEntry::docString() const
*/
const char * Teuchos_ParameterEntry_docString ( 
  CT_Teuchos_ParameterEntry_ID_t selfID );

/*@}*/

/*! @name ParameterEntry operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   ParameterEntry& Teuchos::ParameterEntry::operator=(const ParameterEntry& source)
*/
void Teuchos_ParameterEntry_Assign ( 
  CT_Teuchos_ParameterEntry_ID_t selfID, 
  CT_Teuchos_ParameterEntry_ID_t sourceID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_PARAMETERENTRY_H */

