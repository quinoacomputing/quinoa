#ifndef CTEUCHOS_ANY_H
#define CTEUCHOS_ANY_H

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


/*! @file CTeuchos_any.h
 * @brief Wrappers for Teuchos::any */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name any constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::any::any()
*/
CT_Teuchos_any_ID_t Teuchos_any_Create (  );

/*! @brief Wrapper for 
   template<typename ValueType> explicit Teuchos::any::any(const ValueType & value)
*/
CT_Teuchos_any_ID_t Teuchos_any_Create_double ( double value );

/*! @brief Wrapper for 
   template<typename ValueType> explicit Teuchos::any::any(const ValueType & value)
*/
CT_Teuchos_any_ID_t Teuchos_any_Create_int ( int value );

/*! @brief Wrapper for 
   Teuchos::any::any(const any & other)
*/
CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( 
  CT_Teuchos_any_ID_t otherID );

/*@}*/

/*! @name any destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::any::~any()
*/
void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID );

/*@}*/

/*! @name any member wrappers */
/*@{*/

/*! @brief Wrapper for 
   any & Teuchos::any::swap(any & rhs)
*/
CT_Teuchos_any_ID_t Teuchos_any_swap ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

/*! @brief Wrapper for 
   bool Teuchos::any::empty() const
*/
boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID );

/*! @brief Wrapper for 
   std::string Teuchos::any::typeName() const
*/
const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID );

/*! @brief Wrapper for 
   bool Teuchos::any::same( const any &other ) const
*/
boolean Teuchos_any_same ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID );

/*@}*/

/*! @name any operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   any & Teuchos::any::operator=(const any & rhs)
*/
void Teuchos_any_Assign ( 
  CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_ANY_H */

