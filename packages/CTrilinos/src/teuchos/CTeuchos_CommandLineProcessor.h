#ifndef CTEUCHOS_COMMANDLINEPROCESSOR_H
#define CTEUCHOS_COMMANDLINEPROCESSOR_H

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


/*! @file CTeuchos_CommandLineProcessor.h
 * @brief Wrappers for Teuchos::CommandLineProcessor */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name CommandLineProcessor constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::CommandLineProcessor::CommandLineProcessor( bool throwExceptions = true ,bool recogniseAllOptions = true ,bool addOutputSetupOptions = false )
*/
CT_Teuchos_CommandLineProcessor_ID_t Teuchos_CommandLineProcessor_Create ( 
  boolean throwExceptions, boolean recogniseAllOptions, 
  boolean addOutputSetupOptions );

/*@}*/

/*! @name CommandLineProcessor destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Teuchos::CommandLineProcessor::~CommandLineProcessor()
*/
void Teuchos_CommandLineProcessor_Destroy ( 
  CT_Teuchos_CommandLineProcessor_ID_t * selfID );

/*@}*/

/*! @name CommandLineProcessor member wrappers */
/*@{*/

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::throwExceptions( const bool & throwExceptions )
*/
void Teuchos_CommandLineProcessor_throwExceptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const boolean throwExceptions );

/*! @brief Wrapper for 
   bool Teuchos::CommandLineProcessor::throwExceptions() const
*/
boolean Teuchos_CommandLineProcessor_throwExceptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID );

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::recogniseAllOptions( const bool & recogniseAllOptions )
*/
void Teuchos_CommandLineProcessor_recogniseAllOptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const boolean recogniseAllOptions );

/*! @brief Wrapper for 
   bool Teuchos::CommandLineProcessor::recogniseAllOptions() const
*/
boolean Teuchos_CommandLineProcessor_recogniseAllOptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID );

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::addOutputSetupOptions( const bool &addOutputSetupOptions )
*/
void Teuchos_CommandLineProcessor_addOutputSetupOptions_set ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const boolean addOutputSetupOptions );

/*! @brief Wrapper for 
   bool Teuchos::CommandLineProcessor::addOutputSetupOptions() const
*/
boolean Teuchos_CommandLineProcessor_addOutputSetupOptions_get ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID );

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::setDocString( const char doc_string[] )
*/
void Teuchos_CommandLineProcessor_setDocString ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char doc_string[] );

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::setOption( const char option_true[] ,const char option_false[] ,bool *option_val ,const char documentation[] = NULL )
*/
void Teuchos_CommandLineProcessor_setOption_bool ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_true[], const char option_false[], 
  boolean * option_val, const char documentation[] );

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::setOption( const char option_name[] ,int *option_val ,const char documentation[] = NULL ,const bool required = false )
*/
void Teuchos_CommandLineProcessor_setOption_int ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], int * option_val, 
  const char documentation[], const boolean required );

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::setOption( const char option_name[] ,double *option_val ,const char documentation[] = NULL ,const bool required = false )
*/
void Teuchos_CommandLineProcessor_setOption_double ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], double * option_val, 
  const char documentation[], const boolean required );

/*! @brief Wrapper for 
   void Teuchos::CommandLineProcessor::setOption( const char option_name[] ,std::string *option_val ,const char documentation[] = NULL ,const bool required = false )
*/
void Teuchos_CommandLineProcessor_setOption_str ( 
  CT_Teuchos_CommandLineProcessor_ID_t selfID, 
  const char option_name[], char * option_val[], 
  const char documentation[], const boolean required );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CTEUCHOS_COMMANDLINEPROCESSOR_H */

