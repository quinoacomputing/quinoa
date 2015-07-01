#ifndef CAMESOS_H
#define CAMESOS_H

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


#ifdef HAVE_CTRILINOS_AMESOS



/*! @file CAmesos.h
 * @brief Wrappers for Amesos */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name Amesos constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Amesos::Amesos()
*/
CT_Amesos_ID_t Amesos_Create (  );

/*@}*/

/*! @name Amesos destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Amesos::~Amesos()
*/
void Amesos_Destroy ( CT_Amesos_ID_t * selfID );

/*@}*/

/*! @name Amesos member wrappers */
/*@{*/

/*! @brief Wrapper for 
   Amesos_BaseSolver *Amesos::Create(const char *ClassType, const Epetra_LinearProblem& LinearProblem )
*/
CT_Amesos_BaseSolver_ID_t Amesos_CreateSolver ( 
  CT_Amesos_ID_t selfID, const char * ClassType, 
  CT_Epetra_LinearProblem_ID_t LinearProblemID );

/*! @brief Wrapper for 
   bool Amesos::Query(const char * ClassType)
*/
boolean Amesos_Query ( 
  CT_Amesos_ID_t selfID, const char * ClassType );

/*@}*/

/*! @name Amesos static function wrappers */
/*@{*/

/*! @brief Wrapper for 
   static Teuchos::ParameterList Amesos::GetValidParameters()
*/
CT_Teuchos_ParameterList_ID_t Amesos_GetValidParameters (  );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_CTRILINOS_AMESOS */
#endif /* CAMESOS_H */

