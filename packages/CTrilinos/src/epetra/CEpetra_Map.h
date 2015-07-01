#ifndef CEPETRA_MAP_H
#define CEPETRA_MAP_H

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


/*! @file CEpetra_Map.h
 * @brief Wrappers for Epetra_Map */

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
CT_Epetra_Map_ID_t Epetra_Map_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Map_Generalize ( 
  CT_Epetra_Map_ID_t id );

/*@}*/

/*! @name Epetra_Map constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Map::Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm)
*/
CT_Epetra_Map_ID_t Epetra_Map_Create ( 
  int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

/*! @brief Wrapper for 
   Epetra_Map::Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm)
*/
CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/*! @brief Wrapper for 
   Epetra_Map::Epetra_Map(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int IndexBase, const Epetra_Comm& Comm)
*/
CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/*! @brief Wrapper for 
   Epetra_Map::Epetra_Map(const Epetra_Map& map)
*/
CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

/*@}*/

/*! @name Epetra_Map destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Map::~Epetra_Map(void)
*/
void Epetra_Map_Destroy ( CT_Epetra_Map_ID_t * selfID );

/*@}*/

/*! @name Epetra_Map operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Map & Epetra_Map::operator=(const Epetra_Map & map)
*/
void Epetra_Map_Assign ( 
  CT_Epetra_Map_ID_t selfID, CT_Epetra_Map_ID_t mapID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MAP_H */

