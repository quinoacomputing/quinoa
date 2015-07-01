#ifndef CEPETRA_DIRECTORY_H
#define CEPETRA_DIRECTORY_H

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


/*! @file CEpetra_Directory.h
 * @brief Wrappers for Epetra_Directory */

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
CT_Epetra_Directory_ID_t Epetra_Directory_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Directory_Generalize ( 
  CT_Epetra_Directory_ID_t id );

/*@}*/

/*! @name Epetra_Directory destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Directory::~Epetra_Directory()
*/
void Epetra_Directory_Destroy ( CT_Epetra_Directory_ID_t * selfID );

/*@}*/

/*! @name Epetra_Directory member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual int Epetra_Directory::GetDirectoryEntries( const Epetra_BlockMap& Map, const int NumEntries, const int * GlobalEntries, int * Procs, int * LocalEntries, int * EntrySizes, bool high_rank_sharing_procs=false) const = 0
*/
int Epetra_Directory_GetDirectoryEntries ( 
  CT_Epetra_Directory_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID, 
  const int NumEntries, const int * GlobalEntries, int * Procs, 
  int * LocalEntries, int * EntrySizes, 
  boolean high_rank_sharing_procs );

/*! @brief Wrapper for 
   virtual bool Epetra_Directory::GIDsAllUniquelyOwned() const = 0
*/
boolean Epetra_Directory_GIDsAllUniquelyOwned ( 
  CT_Epetra_Directory_ID_t selfID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DIRECTORY_H */

