#ifndef CEPETRA_BLOCKMAP_H
#define CEPETRA_BLOCKMAP_H

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


/*! @file CEpetra_BlockMap.h
 * @brief Wrappers for Epetra_BlockMap */

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
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_BlockMap_Generalize ( 
  CT_Epetra_BlockMap_ID_t id );

/*@}*/

/*! @name Epetra_BlockMap constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm)
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( 
  int NumGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/*! @brief Wrapper for 
   Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm)
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int ElementSize, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );

/*! @brief Wrapper for 
   Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm)
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID );

/*! @brief Wrapper for 
   Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, const int *ElementSizeList, int IndexBase, const Epetra_Comm& Comm)
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, const int * ElementSizeList, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID );

/*! @brief Wrapper for 
   Epetra_BlockMap::Epetra_BlockMap(const Epetra_BlockMap& map)
*/
CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( 
  CT_Epetra_BlockMap_ID_t mapID );

/*@}*/

/*! @name Epetra_BlockMap destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_BlockMap::~Epetra_BlockMap(void)
*/
void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID );

/*@}*/

/*! @name Epetra_BlockMap member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Epetra_BlockMap::RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList) const
*/
int Epetra_BlockMap_RemoteIDList ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList );

/*! @brief Wrapper for 
   int Epetra_BlockMap::RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList, int * SizeList) const
*/
int Epetra_BlockMap_RemoteIDList_WithSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList, int * SizeList );

/*! @brief Wrapper for 
   int Epetra_BlockMap::LID(int GID) const
*/
int Epetra_BlockMap_LID ( CT_Epetra_BlockMap_ID_t selfID, int GID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::GID(int LID) const
*/
int Epetra_BlockMap_GID ( CT_Epetra_BlockMap_ID_t selfID, int LID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::FindLocalElementID(int PointID, int & ElementID, int & ElementOffset) const
*/
int Epetra_BlockMap_FindLocalElementID ( 
  CT_Epetra_BlockMap_ID_t selfID, int PointID, int * ElementID, 
  int * ElementOffset );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::MyGID(int GID_in) const
*/
boolean Epetra_BlockMap_MyGID ( 
  CT_Epetra_BlockMap_ID_t selfID, int GID_in );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::MyLID(int LID_in) const
*/
boolean Epetra_BlockMap_MyLID ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID_in );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MinAllGID() const
*/
int Epetra_BlockMap_MinAllGID ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MaxAllGID() const
*/
int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MinMyGID() const
*/
int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MaxMyGID() const
*/
int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MinLID() const
*/
int Epetra_BlockMap_MinLID ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MaxLID() const
*/
int Epetra_BlockMap_MaxLID ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::NumGlobalElements() const
*/
int Epetra_BlockMap_NumGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::NumMyElements() const
*/
int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MyGlobalElements(int * MyGlobalElementList) const
*/
int Epetra_BlockMap_MyGlobalElements_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * MyGlobalElementList );

/*! @brief Wrapper for 
   int Epetra_BlockMap::ElementSize() const
*/
int Epetra_BlockMap_ElementSize_Const ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::ElementSize(int LID) const
*/
int Epetra_BlockMap_ElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::FirstPointInElement(int LID) const
*/
int Epetra_BlockMap_FirstPointInElement ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::IndexBase() const
*/
int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::NumGlobalPoints() const
*/
int Epetra_BlockMap_NumGlobalPoints ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::NumMyPoints() const
*/
int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MinMyElementSize() const
*/
int Epetra_BlockMap_MinMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MaxMyElementSize() const
*/
int Epetra_BlockMap_MaxMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MinElementSize() const
*/
int Epetra_BlockMap_MinElementSize ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::MaxElementSize() const
*/
int Epetra_BlockMap_MaxElementSize ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::UniqueGIDs() const
*/
boolean Epetra_BlockMap_UniqueGIDs ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::ConstantElementSize() const
*/
boolean Epetra_BlockMap_ConstantElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::SameAs(const Epetra_BlockMap & Map) const
*/
boolean Epetra_BlockMap_SameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::PointSameAs(const Epetra_BlockMap & Map) const
*/
boolean Epetra_BlockMap_PointSameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::LinearMap() const
*/
boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::DistributedGlobal() const
*/
boolean Epetra_BlockMap_DistributedGlobal ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_BlockMap::MyGlobalElements() const
*/
int * Epetra_BlockMap_MyGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_BlockMap::FirstPointInElementList() const
*/
int * Epetra_BlockMap_FirstPointInElementList ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_BlockMap::ElementSizeList() const
*/
int * Epetra_BlockMap_ElementSizeList ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int * Epetra_BlockMap::PointToElementList() const
*/
int * Epetra_BlockMap_PointToElementList ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_BlockMap::ElementSizeList(int * ElementSizeList)const
*/
int Epetra_BlockMap_ElementSizeList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * ElementSizeList );

/*! @brief Wrapper for 
   int Epetra_BlockMap::FirstPointInElementList(int * FirstPointInElementList)const
*/
int Epetra_BlockMap_FirstPointInElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * FirstPointInElementList );

/*! @brief Wrapper for 
   int Epetra_BlockMap::PointToElementList(int * PointToElementList) const
*/
int Epetra_BlockMap_PointToElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * PointToElementList );

/*! @brief Wrapper for 
   const Epetra_Comm & Epetra_BlockMap::Comm() const
*/
CT_Epetra_Comm_ID_t Epetra_BlockMap_Comm ( 
  CT_Epetra_BlockMap_ID_t selfID );

/*! @brief Wrapper for 
   bool Epetra_BlockMap::IsOneToOne() const
*/
boolean Epetra_BlockMap_IsOneToOne ( CT_Epetra_BlockMap_ID_t selfID );

/*@}*/

/*! @name Epetra_BlockMap operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_BlockMap & Epetra_BlockMap::operator=(const Epetra_BlockMap & map)
*/
void Epetra_BlockMap_Assign ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_BLOCKMAP_H */

