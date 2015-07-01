
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

#include "CTrilinos_enums.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "Epetra_BlockMap.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_Comm_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_BlockMap */
Table<Epetra_BlockMap>& tableOfBlockMaps()
{
    static Table<Epetra_BlockMap> loc_tableOfBlockMaps(CT_Epetra_BlockMap_ID);
    return loc_tableOfBlockMaps;
}


} // namespace


//
// Definitions from CEpetra_BlockMap.h
//


extern "C" {


CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_BlockMap_Generalize ( 
  CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id);
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( 
  int NumGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewBlockMap(new Epetra_BlockMap(NumGlobalElements, 
        ElementSize, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( 
  int NumGlobalElements, int NumMyElements, int ElementSize, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewBlockMap(new Epetra_BlockMap(NumGlobalElements, 
        NumMyElements, ElementSize, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, int ElementSize, int IndexBase, 
  CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewBlockMap(new Epetra_BlockMap(NumGlobalElements, 
        NumMyElements, MyGlobalElements, ElementSize, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( 
  int NumGlobalElements, int NumMyElements, 
  const int * MyGlobalElements, const int * ElementSizeList, 
  int IndexBase, CT_Epetra_Comm_ID_t CommID )
{
    const Teuchos::RCP<const Epetra_Comm> Comm = CEpetra::getConstComm(CommID);
    return CEpetra::storeNewBlockMap(new Epetra_BlockMap(NumGlobalElements, 
        NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, *Comm));
}

CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( 
  CT_Epetra_BlockMap_ID_t mapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> map = CEpetra::getConstBlockMap(
        mapID);
    return CEpetra::storeNewBlockMap(new Epetra_BlockMap(*map));
}

void Epetra_BlockMap_Destroy ( CT_Epetra_BlockMap_ID_t * selfID )
{
    CEpetra::removeBlockMap(selfID);
}

int Epetra_BlockMap_RemoteIDList ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList )
{
    return CEpetra::getConstBlockMap(selfID)->RemoteIDList(NumIDs, GIDList, 
        PIDList, LIDList);
}

int Epetra_BlockMap_RemoteIDList_WithSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int NumIDs, const int * GIDList, 
  int * PIDList, int * LIDList, int * SizeList )
{
    return CEpetra::getConstBlockMap(selfID)->RemoteIDList(NumIDs, GIDList, 
        PIDList, LIDList, SizeList);
}

int Epetra_BlockMap_LID ( CT_Epetra_BlockMap_ID_t selfID, int GID )
{
    return CEpetra::getConstBlockMap(selfID)->LID(GID);
}

int Epetra_BlockMap_GID ( CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CEpetra::getConstBlockMap(selfID)->GID(LID);
}

int Epetra_BlockMap_FindLocalElementID ( 
  CT_Epetra_BlockMap_ID_t selfID, int PointID, int * ElementID, 
  int * ElementOffset )
{
    return CEpetra::getConstBlockMap(selfID)->FindLocalElementID(PointID, 
        *ElementID, *ElementOffset);
}

boolean Epetra_BlockMap_MyGID ( 
  CT_Epetra_BlockMap_ID_t selfID, int GID_in )
{
    return ((CEpetra::getConstBlockMap(selfID)->MyGID(GID_in)) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_MyLID ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID_in )
{
    return ((CEpetra::getConstBlockMap(selfID)->MyLID(LID_in)) ? TRUE : FALSE);
}

int Epetra_BlockMap_MinAllGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinAllGID();
}

int Epetra_BlockMap_MaxAllGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxAllGID();
}

int Epetra_BlockMap_MinMyGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinMyGID();
}

int Epetra_BlockMap_MaxMyGID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxMyGID();
}

int Epetra_BlockMap_MinLID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinLID();
}

int Epetra_BlockMap_MaxLID ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxLID();
}

int Epetra_BlockMap_NumGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumGlobalElements();
}

int Epetra_BlockMap_NumMyElements ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumMyElements();
}

int Epetra_BlockMap_MyGlobalElements_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * MyGlobalElementList )
{
    return CEpetra::getConstBlockMap(selfID)->MyGlobalElements(
        MyGlobalElementList);
}

int Epetra_BlockMap_ElementSize_Const ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSize();
}

int Epetra_BlockMap_ElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSize(LID);
}

int Epetra_BlockMap_FirstPointInElement ( 
  CT_Epetra_BlockMap_ID_t selfID, int LID )
{
    return CEpetra::getConstBlockMap(selfID)->FirstPointInElement(LID);
}

int Epetra_BlockMap_IndexBase ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->IndexBase();
}

int Epetra_BlockMap_NumGlobalPoints ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumGlobalPoints();
}

int Epetra_BlockMap_NumMyPoints ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->NumMyPoints();
}

int Epetra_BlockMap_MinMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinMyElementSize();
}

int Epetra_BlockMap_MaxMyElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxMyElementSize();
}

int Epetra_BlockMap_MinElementSize ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MinElementSize();
}

int Epetra_BlockMap_MaxElementSize ( CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MaxElementSize();
}

boolean Epetra_BlockMap_UniqueGIDs ( CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CEpetra::getConstBlockMap(selfID)->UniqueGIDs()) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_ConstantElementSize ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CEpetra::getConstBlockMap(
        selfID)->ConstantElementSize()) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_SameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return ((CEpetra::getConstBlockMap(selfID)->SameAs(*Map)) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_PointSameAs ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return ((CEpetra::getConstBlockMap(selfID)->PointSameAs(
        *Map)) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_LinearMap ( CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CEpetra::getConstBlockMap(selfID)->LinearMap()) ? TRUE : FALSE);
}

boolean Epetra_BlockMap_DistributedGlobal ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CEpetra::getConstBlockMap(
        selfID)->DistributedGlobal()) ? TRUE : FALSE);
}

int * Epetra_BlockMap_MyGlobalElements ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->MyGlobalElements();
}

int * Epetra_BlockMap_FirstPointInElementList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->FirstPointInElementList();
}

int * Epetra_BlockMap_ElementSizeList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSizeList();
}

int * Epetra_BlockMap_PointToElementList ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::getConstBlockMap(selfID)->PointToElementList();
}

int Epetra_BlockMap_ElementSizeList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * ElementSizeList )
{
    return CEpetra::getConstBlockMap(selfID)->ElementSizeList(ElementSizeList);
}

int Epetra_BlockMap_FirstPointInElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * FirstPointInElementList )
{
    return CEpetra::getConstBlockMap(selfID)->FirstPointInElementList(
        FirstPointInElementList);
}

int Epetra_BlockMap_PointToElementList_Fill ( 
  CT_Epetra_BlockMap_ID_t selfID, int * PointToElementList )
{
    return CEpetra::getConstBlockMap(selfID)->PointToElementList(
        PointToElementList);
}

CT_Epetra_Comm_ID_t Epetra_BlockMap_Comm ( 
  CT_Epetra_BlockMap_ID_t selfID )
{
    return CEpetra::storeConstComm(&( CEpetra::getConstBlockMap(
        selfID)->Comm() ));
}

boolean Epetra_BlockMap_IsOneToOne ( CT_Epetra_BlockMap_ID_t selfID )
{
    return ((CEpetra::getConstBlockMap(selfID)->IsOneToOne()) ? TRUE : FALSE);
}

void Epetra_BlockMap_Assign ( 
  CT_Epetra_BlockMap_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID )
{
    Epetra_BlockMap& self = *( CEpetra::getBlockMap(selfID) );

    const Teuchos::RCP<const Epetra_BlockMap> map = CEpetra::getConstBlockMap(
        mapID);
    self = *map;
}


} // extern "C"


//
// Definitions from CEpetra_BlockMap_Cpp.hpp
//


/* get Epetra_BlockMap from non-const table using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    if (tableOfBlockMaps().isType(id.table))
        return tableOfBlockMaps().get<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
}

/* get Epetra_BlockMap from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_BlockMap>
CEpetra::getBlockMap( CTrilinos_Universal_ID_t id )
{
    if (tableOfBlockMaps().isType(id.table))
        return tableOfBlockMaps().get<Epetra_BlockMap>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_BlockMap>(id);
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    if (tableOfBlockMaps().isType(id.table))
        return tableOfBlockMaps().getConst<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_BlockMap>
CEpetra::getConstBlockMap( CTrilinos_Universal_ID_t id )
{
    if (tableOfBlockMaps().isType(id.table))
        return tableOfBlockMaps().getConst<Epetra_BlockMap>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_BlockMap>(id);
}

/* store Epetra_BlockMap (owned) in non-const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeNewBlockMap( Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store<Epetra_BlockMap>(pobj, true));
}

/* store Epetra_BlockMap in non-const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeBlockMap( Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store<Epetra_BlockMap>(pobj, false));
}

/* store const Epetra_BlockMap in const table */
CT_Epetra_BlockMap_ID_t
CEpetra::storeConstBlockMap( const Epetra_BlockMap *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tableOfBlockMaps().store<Epetra_BlockMap>(pobj, false));
}

/* remove Epetra_BlockMap from table using CT_Epetra_BlockMap_ID */
void
CEpetra::removeBlockMap( CT_Epetra_BlockMap_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(*id);
    if (tableOfBlockMaps().isType(aid.table))
        tableOfBlockMaps().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(aid);
}

/* remove Epetra_BlockMap from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeBlockMap( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfBlockMaps().isType(aid->table))
        tableOfBlockMaps().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_BlockMap table */
void
CEpetra::purgeBlockMap(  )
{
    tableOfBlockMaps().purge();
}

/* store Epetra_BlockMap in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasBlockMap( const Teuchos::RCP< Epetra_BlockMap > & robj )
{
    return tableOfBlockMaps().alias(robj);
}

/* store const Epetra_BlockMap in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstBlockMap( const Teuchos::RCP< const Epetra_BlockMap > & robj )
{
    return tableOfBlockMaps().alias(robj);
}



