
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
#include "CEpetra_MultiVector.h"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_MultiVector */
Table<Epetra_MultiVector>& tableOfMultiVectors()
{
    static Table<Epetra_MultiVector> loc_tableOfMultiVectors(CT_Epetra_MultiVector_ID);
    return loc_tableOfMultiVectors;
}


} // namespace


//
// Definitions from CEpetra_MultiVector.h
//


extern "C" {


CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Degeneralize ( 
  CTrilinos_Universal_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(id);
}

CTrilinos_Universal_ID_t Epetra_MultiVector_Generalize ( 
  CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id);
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeNewMultiVector(new Epetra_MultiVector(*Map, 
        NumVectors, ((zeroOut) != FALSE ? true : false)));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( 
  CT_Epetra_MultiVector_ID_t SourceID )
{
    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CEpetra::getConstMultiVector(SourceID);
    return CEpetra::storeNewMultiVector(new Epetra_MultiVector(*Source));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double * A, int MyLDA, int NumVectors )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeNewMultiVector(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Map, A, MyLDA, NumVectors));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_BlockMap_ID_t MapID, 
  double ** ArrayOfPointers, int NumVectors )
{
    const Teuchos::RCP<const Epetra_BlockMap> Map = CEpetra::getConstBlockMap(
        MapID);
    return CEpetra::storeNewMultiVector(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Map, ArrayOfPointers, NumVectors));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int * Indices, int NumVectors )
{
    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CEpetra::getConstMultiVector(SourceID);
    return CEpetra::storeNewMultiVector(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Source, Indices, NumVectors));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( 
  CT_Epetra_DataAccess_E_t CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int StartIndex, int NumVectors )
{
    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CEpetra::getConstMultiVector(SourceID);
    return CEpetra::storeNewMultiVector(new Epetra_MultiVector(
        (Epetra_DataAccess) CV, *Source, StartIndex, NumVectors));
}

void Epetra_MultiVector_Destroy ( 
  CT_Epetra_MultiVector_ID_t * selfID )
{
    CEpetra::removeMultiVector(selfID);
}

int Epetra_MultiVector_ReplaceGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceGlobalValue(GlobalRow, 
        VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceGlobalValue(GlobalBlockRow, 
        BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoGlobalValue(GlobalRow, 
        VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoGlobalValue(GlobalBlockRow, 
        BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceMyValue(MyRow, VectorIndex, 
        ScalarValue);
}

int Epetra_MultiVector_ReplaceMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceMyValue(MyBlockRow, 
        BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoMyValue(MyRow, VectorIndex, 
        ScalarValue);
}

int Epetra_MultiVector_SumIntoMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoMyValue(MyBlockRow, 
        BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_PutScalar ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarConstant )
{
    return CEpetra::getMultiVector(selfID)->PutScalar(ScalarConstant);
}

int Epetra_MultiVector_Random ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->Random();
}

int Epetra_MultiVector_ExtractCopy_Fill2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double * A, int MyLDA )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractCopy(A, MyLDA);
}

int Epetra_MultiVector_ExtractCopy_FillAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** ArrayOfPointers )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractCopy(ArrayOfPointers);
}

int Epetra_MultiVector_ExtractView_Set2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** A, int * MyLDA )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractView(A, MyLDA);
}

int Epetra_MultiVector_ExtractView_SetAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double *** ArrayOfPointers )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractView(ArrayOfPointers);
}

int Epetra_MultiVector_Dot ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID, 
  double * Result )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    return CEpetra::getConstMultiVector(selfID)->Dot(*A, Result);
}

int Epetra_MultiVector_Abs ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    return CEpetra::getMultiVector(selfID)->Abs(*A);
}

int Epetra_MultiVector_Reciprocal ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    return CEpetra::getMultiVector(selfID)->Reciprocal(*A);
}

int Epetra_MultiVector_Scale_Self ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->Scale(ScalarValue);
}

int Epetra_MultiVector_Scale ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    return CEpetra::getMultiVector(selfID)->Scale(ScalarA, *A);
}

int Epetra_MultiVector_Update_WithA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    return CEpetra::getMultiVector(selfID)->Update(ScalarA, *A, ScalarThis);
}

int Epetra_MultiVector_Update_WithAB ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarB, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CEpetra::getConstMultiVector(BID);
    return CEpetra::getMultiVector(selfID)->Update(ScalarA, *A, ScalarB, *B, 
        ScalarThis);
}

int Epetra_MultiVector_Norm1 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->Norm1(Result);
}

int Epetra_MultiVector_Norm2 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->Norm2(Result);
}

int Epetra_MultiVector_NormInf ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->NormInf(Result);
}

int Epetra_MultiVector_NormWeighted ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t WeightsID, double * Result )
{
    const Teuchos::RCP<const Epetra_MultiVector> Weights = 
        CEpetra::getConstMultiVector(WeightsID);
    return CEpetra::getConstMultiVector(selfID)->NormWeighted(*Weights, 
        Result);
}

int Epetra_MultiVector_MinValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->MinValue(Result);
}

int Epetra_MultiVector_MaxValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->MaxValue(Result);
}

int Epetra_MultiVector_MeanValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->MeanValue(Result);
}

int Epetra_MultiVector_Multiply_Matrix ( 
  CT_Epetra_MultiVector_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CEpetra::getConstMultiVector(BID);
    return CEpetra::getMultiVector(selfID)->Multiply(TransA, TransB, ScalarAB, 
        *A, *B, ScalarThis);
}

int Epetra_MultiVector_Multiply_ByEl ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CEpetra::getConstMultiVector(BID);
    return CEpetra::getMultiVector(selfID)->Multiply(ScalarAB, *A, *B, 
        ScalarThis);
}

int Epetra_MultiVector_ReciprocalMultiply ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    const Teuchos::RCP<const Epetra_MultiVector> A = 
        CEpetra::getConstMultiVector(AID);
    const Teuchos::RCP<const Epetra_MultiVector> B = 
        CEpetra::getConstMultiVector(BID);
    return CEpetra::getMultiVector(selfID)->ReciprocalMultiply(ScalarAB, *A, 
        *B, ScalarThis);
}

int Epetra_MultiVector_SetSeed ( 
  CT_Epetra_MultiVector_ID_t selfID, unsigned int Seed_in )
{
    return CEpetra::getMultiVector(selfID)->SetSeed(Seed_in);
}

unsigned int Epetra_MultiVector_Seed ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->Seed();
}

int Epetra_MultiVector_NumVectors ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->NumVectors();
}

int Epetra_MultiVector_MyLength ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->MyLength();
}

int Epetra_MultiVector_GlobalLength ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->GlobalLength();
}

int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->Stride();
}

boolean Epetra_MultiVector_ConstantStride ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return ((CEpetra::getConstMultiVector(
        selfID)->ConstantStride()) ? TRUE : FALSE);
}

int Epetra_MultiVector_ReplaceMap ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID )
{
    const Teuchos::RCP<const Epetra_BlockMap> map = CEpetra::getConstBlockMap(
        mapID);
    return CEpetra::getMultiVector(selfID)->ReplaceMap(*map);
}

void Epetra_MultiVector_Assign ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t SourceID )
{
    Epetra_MultiVector& self = *( CEpetra::getMultiVector(selfID) );

    const Teuchos::RCP<const Epetra_MultiVector> Source = 
        CEpetra::getConstMultiVector(SourceID);
    self = *Source;
}

double * Epetra_MultiVector_getArray ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    const Epetra_MultiVector& self = *( CEpetra::getConstMultiVector(
        selfID) );

    return self[i];
}

CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    const Epetra_MultiVector& self = *( CEpetra::getConstMultiVector(
        selfID) );

    return CEpetra::storeConstVector(self(i));
}


} // extern "C"


//
// Definitions from CEpetra_MultiVector_Cpp.hpp
//


/* get Epetra_MultiVector from non-const table using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    if (tableOfMultiVectors().isType(id.table))
        return tableOfMultiVectors().get<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
    else
        return CTrilinos::TableRepos::get<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
}

/* get Epetra_MultiVector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CTrilinos_Universal_ID_t id )
{
    if (tableOfMultiVectors().isType(id.table))
        return tableOfMultiVectors().get<Epetra_MultiVector>(id);
    else
        return CTrilinos::TableRepos::get<Epetra_MultiVector>(id);
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    if (tableOfMultiVectors().isType(id.table))
        return tableOfMultiVectors().getConst<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
    else
        return CTrilinos::TableRepos::getConst<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CTrilinos_Universal_ID_t id )
{
    if (tableOfMultiVectors().isType(id.table))
        return tableOfMultiVectors().getConst<Epetra_MultiVector>(id);
    else
        return CTrilinos::TableRepos::getConst<Epetra_MultiVector>(id);
}

/* store Epetra_MultiVector (owned) in non-const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeNewMultiVector( Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store<Epetra_MultiVector>(pobj, true));
}

/* store Epetra_MultiVector in non-const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeMultiVector( Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store<Epetra_MultiVector>(pobj, false));
}

/* store const Epetra_MultiVector in const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeConstMultiVector( const Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store<Epetra_MultiVector>(pobj, false));
}

/* remove Epetra_MultiVector from table using CT_Epetra_MultiVector_ID */
void
CEpetra::removeMultiVector( CT_Epetra_MultiVector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(*id);
    if (tableOfMultiVectors().isType(aid.table))
        tableOfMultiVectors().remove(&aid);
    else
        CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(aid);
}

/* remove Epetra_MultiVector from table using CTrilinos_Universal_ID_t */
void
CEpetra::removeMultiVector( CTrilinos_Universal_ID_t *aid )
{
    if (tableOfMultiVectors().isType(aid->table))
        tableOfMultiVectors().remove(aid);
    else
        CTrilinos::TableRepos::remove(aid);
}

/* purge Epetra_MultiVector table */
void
CEpetra::purgeMultiVector(  )
{
    tableOfMultiVectors().purge();
}

/* store Epetra_MultiVector in non-const table */
CTrilinos_Universal_ID_t
CEpetra::aliasMultiVector( const Teuchos::RCP< Epetra_MultiVector > & robj )
{
    return tableOfMultiVectors().alias(robj);
}

/* store const Epetra_MultiVector in const table */
CTrilinos_Universal_ID_t
CEpetra::aliasConstMultiVector( const Teuchos::RCP< const Epetra_MultiVector > & robj )
{
    return tableOfMultiVectors().alias(robj);
}



