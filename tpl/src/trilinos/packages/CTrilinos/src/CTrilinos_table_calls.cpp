
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


/*! @file CTrilinos_table_calls.cpp
 * @brief Calls to pull class instances out of the tables. */


#include "CTrilinos_config.h"
#include "CTrilinos_Table.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "CTrilinos_utils_templ.hpp"

//#include "CEpetra_Distributor_Cpp.hpp"
//#include "CEpetra_SerialComm_Cpp.hpp"
//#include "CEpetra_BLAS_Cpp.hpp"
//#include "CEpetra_Comm_Cpp.hpp"
//#include "CEpetra_Operator_Cpp.hpp"
//#include "CEpetra_MultiVector_Cpp.hpp"
//#include "CEpetra_OffsetIndex_Cpp.hpp"
//#include "CEpetra_Object_Cpp.hpp"
//#include "CEpetra_RowMatrix_Cpp.hpp"
//#include "CEpetra_CompObject_Cpp.hpp"
//#include "CEpetra_Directory_Cpp.hpp"
//#include "CEpetra_Flops_Cpp.hpp"
//#include "CEpetra_SrcDistObject_Cpp.hpp"
#ifdef HAVE_MPI
//#include "CEpetra_MpiComm_Cpp.hpp"
#endif /* HAVE_MPI */
//#include "CEpetra_CrsMatrix_Cpp.hpp"
//#include "CEpetra_CrsGraph_Cpp.hpp"
//#include "CEpetra_DistObject_Cpp.hpp"
//#include "CEpetra_Vector_Cpp.hpp"
//#include "CEpetra_Export_Cpp.hpp"
//#include "CEpetra_Map_Cpp.hpp"
//#include "CEpetra_BlockMap_Cpp.hpp"
//#include "CEpetra_Import_Cpp.hpp"
//#include "CEpetra_Time_Cpp.hpp"
//#include "CEpetra_JadMatrix_Cpp.hpp"
//#include "CEpetra_LinearProblem_Cpp.hpp"
//#include "CEpetra_LAPACK_Cpp.hpp"
//#include "CTeuchos_ParameterList_Cpp.hpp"
#ifdef HAVE_CTRILINOS_AMESOS
//#include "CAmesos_BaseSolver_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AMESOS */
//#include "CEpetra_FECrsMatrix_Cpp.hpp"
//#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
//#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#ifdef HAVE_CTRILINOS_AZTECOO
//#include "CAztecOO_StatusTest_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
//#include "CAztecOO_StatusTestCombo_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
//#include "CAztecOO_StatusTestMaxIters_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
//#include "CAztecOO_StatusTestResNorm_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
//#include "CIfpack_Preconditioner_Cpp.hpp"
#endif /* HAVE_CTRILINOS_IFPACK */
//#include "CEpetra_SerialDenseVector_Cpp.hpp"


//
// Definitions for Epetra_Distributor
//

namespace CEpetra {

/* get Epetra_Distributor from non-const table using CT_Epetra_Distributor_ID */
const Teuchos::RCP<Epetra_Distributor>
getDistributor( CT_Epetra_Distributor_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
}

/* get Epetra_Distributor from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Distributor>
getDistributor( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Distributor>(id);
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CT_Epetra_Distributor_ID */
const Teuchos::RCP<const Epetra_Distributor>
getConstDistributor( CT_Epetra_Distributor_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Distributor>(
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(id));
}

/* get const Epetra_Distributor from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Distributor>
getConstDistributor( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Distributor>(id);
}

/* store Epetra_Distributor (owned) in non-const table */
CT_Epetra_Distributor_ID_t
storeNewDistributor( Epetra_Distributor *pobj )
{
    CTrilinos::Table<Epetra_Distributor> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Distributor in non-const table */
CT_Epetra_Distributor_ID_t
storeDistributor( Epetra_Distributor *pobj )
{
    CTrilinos::Table<Epetra_Distributor> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Distributor in const table */
CT_Epetra_Distributor_ID_t
storeConstDistributor( const Epetra_Distributor *pobj )
{
    CTrilinos::Table<Epetra_Distributor> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Distributor from table using CT_Epetra_Distributor_ID */
void
removeDistributor( CT_Epetra_Distributor_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Distributor_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Distributor_ID_t>(aid);
}

/* purge Epetra_Distributor table */
void
purgeDistributor(  )
{
    CTrilinos::TableRepos::purge<Epetra_Distributor>();
}

}


//
// Definitions for Epetra_SerialComm
//

namespace CEpetra {

/* get Epetra_SerialComm from non-const table using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<Epetra_SerialComm>
getSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get Epetra_SerialComm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SerialComm>
getSerialComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SerialComm>(id);
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CT_Epetra_SerialComm_ID */
const Teuchos::RCP<const Epetra_SerialComm>
getConstSerialComm( CT_Epetra_SerialComm_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SerialComm>(
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(id));
}

/* get const Epetra_SerialComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SerialComm>
getConstSerialComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SerialComm>(id);
}

/* store Epetra_SerialComm (owned) in non-const table */
CT_Epetra_SerialComm_ID_t
storeNewSerialComm( Epetra_SerialComm *pobj )
{
    CTrilinos::Table<Epetra_SerialComm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_SerialComm in non-const table */
CT_Epetra_SerialComm_ID_t
storeSerialComm( Epetra_SerialComm *pobj )
{
    CTrilinos::Table<Epetra_SerialComm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_SerialComm in const table */
CT_Epetra_SerialComm_ID_t
storeConstSerialComm( const Epetra_SerialComm *pobj )
{
    CTrilinos::Table<Epetra_SerialComm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_SerialComm from table using CT_Epetra_SerialComm_ID */
void
removeSerialComm( CT_Epetra_SerialComm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SerialComm_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SerialComm_ID_t>(aid);
}

/* purge Epetra_SerialComm table */
void
purgeSerialComm(  )
{
    CTrilinos::TableRepos::purge<Epetra_SerialComm>();
}

}


//
// Definitions for Epetra_BLAS
//

namespace CEpetra {

/* get Epetra_BLAS from non-const table using CT_Epetra_BLAS_ID */
const Teuchos::RCP<Epetra_BLAS>
getBLAS( CT_Epetra_BLAS_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
}

/* get Epetra_BLAS from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_BLAS>
getBLAS( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_BLAS>(id);
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CT_Epetra_BLAS_ID */
const Teuchos::RCP<const Epetra_BLAS>
getConstBLAS( CT_Epetra_BLAS_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_BLAS>(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_BLAS>
getConstBLAS( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_BLAS>(id);
}

/* store Epetra_BLAS (owned) in non-const table */
CT_Epetra_BLAS_ID_t
storeNewBLAS( Epetra_BLAS *pobj )
{
    CTrilinos::Table<Epetra_BLAS> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_BLAS in non-const table */
CT_Epetra_BLAS_ID_t
storeBLAS( Epetra_BLAS *pobj )
{
    CTrilinos::Table<Epetra_BLAS> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_BLAS in const table */
CT_Epetra_BLAS_ID_t
storeConstBLAS( const Epetra_BLAS *pobj )
{
    CTrilinos::Table<Epetra_BLAS> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_BLAS from table using CT_Epetra_BLAS_ID */
void
removeBLAS( CT_Epetra_BLAS_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(aid);
}

/* purge Epetra_BLAS table */
void
purgeBLAS(  )
{
    CTrilinos::TableRepos::purge<Epetra_BLAS>();
}

}


//
// Definitions for Epetra_Comm
//

namespace CEpetra {

/* get Epetra_Comm from non-const table using CT_Epetra_Comm_ID */
const Teuchos::RCP<Epetra_Comm>
getComm( CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get Epetra_Comm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Comm>
getComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Comm>(id);
}

/* get const Epetra_Comm from either the const or non-const table
 * using CT_Epetra_Comm_ID */
const Teuchos::RCP<const Epetra_Comm>
getConstComm( CT_Epetra_Comm_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Comm>(
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(id));
}

/* get const Epetra_Comm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Comm>
getConstComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Comm>(id);
}

/* store Epetra_Comm (owned) in non-const table */
CT_Epetra_Comm_ID_t
storeNewComm( Epetra_Comm *pobj )
{
    CTrilinos::Table<Epetra_Comm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Comm in non-const table */
CT_Epetra_Comm_ID_t
storeComm( Epetra_Comm *pobj )
{
    CTrilinos::Table<Epetra_Comm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Comm in const table */
CT_Epetra_Comm_ID_t
storeConstComm( const Epetra_Comm *pobj )
{
    CTrilinos::Table<Epetra_Comm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Comm from table using CT_Epetra_Comm_ID */
void
removeComm( CT_Epetra_Comm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Comm_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Comm_ID_t>(aid);
}

/* purge Epetra_Comm table */
void
purgeComm(  )
{
    CTrilinos::TableRepos::purge<Epetra_Comm>();
}

}


//
// Definitions for Epetra_Operator
//

namespace CEpetra {

/* get Epetra_Operator from non-const table using CT_Epetra_Operator_ID */
const Teuchos::RCP<Epetra_Operator>
getOperator( CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
}

/* get Epetra_Operator from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Operator>
getOperator( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Operator>(id);
}

/* get const Epetra_Operator from either the const or non-const table
 * using CT_Epetra_Operator_ID */
const Teuchos::RCP<const Epetra_Operator>
getConstOperator( CT_Epetra_Operator_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Operator>(
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(id));
}

/* get const Epetra_Operator from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Operator>
getConstOperator( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Operator>(id);
}

/* store Epetra_Operator (owned) in non-const table */
CT_Epetra_Operator_ID_t
storeNewOperator( Epetra_Operator *pobj )
{
    CTrilinos::Table<Epetra_Operator> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Operator in non-const table */
CT_Epetra_Operator_ID_t
storeOperator( Epetra_Operator *pobj )
{
    CTrilinos::Table<Epetra_Operator> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Operator in const table */
CT_Epetra_Operator_ID_t
storeConstOperator( const Epetra_Operator *pobj )
{
    CTrilinos::Table<Epetra_Operator> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Operator from table using CT_Epetra_Operator_ID */
void
removeOperator( CT_Epetra_Operator_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Operator_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Operator_ID_t>(aid);
}

/* purge Epetra_Operator table */
void
purgeOperator(  )
{
    CTrilinos::TableRepos::purge<Epetra_Operator>();
}

}


//
// Definitions for Epetra_MultiVector
//

namespace CEpetra {

/* get Epetra_MultiVector from non-const table using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<Epetra_MultiVector>
getMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
}

/* get Epetra_MultiVector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_MultiVector>
getMultiVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_MultiVector>(id);
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<const Epetra_MultiVector>
getConstMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_MultiVector>(
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id));
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_MultiVector>
getConstMultiVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_MultiVector>(id);
}

/* store Epetra_MultiVector (owned) in non-const table */
CT_Epetra_MultiVector_ID_t
storeNewMultiVector( Epetra_MultiVector *pobj )
{
    CTrilinos::Table<Epetra_MultiVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_MultiVector in non-const table */
CT_Epetra_MultiVector_ID_t
storeMultiVector( Epetra_MultiVector *pobj )
{
    CTrilinos::Table<Epetra_MultiVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_MultiVector in const table */
CT_Epetra_MultiVector_ID_t
storeConstMultiVector( const Epetra_MultiVector *pobj )
{
    CTrilinos::Table<Epetra_MultiVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_MultiVector from table using CT_Epetra_MultiVector_ID */
void
removeMultiVector( CT_Epetra_MultiVector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(aid);
}

/* purge Epetra_MultiVector table */
void
purgeMultiVector(  )
{
    CTrilinos::TableRepos::purge<Epetra_MultiVector>();
}

}


//
// Definitions for Epetra_OffsetIndex
//

namespace CEpetra {

/* get Epetra_OffsetIndex from non-const table using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<Epetra_OffsetIndex>
getOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
}

/* get Epetra_OffsetIndex from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_OffsetIndex>
getOffsetIndex( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_OffsetIndex>(id);
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CT_Epetra_OffsetIndex_ID */
const Teuchos::RCP<const Epetra_OffsetIndex>
getConstOffsetIndex( CT_Epetra_OffsetIndex_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_OffsetIndex>(
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(id));
}

/* get const Epetra_OffsetIndex from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_OffsetIndex>
getConstOffsetIndex( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_OffsetIndex>(id);
}

/* store Epetra_OffsetIndex (owned) in non-const table */
CT_Epetra_OffsetIndex_ID_t
storeNewOffsetIndex( Epetra_OffsetIndex *pobj )
{
    CTrilinos::Table<Epetra_OffsetIndex> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_OffsetIndex in non-const table */
CT_Epetra_OffsetIndex_ID_t
storeOffsetIndex( Epetra_OffsetIndex *pobj )
{
    CTrilinos::Table<Epetra_OffsetIndex> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_OffsetIndex in const table */
CT_Epetra_OffsetIndex_ID_t
storeConstOffsetIndex( const Epetra_OffsetIndex *pobj )
{
    CTrilinos::Table<Epetra_OffsetIndex> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_OffsetIndex from table using CT_Epetra_OffsetIndex_ID */
void
removeOffsetIndex( CT_Epetra_OffsetIndex_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_OffsetIndex_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_OffsetIndex_ID_t>(aid);
}

/* purge Epetra_OffsetIndex table */
void
purgeOffsetIndex(  )
{
    CTrilinos::TableRepos::purge<Epetra_OffsetIndex>();
}

}


//
// Definitions for Epetra_Object
//

namespace CEpetra {

/* get Epetra_Object from non-const table using CT_Epetra_Object_ID */
const Teuchos::RCP<Epetra_Object>
getObject( CT_Epetra_Object_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
}

/* get Epetra_Object from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Object>
getObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Object>(id);
}

/* get const Epetra_Object from either the const or non-const table
 * using CT_Epetra_Object_ID */
const Teuchos::RCP<const Epetra_Object>
getConstObject( CT_Epetra_Object_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Object>(
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(id));
}

/* get const Epetra_Object from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Object>
getConstObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Object>(id);
}

/* store Epetra_Object (owned) in non-const table */
CT_Epetra_Object_ID_t
storeNewObject( Epetra_Object *pobj )
{
    CTrilinos::Table<Epetra_Object> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Object in non-const table */
CT_Epetra_Object_ID_t
storeObject( Epetra_Object *pobj )
{
    CTrilinos::Table<Epetra_Object> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Object in const table */
CT_Epetra_Object_ID_t
storeConstObject( const Epetra_Object *pobj )
{
    CTrilinos::Table<Epetra_Object> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Object_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Object from table using CT_Epetra_Object_ID */
void
removeObject( CT_Epetra_Object_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Object_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Object_ID_t>(aid);
}

/* purge Epetra_Object table */
void
purgeObject(  )
{
    CTrilinos::TableRepos::purge<Epetra_Object>();
}

}


//
// Definitions for Epetra_RowMatrix
//

namespace CEpetra {

/* get Epetra_RowMatrix from non-const table using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<Epetra_RowMatrix>
getRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
}

/* get Epetra_RowMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_RowMatrix>
getRowMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_RowMatrix>(id);
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CT_Epetra_RowMatrix_ID */
const Teuchos::RCP<const Epetra_RowMatrix>
getConstRowMatrix( CT_Epetra_RowMatrix_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_RowMatrix>(
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(id));
}

/* get const Epetra_RowMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_RowMatrix>
getConstRowMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_RowMatrix>(id);
}

/* store Epetra_RowMatrix (owned) in non-const table */
CT_Epetra_RowMatrix_ID_t
storeNewRowMatrix( Epetra_RowMatrix *pobj )
{
    CTrilinos::Table<Epetra_RowMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_RowMatrix in non-const table */
CT_Epetra_RowMatrix_ID_t
storeRowMatrix( Epetra_RowMatrix *pobj )
{
    CTrilinos::Table<Epetra_RowMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_RowMatrix in const table */
CT_Epetra_RowMatrix_ID_t
storeConstRowMatrix( const Epetra_RowMatrix *pobj )
{
    CTrilinos::Table<Epetra_RowMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_RowMatrix from table using CT_Epetra_RowMatrix_ID */
void
removeRowMatrix( CT_Epetra_RowMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_RowMatrix_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_RowMatrix_ID_t>(aid);
}

/* purge Epetra_RowMatrix table */
void
purgeRowMatrix(  )
{
    CTrilinos::TableRepos::purge<Epetra_RowMatrix>();
}

}


//
// Definitions for Epetra_CompObject
//

namespace CEpetra {

/* get Epetra_CompObject from non-const table using CT_Epetra_CompObject_ID */
const Teuchos::RCP<Epetra_CompObject>
getCompObject( CT_Epetra_CompObject_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
}

/* get Epetra_CompObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CompObject>
getCompObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_CompObject>(id);
}

/* get const Epetra_CompObject from either the const or non-const table
 * using CT_Epetra_CompObject_ID */
const Teuchos::RCP<const Epetra_CompObject>
getConstCompObject( CT_Epetra_CompObject_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_CompObject>(
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(id));
}

/* get const Epetra_CompObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CompObject>
getConstCompObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_CompObject>(id);
}

/* store Epetra_CompObject (owned) in non-const table */
CT_Epetra_CompObject_ID_t
storeNewCompObject( Epetra_CompObject *pobj )
{
    CTrilinos::Table<Epetra_CompObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_CompObject in non-const table */
CT_Epetra_CompObject_ID_t
storeCompObject( Epetra_CompObject *pobj )
{
    CTrilinos::Table<Epetra_CompObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_CompObject in const table */
CT_Epetra_CompObject_ID_t
storeConstCompObject( const Epetra_CompObject *pobj )
{
    CTrilinos::Table<Epetra_CompObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_CompObject from table using CT_Epetra_CompObject_ID */
void
removeCompObject( CT_Epetra_CompObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CompObject_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CompObject_ID_t>(aid);
}

/* purge Epetra_CompObject table */
void
purgeCompObject(  )
{
    CTrilinos::TableRepos::purge<Epetra_CompObject>();
}

}


//
// Definitions for Epetra_Directory
//

namespace CEpetra {

/* get Epetra_Directory from non-const table using CT_Epetra_Directory_ID */
const Teuchos::RCP<Epetra_Directory>
getDirectory( CT_Epetra_Directory_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
}

/* get Epetra_Directory from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Directory>
getDirectory( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Directory>(id);
}

/* get const Epetra_Directory from either the const or non-const table
 * using CT_Epetra_Directory_ID */
const Teuchos::RCP<const Epetra_Directory>
getConstDirectory( CT_Epetra_Directory_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Directory>(
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(id));
}

/* get const Epetra_Directory from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Directory>
getConstDirectory( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Directory>(id);
}

/* store Epetra_Directory (owned) in non-const table */
CT_Epetra_Directory_ID_t
storeNewDirectory( Epetra_Directory *pobj )
{
    CTrilinos::Table<Epetra_Directory> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Directory in non-const table */
CT_Epetra_Directory_ID_t
storeDirectory( Epetra_Directory *pobj )
{
    CTrilinos::Table<Epetra_Directory> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Directory in const table */
CT_Epetra_Directory_ID_t
storeConstDirectory( const Epetra_Directory *pobj )
{
    CTrilinos::Table<Epetra_Directory> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Directory from table using CT_Epetra_Directory_ID */
void
removeDirectory( CT_Epetra_Directory_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Directory_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Directory_ID_t>(aid);
}

/* purge Epetra_Directory table */
void
purgeDirectory(  )
{
    CTrilinos::TableRepos::purge<Epetra_Directory>();
}

}


//
// Definitions for Epetra_Flops
//

namespace CEpetra {

/* get Epetra_Flops from non-const table using CT_Epetra_Flops_ID */
const Teuchos::RCP<Epetra_Flops>
getFlops( CT_Epetra_Flops_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Flops>(
        CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id));
}

/* get Epetra_Flops from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Flops>
getFlops( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Flops>(id);
}

/* get const Epetra_Flops from either the const or non-const table
 * using CT_Epetra_Flops_ID */
const Teuchos::RCP<const Epetra_Flops>
getConstFlops( CT_Epetra_Flops_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Flops>(
        CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(id));
}

/* get const Epetra_Flops from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Flops>
getConstFlops( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Flops>(id);
}

/* store Epetra_Flops (owned) in non-const table */
CT_Epetra_Flops_ID_t
storeNewFlops( Epetra_Flops *pobj )
{
    CTrilinos::Table<Epetra_Flops> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Flops in non-const table */
CT_Epetra_Flops_ID_t
storeFlops( Epetra_Flops *pobj )
{
    CTrilinos::Table<Epetra_Flops> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Flops in const table */
CT_Epetra_Flops_ID_t
storeConstFlops( const Epetra_Flops *pobj )
{
    CTrilinos::Table<Epetra_Flops> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Flops from table using CT_Epetra_Flops_ID */
void
removeFlops( CT_Epetra_Flops_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Flops_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Flops_ID_t>(aid);
}

/* purge Epetra_Flops table */
void
purgeFlops(  )
{
    CTrilinos::TableRepos::purge<Epetra_Flops>();
}

}


//
// Definitions for Epetra_SrcDistObject
//

namespace CEpetra {

/* get Epetra_SrcDistObject from non-const table using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<Epetra_SrcDistObject>
getSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
}

/* get Epetra_SrcDistObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SrcDistObject>
getSrcDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SrcDistObject>(id);
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CT_Epetra_SrcDistObject_ID */
const Teuchos::RCP<const Epetra_SrcDistObject>
getConstSrcDistObject( CT_Epetra_SrcDistObject_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SrcDistObject>(
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(id));
}

/* get const Epetra_SrcDistObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SrcDistObject>
getConstSrcDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SrcDistObject>(id);
}

/* store Epetra_SrcDistObject (owned) in non-const table */
CT_Epetra_SrcDistObject_ID_t
storeNewSrcDistObject( Epetra_SrcDistObject *pobj )
{
    CTrilinos::Table<Epetra_SrcDistObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_SrcDistObject in non-const table */
CT_Epetra_SrcDistObject_ID_t
storeSrcDistObject( Epetra_SrcDistObject *pobj )
{
    CTrilinos::Table<Epetra_SrcDistObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_SrcDistObject in const table */
CT_Epetra_SrcDistObject_ID_t
storeConstSrcDistObject( const Epetra_SrcDistObject *pobj )
{
    CTrilinos::Table<Epetra_SrcDistObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_SrcDistObject from table using CT_Epetra_SrcDistObject_ID */
void
removeSrcDistObject( CT_Epetra_SrcDistObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SrcDistObject_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SrcDistObject_ID_t>(aid);
}

/* purge Epetra_SrcDistObject table */
void
purgeSrcDistObject(  )
{
    CTrilinos::TableRepos::purge<Epetra_SrcDistObject>();
}

}


//
// Definitions for Epetra_MpiComm
//

#ifdef HAVE_MPI

namespace CEpetra {

/* get Epetra_MpiComm from non-const table using CT_Epetra_MpiComm_ID */
const Teuchos::RCP<Epetra_MpiComm>
getMpiComm( CT_Epetra_MpiComm_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_MpiComm>(
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(id));
}

/* get Epetra_MpiComm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_MpiComm>
getMpiComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_MpiComm>(id);
}

/* get const Epetra_MpiComm from either the const or non-const table
 * using CT_Epetra_MpiComm_ID */
const Teuchos::RCP<const Epetra_MpiComm>
getConstMpiComm( CT_Epetra_MpiComm_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_MpiComm>(
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(id));
}

/* get const Epetra_MpiComm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_MpiComm>
getConstMpiComm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_MpiComm>(id);
}

/* store Epetra_MpiComm (owned) in non-const table */
CT_Epetra_MpiComm_ID_t
storeNewMpiComm( Epetra_MpiComm *pobj )
{
    CTrilinos::Table<Epetra_MpiComm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_MpiComm in non-const table */
CT_Epetra_MpiComm_ID_t
storeMpiComm( Epetra_MpiComm *pobj )
{
    CTrilinos::Table<Epetra_MpiComm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_MpiComm in const table */
CT_Epetra_MpiComm_ID_t
storeConstMpiComm( const Epetra_MpiComm *pobj )
{
    CTrilinos::Table<Epetra_MpiComm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_MpiComm from table using CT_Epetra_MpiComm_ID */
void
removeMpiComm( CT_Epetra_MpiComm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_MpiComm_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_MpiComm_ID_t>(aid);
}

/* purge Epetra_MpiComm table */
void
purgeMpiComm(  )
{
    CTrilinos::TableRepos::purge<Epetra_MpiComm>();
}

}

#endif /* HAVE_MPI */


//
// Definitions for Epetra_CrsMatrix
//

namespace CEpetra {

/* get Epetra_CrsMatrix from non-const table using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<Epetra_CrsMatrix>
getCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
}

/* get Epetra_CrsMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsMatrix>
getCrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_CrsMatrix>(id);
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CT_Epetra_CrsMatrix_ID */
const Teuchos::RCP<const Epetra_CrsMatrix>
getConstCrsMatrix( CT_Epetra_CrsMatrix_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_CrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(id));
}

/* get const Epetra_CrsMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsMatrix>
getConstCrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_CrsMatrix>(id);
}

/* store Epetra_CrsMatrix (owned) in non-const table */
CT_Epetra_CrsMatrix_ID_t
storeNewCrsMatrix( Epetra_CrsMatrix *pobj )
{
    CTrilinos::Table<Epetra_CrsMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_CrsMatrix in non-const table */
CT_Epetra_CrsMatrix_ID_t
storeCrsMatrix( Epetra_CrsMatrix *pobj )
{
    CTrilinos::Table<Epetra_CrsMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_CrsMatrix in const table */
CT_Epetra_CrsMatrix_ID_t
storeConstCrsMatrix( const Epetra_CrsMatrix *pobj )
{
    CTrilinos::Table<Epetra_CrsMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_CrsMatrix from table using CT_Epetra_CrsMatrix_ID */
void
removeCrsMatrix( CT_Epetra_CrsMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CrsMatrix_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CrsMatrix_ID_t>(aid);
}

/* purge Epetra_CrsMatrix table */
void
purgeCrsMatrix(  )
{
    CTrilinos::TableRepos::purge<Epetra_CrsMatrix>();
}

}


//
// Definitions for Epetra_CrsGraph
//

namespace CEpetra {

/* get Epetra_CrsGraph from non-const table using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<Epetra_CrsGraph>
getCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
}

/* get Epetra_CrsGraph from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_CrsGraph>
getCrsGraph( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_CrsGraph>(id);
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CT_Epetra_CrsGraph_ID */
const Teuchos::RCP<const Epetra_CrsGraph>
getConstCrsGraph( CT_Epetra_CrsGraph_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_CrsGraph>(
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(id));
}

/* get const Epetra_CrsGraph from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_CrsGraph>
getConstCrsGraph( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_CrsGraph>(id);
}

/* store Epetra_CrsGraph (owned) in non-const table */
CT_Epetra_CrsGraph_ID_t
storeNewCrsGraph( Epetra_CrsGraph *pobj )
{
    CTrilinos::Table<Epetra_CrsGraph> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_CrsGraph in non-const table */
CT_Epetra_CrsGraph_ID_t
storeCrsGraph( Epetra_CrsGraph *pobj )
{
    CTrilinos::Table<Epetra_CrsGraph> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_CrsGraph in const table */
CT_Epetra_CrsGraph_ID_t
storeConstCrsGraph( const Epetra_CrsGraph *pobj )
{
    CTrilinos::Table<Epetra_CrsGraph> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_CrsGraph from table using CT_Epetra_CrsGraph_ID */
void
removeCrsGraph( CT_Epetra_CrsGraph_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_CrsGraph_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_CrsGraph_ID_t>(aid);
}

/* purge Epetra_CrsGraph table */
void
purgeCrsGraph(  )
{
    CTrilinos::TableRepos::purge<Epetra_CrsGraph>();
}

}


//
// Definitions for Epetra_DistObject
//

namespace CEpetra {

/* get Epetra_DistObject from non-const table using CT_Epetra_DistObject_ID */
const Teuchos::RCP<Epetra_DistObject>
getDistObject( CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
}

/* get Epetra_DistObject from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_DistObject>
getDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_DistObject>(id);
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CT_Epetra_DistObject_ID */
const Teuchos::RCP<const Epetra_DistObject>
getConstDistObject( CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_DistObject>(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
}

/* get const Epetra_DistObject from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_DistObject>
getConstDistObject( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_DistObject>(id);
}

/* store Epetra_DistObject (owned) in non-const table */
CT_Epetra_DistObject_ID_t
storeNewDistObject( Epetra_DistObject *pobj )
{
    CTrilinos::Table<Epetra_DistObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_DistObject in non-const table */
CT_Epetra_DistObject_ID_t
storeDistObject( Epetra_DistObject *pobj )
{
    CTrilinos::Table<Epetra_DistObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_DistObject in const table */
CT_Epetra_DistObject_ID_t
storeConstDistObject( const Epetra_DistObject *pobj )
{
    CTrilinos::Table<Epetra_DistObject> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_DistObject from table using CT_Epetra_DistObject_ID */
void
removeDistObject( CT_Epetra_DistObject_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(aid);
}

/* purge Epetra_DistObject table */
void
purgeDistObject(  )
{
    CTrilinos::TableRepos::purge<Epetra_DistObject>();
}

}


//
// Definitions for Epetra_Vector
//

namespace CEpetra {

/* get Epetra_Vector from non-const table using CT_Epetra_Vector_ID */
const Teuchos::RCP<Epetra_Vector>
getVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Vector>(
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id));
}

/* get Epetra_Vector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Vector>
getVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Vector>(id);
}

/* get const Epetra_Vector from either the const or non-const table
 * using CT_Epetra_Vector_ID */
const Teuchos::RCP<const Epetra_Vector>
getConstVector( CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Vector>(
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id));
}

/* get const Epetra_Vector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Vector>
getConstVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Vector>(id);
}

/* store Epetra_Vector (owned) in non-const table */
CT_Epetra_Vector_ID_t
storeNewVector( Epetra_Vector *pobj )
{
    CTrilinos::Table<Epetra_Vector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Vector in non-const table */
CT_Epetra_Vector_ID_t
storeVector( Epetra_Vector *pobj )
{
    CTrilinos::Table<Epetra_Vector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Vector in const table */
CT_Epetra_Vector_ID_t
storeConstVector( const Epetra_Vector *pobj )
{
    CTrilinos::Table<Epetra_Vector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Vector from table using CT_Epetra_Vector_ID */
void
removeVector( CT_Epetra_Vector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(aid);
}

/* purge Epetra_Vector table */
void
purgeVector(  )
{
    CTrilinos::TableRepos::purge<Epetra_Vector>();
}

}


//
// Definitions for Epetra_Export
//

namespace CEpetra {

/* get Epetra_Export from non-const table using CT_Epetra_Export_ID */
const Teuchos::RCP<Epetra_Export>
getExport( CT_Epetra_Export_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
}

/* get Epetra_Export from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Export>
getExport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Export>(id);
}

/* get const Epetra_Export from either the const or non-const table
 * using CT_Epetra_Export_ID */
const Teuchos::RCP<const Epetra_Export>
getConstExport( CT_Epetra_Export_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Export>(
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(id));
}

/* get const Epetra_Export from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Export>
getConstExport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Export>(id);
}

/* store Epetra_Export (owned) in non-const table */
CT_Epetra_Export_ID_t
storeNewExport( Epetra_Export *pobj )
{
    CTrilinos::Table<Epetra_Export> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Export in non-const table */
CT_Epetra_Export_ID_t
storeExport( Epetra_Export *pobj )
{
    CTrilinos::Table<Epetra_Export> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Export in const table */
CT_Epetra_Export_ID_t
storeConstExport( const Epetra_Export *pobj )
{
    CTrilinos::Table<Epetra_Export> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Export_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Export from table using CT_Epetra_Export_ID */
void
removeExport( CT_Epetra_Export_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Export_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Export_ID_t>(aid);
}

/* purge Epetra_Export table */
void
purgeExport(  )
{
    CTrilinos::TableRepos::purge<Epetra_Export>();
}

}


//
// Definitions for Epetra_Map
//

namespace CEpetra {

/* get Epetra_Map from non-const table using CT_Epetra_Map_ID */
const Teuchos::RCP<Epetra_Map>
getMap( CT_Epetra_Map_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
}

/* get Epetra_Map from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Map>
getMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Map>(id);
}

/* get const Epetra_Map from either the const or non-const table
 * using CT_Epetra_Map_ID */
const Teuchos::RCP<const Epetra_Map>
getConstMap( CT_Epetra_Map_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Map>(
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(id));
}

/* get const Epetra_Map from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Map>
getConstMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Map>(id);
}

/* store Epetra_Map (owned) in non-const table */
CT_Epetra_Map_ID_t
storeNewMap( Epetra_Map *pobj )
{
    CTrilinos::Table<Epetra_Map> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Map in non-const table */
CT_Epetra_Map_ID_t
storeMap( Epetra_Map *pobj )
{
    CTrilinos::Table<Epetra_Map> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Map in const table */
CT_Epetra_Map_ID_t
storeConstMap( const Epetra_Map *pobj )
{
    CTrilinos::Table<Epetra_Map> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Map_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Map from table using CT_Epetra_Map_ID */
void
removeMap( CT_Epetra_Map_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Map_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Map_ID_t>(aid);
}

/* purge Epetra_Map table */
void
purgeMap(  )
{
    CTrilinos::TableRepos::purge<Epetra_Map>();
}

}


//
// Definitions for Epetra_BlockMap
//

namespace CEpetra {

/* get Epetra_BlockMap from non-const table using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<Epetra_BlockMap>
getBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
}

/* get Epetra_BlockMap from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_BlockMap>
getBlockMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_BlockMap>(id);
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CT_Epetra_BlockMap_ID */
const Teuchos::RCP<const Epetra_BlockMap>
getConstBlockMap( CT_Epetra_BlockMap_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_BlockMap>(
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(id));
}

/* get const Epetra_BlockMap from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_BlockMap>
getConstBlockMap( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_BlockMap>(id);
}

/* store Epetra_BlockMap (owned) in non-const table */
CT_Epetra_BlockMap_ID_t
storeNewBlockMap( Epetra_BlockMap *pobj )
{
    CTrilinos::Table<Epetra_BlockMap> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_BlockMap in non-const table */
CT_Epetra_BlockMap_ID_t
storeBlockMap( Epetra_BlockMap *pobj )
{
    CTrilinos::Table<Epetra_BlockMap> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_BlockMap in const table */
CT_Epetra_BlockMap_ID_t
storeConstBlockMap( const Epetra_BlockMap *pobj )
{
    CTrilinos::Table<Epetra_BlockMap> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_BlockMap from table using CT_Epetra_BlockMap_ID */
void
removeBlockMap( CT_Epetra_BlockMap_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_BlockMap_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_BlockMap_ID_t>(aid);
}

/* purge Epetra_BlockMap table */
void
purgeBlockMap(  )
{
    CTrilinos::TableRepos::purge<Epetra_BlockMap>();
}

}


//
// Definitions for Epetra_Import
//

namespace CEpetra {

/* get Epetra_Import from non-const table using CT_Epetra_Import_ID */
const Teuchos::RCP<Epetra_Import>
getImport( CT_Epetra_Import_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Import>(
        CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id));
}

/* get Epetra_Import from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Import>
getImport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Import>(id);
}

/* get const Epetra_Import from either the const or non-const table
 * using CT_Epetra_Import_ID */
const Teuchos::RCP<const Epetra_Import>
getConstImport( CT_Epetra_Import_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Import>(
        CTrilinos::abstractType<CT_Epetra_Import_ID_t>(id));
}

/* get const Epetra_Import from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Import>
getConstImport( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Import>(id);
}

/* store Epetra_Import (owned) in non-const table */
CT_Epetra_Import_ID_t
storeNewImport( Epetra_Import *pobj )
{
    CTrilinos::Table<Epetra_Import> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Import in non-const table */
CT_Epetra_Import_ID_t
storeImport( Epetra_Import *pobj )
{
    CTrilinos::Table<Epetra_Import> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Import in const table */
CT_Epetra_Import_ID_t
storeConstImport( const Epetra_Import *pobj )
{
    CTrilinos::Table<Epetra_Import> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Import_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Import from table using CT_Epetra_Import_ID */
void
removeImport( CT_Epetra_Import_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Import_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Import_ID_t>(aid);
}

/* purge Epetra_Import table */
void
purgeImport(  )
{
    CTrilinos::TableRepos::purge<Epetra_Import>();
}

}


//
// Definitions for Epetra_Time
//

namespace CEpetra {

/* get Epetra_Time from non-const table using CT_Epetra_Time_ID */
const Teuchos::RCP<Epetra_Time>
getTime( CT_Epetra_Time_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
}

/* get Epetra_Time from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_Time>
getTime( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_Time>(id);
}

/* get const Epetra_Time from either the const or non-const table
 * using CT_Epetra_Time_ID */
const Teuchos::RCP<const Epetra_Time>
getConstTime( CT_Epetra_Time_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Time>(
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(id));
}

/* get const Epetra_Time from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_Time>
getConstTime( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_Time>(id);
}

/* store Epetra_Time (owned) in non-const table */
CT_Epetra_Time_ID_t
storeNewTime( Epetra_Time *pobj )
{
    CTrilinos::Table<Epetra_Time> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_Time in non-const table */
CT_Epetra_Time_ID_t
storeTime( Epetra_Time *pobj )
{
    CTrilinos::Table<Epetra_Time> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_Time in const table */
CT_Epetra_Time_ID_t
storeConstTime( const Epetra_Time *pobj )
{
    CTrilinos::Table<Epetra_Time> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_Time_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_Time from table using CT_Epetra_Time_ID */
void
removeTime( CT_Epetra_Time_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_Time_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_Time_ID_t>(aid);
}

/* purge Epetra_Time table */
void
purgeTime(  )
{
    CTrilinos::TableRepos::purge<Epetra_Time>();
}

}


//
// Definitions for Epetra_JadMatrix
//

namespace CEpetra {

/* get Epetra_JadMatrix from non-const table using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<Epetra_JadMatrix>
getJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
}

/* get Epetra_JadMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_JadMatrix>
getJadMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_JadMatrix>(id);
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CT_Epetra_JadMatrix_ID */
const Teuchos::RCP<const Epetra_JadMatrix>
getConstJadMatrix( CT_Epetra_JadMatrix_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_JadMatrix>(
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(id));
}

/* get const Epetra_JadMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_JadMatrix>
getConstJadMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_JadMatrix>(id);
}

/* store Epetra_JadMatrix (owned) in non-const table */
CT_Epetra_JadMatrix_ID_t
storeNewJadMatrix( Epetra_JadMatrix *pobj )
{
    CTrilinos::Table<Epetra_JadMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_JadMatrix in non-const table */
CT_Epetra_JadMatrix_ID_t
storeJadMatrix( Epetra_JadMatrix *pobj )
{
    CTrilinos::Table<Epetra_JadMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_JadMatrix in const table */
CT_Epetra_JadMatrix_ID_t
storeConstJadMatrix( const Epetra_JadMatrix *pobj )
{
    CTrilinos::Table<Epetra_JadMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_JadMatrix from table using CT_Epetra_JadMatrix_ID */
void
removeJadMatrix( CT_Epetra_JadMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_JadMatrix_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_JadMatrix_ID_t>(aid);
}

/* purge Epetra_JadMatrix table */
void
purgeJadMatrix(  )
{
    CTrilinos::TableRepos::purge<Epetra_JadMatrix>();
}

}


//
// Definitions for Epetra_LinearProblem
//

namespace CEpetra {

/* get Epetra_LinearProblem from non-const table using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<Epetra_LinearProblem>
getLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get Epetra_LinearProblem from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_LinearProblem>
getLinearProblem( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_LinearProblem>(id);
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CT_Epetra_LinearProblem_ID */
const Teuchos::RCP<const Epetra_LinearProblem>
getConstLinearProblem( CT_Epetra_LinearProblem_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_LinearProblem>(
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(id));
}

/* get const Epetra_LinearProblem from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_LinearProblem>
getConstLinearProblem( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_LinearProblem>(id);
}

/* store Epetra_LinearProblem (owned) in non-const table */
CT_Epetra_LinearProblem_ID_t
storeNewLinearProblem( Epetra_LinearProblem *pobj )
{
    CTrilinos::Table<Epetra_LinearProblem> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_LinearProblem in non-const table */
CT_Epetra_LinearProblem_ID_t
storeLinearProblem( Epetra_LinearProblem *pobj )
{
    CTrilinos::Table<Epetra_LinearProblem> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_LinearProblem in const table */
CT_Epetra_LinearProblem_ID_t
storeConstLinearProblem( const Epetra_LinearProblem *pobj )
{
    CTrilinos::Table<Epetra_LinearProblem> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_LinearProblem from table using CT_Epetra_LinearProblem_ID */
void
removeLinearProblem( CT_Epetra_LinearProblem_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_LinearProblem_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_LinearProblem_ID_t>(aid);
}

/* purge Epetra_LinearProblem table */
void
purgeLinearProblem(  )
{
    CTrilinos::TableRepos::purge<Epetra_LinearProblem>();
}

}


//
// Definitions for Epetra_LAPACK
//

namespace CEpetra {

/* get Epetra_LAPACK from non-const table using CT_Epetra_LAPACK_ID */
const Teuchos::RCP<Epetra_LAPACK>
getLAPACK( CT_Epetra_LAPACK_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
}

/* get Epetra_LAPACK from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_LAPACK>
getLAPACK( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_LAPACK>(id);
}

/* get const Epetra_LAPACK from either the const or non-const table
 * using CT_Epetra_LAPACK_ID */
const Teuchos::RCP<const Epetra_LAPACK>
getConstLAPACK( CT_Epetra_LAPACK_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_LAPACK>(
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(id));
}

/* get const Epetra_LAPACK from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_LAPACK>
getConstLAPACK( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_LAPACK>(id);
}

/* store Epetra_LAPACK (owned) in non-const table */
CT_Epetra_LAPACK_ID_t
storeNewLAPACK( Epetra_LAPACK *pobj )
{
    CTrilinos::Table<Epetra_LAPACK> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_LAPACK in non-const table */
CT_Epetra_LAPACK_ID_t
storeLAPACK( Epetra_LAPACK *pobj )
{
    CTrilinos::Table<Epetra_LAPACK> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_LAPACK in const table */
CT_Epetra_LAPACK_ID_t
storeConstLAPACK( const Epetra_LAPACK *pobj )
{
    CTrilinos::Table<Epetra_LAPACK> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_LAPACK from table using CT_Epetra_LAPACK_ID */
void
removeLAPACK( CT_Epetra_LAPACK_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_LAPACK_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_LAPACK_ID_t>(aid);
}

/* purge Epetra_LAPACK table */
void
purgeLAPACK(  )
{
    CTrilinos::TableRepos::purge<Epetra_LAPACK>();
}

}


//
// Definitions for ParameterList
//

namespace CTeuchos {

/* get Teuchos::ParameterList from non-const table using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<Teuchos::ParameterList>
getParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::TableRepos::get<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
}

/* get Teuchos::ParameterList from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Teuchos::ParameterList>
getParameterList( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Teuchos::ParameterList>(id);
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CT_Teuchos_ParameterList_ID */
const Teuchos::RCP<const Teuchos::ParameterList>
getConstParameterList( CT_Teuchos_ParameterList_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Teuchos::ParameterList>(
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(id));
}

/* get const Teuchos::ParameterList from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Teuchos::ParameterList>
getConstParameterList( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Teuchos::ParameterList>(id);
}

/* store Teuchos::ParameterList (owned) in non-const table */
CT_Teuchos_ParameterList_ID_t
storeNewParameterList( Teuchos::ParameterList *pobj )
{
    CTrilinos::Table<Teuchos::ParameterList> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tab->store(pobj, true));
}

/* store Teuchos::ParameterList in non-const table */
CT_Teuchos_ParameterList_ID_t
storeParameterList( Teuchos::ParameterList *pobj )
{
    CTrilinos::Table<Teuchos::ParameterList> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tab->store(pobj, false));
}

/* store const Teuchos::ParameterList in const table */
CT_Teuchos_ParameterList_ID_t
storeConstParameterList( const Teuchos::ParameterList *pobj )
{
    CTrilinos::Table<Teuchos::ParameterList> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(
        tab->store(pobj, false));
}

/* remove Teuchos::ParameterList from table using CT_Teuchos_ParameterList_ID */
void
removeParameterList( CT_Teuchos_ParameterList_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Teuchos_ParameterList_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Teuchos_ParameterList_ID_t>(aid);
}

/* purge Teuchos::ParameterList table */
void
purgeParameterList(  )
{
    CTrilinos::TableRepos::purge<Teuchos::ParameterList>();
}

}


//
// Definitions for Amesos_BaseSolver
//

#ifdef HAVE_CTRILINOS_AMESOS

namespace CAmesos {

/* get Amesos_BaseSolver from non-const table using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<Amesos_BaseSolver>
getBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    return CTrilinos::TableRepos::get<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
}

/* get Amesos_BaseSolver from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Amesos_BaseSolver>
getBaseSolver( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Amesos_BaseSolver>(id);
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CT_Amesos_BaseSolver_ID */
const Teuchos::RCP<const Amesos_BaseSolver>
getConstBaseSolver( CT_Amesos_BaseSolver_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Amesos_BaseSolver>(
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(id));
}

/* get const Amesos_BaseSolver from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Amesos_BaseSolver>
getConstBaseSolver( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Amesos_BaseSolver>(id);
}

/* store Amesos_BaseSolver (owned) in non-const table */
CT_Amesos_BaseSolver_ID_t
storeNewBaseSolver( Amesos_BaseSolver *pobj )
{
    CTrilinos::Table<Amesos_BaseSolver> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        tab->store(pobj, true));
}

/* store Amesos_BaseSolver in non-const table */
CT_Amesos_BaseSolver_ID_t
storeBaseSolver( Amesos_BaseSolver *pobj )
{
    CTrilinos::Table<Amesos_BaseSolver> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        tab->store(pobj, false));
}

/* store const Amesos_BaseSolver in const table */
CT_Amesos_BaseSolver_ID_t
storeConstBaseSolver( const Amesos_BaseSolver *pobj )
{
    CTrilinos::Table<Amesos_BaseSolver> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(
        tab->store(pobj, false));
}

/* remove Amesos_BaseSolver from table using CT_Amesos_BaseSolver_ID */
void
removeBaseSolver( CT_Amesos_BaseSolver_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Amesos_BaseSolver_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Amesos_BaseSolver_ID_t>(aid);
}

/* purge Amesos_BaseSolver table */
void
purgeBaseSolver(  )
{
    CTrilinos::TableRepos::purge<Amesos_BaseSolver>();
}

}

#endif /* HAVE_CTRILINOS_AMESOS */


//
// Definitions for Epetra_FECrsMatrix
//

namespace CEpetra {

/* get Epetra_FECrsMatrix from non-const table using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<Epetra_FECrsMatrix>
getFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
}

/* get Epetra_FECrsMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_FECrsMatrix>
getFECrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_FECrsMatrix>(id);
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CT_Epetra_FECrsMatrix_ID */
const Teuchos::RCP<const Epetra_FECrsMatrix>
getConstFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_FECrsMatrix>(
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(id));
}

/* get const Epetra_FECrsMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_FECrsMatrix>
getConstFECrsMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_FECrsMatrix>(id);
}

/* store Epetra_FECrsMatrix (owned) in non-const table */
CT_Epetra_FECrsMatrix_ID_t
storeNewFECrsMatrix( Epetra_FECrsMatrix *pobj )
{
    CTrilinos::Table<Epetra_FECrsMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_FECrsMatrix in non-const table */
CT_Epetra_FECrsMatrix_ID_t
storeFECrsMatrix( Epetra_FECrsMatrix *pobj )
{
    CTrilinos::Table<Epetra_FECrsMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_FECrsMatrix in const table */
CT_Epetra_FECrsMatrix_ID_t
storeConstFECrsMatrix( const Epetra_FECrsMatrix *pobj )
{
    CTrilinos::Table<Epetra_FECrsMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_FECrsMatrix from table using CT_Epetra_FECrsMatrix_ID */
void
removeFECrsMatrix( CT_Epetra_FECrsMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_FECrsMatrix_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_FECrsMatrix_ID_t>(aid);
}

/* purge Epetra_FECrsMatrix table */
void
purgeFECrsMatrix(  )
{
    CTrilinos::TableRepos::purge<Epetra_FECrsMatrix>();
}

}


//
// Definitions for Epetra_IntSerialDenseVector
//

namespace CEpetra {

/* get Epetra_IntSerialDenseVector from non-const table using CT_Epetra_IntSerialDenseVector_ID */
const Teuchos::RCP<Epetra_IntSerialDenseVector>
getIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_IntSerialDenseVector>(
        CTrilinos::abstractType<CT_Epetra_IntSerialDenseVector_ID_t>(id));
}

/* get Epetra_IntSerialDenseVector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_IntSerialDenseVector>
getIntSerialDenseVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_IntSerialDenseVector>(id);
}

/* get const Epetra_IntSerialDenseVector from either the const or non-const table
 * using CT_Epetra_IntSerialDenseVector_ID */
const Teuchos::RCP<const Epetra_IntSerialDenseVector>
getConstIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_IntSerialDenseVector>(
        CTrilinos::abstractType<CT_Epetra_IntSerialDenseVector_ID_t>(id));
}

/* get const Epetra_IntSerialDenseVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_IntSerialDenseVector>
getConstIntSerialDenseVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_IntSerialDenseVector>(id);
}

/* store Epetra_IntSerialDenseVector (owned) in non-const table */
CT_Epetra_IntSerialDenseVector_ID_t
storeNewIntSerialDenseVector( Epetra_IntSerialDenseVector *pobj )
{
    CTrilinos::Table<Epetra_IntSerialDenseVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_IntSerialDenseVector in non-const table */
CT_Epetra_IntSerialDenseVector_ID_t
storeIntSerialDenseVector( Epetra_IntSerialDenseVector *pobj )
{
    CTrilinos::Table<Epetra_IntSerialDenseVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_IntSerialDenseVector in const table */
CT_Epetra_IntSerialDenseVector_ID_t
storeConstIntSerialDenseVector( const Epetra_IntSerialDenseVector *pobj )
{
    CTrilinos::Table<Epetra_IntSerialDenseVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_IntSerialDenseVector from table using CT_Epetra_IntSerialDenseVector_ID */
void
removeIntSerialDenseVector( CT_Epetra_IntSerialDenseVector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_IntSerialDenseVector_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_IntSerialDenseVector_ID_t>(aid);
}

/* purge Epetra_IntSerialDenseVector table */
void
purgeIntSerialDenseVector(  )
{
    CTrilinos::TableRepos::purge<Epetra_IntSerialDenseVector>();
}

}


//
// Definitions for Epetra_SerialDenseMatrix
//

namespace CEpetra {

/* get Epetra_SerialDenseMatrix from non-const table using CT_Epetra_SerialDenseMatrix_ID */
const Teuchos::RCP<Epetra_SerialDenseMatrix>
getSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SerialDenseMatrix>(
        CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id));
}

/* get Epetra_SerialDenseMatrix from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SerialDenseMatrix>
getSerialDenseMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SerialDenseMatrix>(id);
}

/* get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CT_Epetra_SerialDenseMatrix_ID */
const Teuchos::RCP<const Epetra_SerialDenseMatrix>
getConstSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SerialDenseMatrix>(
        CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(id));
}

/* get const Epetra_SerialDenseMatrix from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SerialDenseMatrix>
getConstSerialDenseMatrix( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SerialDenseMatrix>(id);
}

/* store Epetra_SerialDenseMatrix (owned) in non-const table */
CT_Epetra_SerialDenseMatrix_ID_t
storeNewSerialDenseMatrix( Epetra_SerialDenseMatrix *pobj )
{
    CTrilinos::Table<Epetra_SerialDenseMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_SerialDenseMatrix in non-const table */
CT_Epetra_SerialDenseMatrix_ID_t
storeSerialDenseMatrix( Epetra_SerialDenseMatrix *pobj )
{
    CTrilinos::Table<Epetra_SerialDenseMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_SerialDenseMatrix in const table */
CT_Epetra_SerialDenseMatrix_ID_t
storeConstSerialDenseMatrix( const Epetra_SerialDenseMatrix *pobj )
{
    CTrilinos::Table<Epetra_SerialDenseMatrix> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_SerialDenseMatrix from table using CT_Epetra_SerialDenseMatrix_ID */
void
removeSerialDenseMatrix( CT_Epetra_SerialDenseMatrix_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SerialDenseMatrix_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SerialDenseMatrix_ID_t>(aid);
}

/* purge Epetra_SerialDenseMatrix table */
void
purgeSerialDenseMatrix(  )
{
    CTrilinos::TableRepos::purge<Epetra_SerialDenseMatrix>();
}

}


//
// Definitions for AztecOO_StatusTest
//

#ifdef HAVE_CTRILINOS_AZTECOO

namespace CAztecOO {

/* get AztecOO_StatusTest from non-const table using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<AztecOO_StatusTest>
getStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get AztecOO_StatusTest from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTest>
getStatusTest( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTest>(id);
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CT_AztecOO_StatusTest_ID */
const Teuchos::RCP<const AztecOO_StatusTest>
getConstStatusTest( CT_AztecOO_StatusTest_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTest>(
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(id));
}

/* get const AztecOO_StatusTest from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTest>
getConstStatusTest( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTest>(id);
}

/* store AztecOO_StatusTest (owned) in non-const table */
CT_AztecOO_StatusTest_ID_t
storeNewStatusTest( AztecOO_StatusTest *pobj )
{
    CTrilinos::Table<AztecOO_StatusTest> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        tab->store(pobj, true));
}

/* store AztecOO_StatusTest in non-const table */
CT_AztecOO_StatusTest_ID_t
storeStatusTest( AztecOO_StatusTest *pobj )
{
    CTrilinos::Table<AztecOO_StatusTest> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        tab->store(pobj, false));
}

/* store const AztecOO_StatusTest in const table */
CT_AztecOO_StatusTest_ID_t
storeConstStatusTest( const AztecOO_StatusTest *pobj )
{
    CTrilinos::Table<AztecOO_StatusTest> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(
        tab->store(pobj, false));
}

/* remove AztecOO_StatusTest from table using CT_AztecOO_StatusTest_ID */
void
removeStatusTest( CT_AztecOO_StatusTest_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTest_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTest_ID_t>(aid);
}

/* purge AztecOO_StatusTest table */
void
purgeStatusTest(  )
{
    CTrilinos::TableRepos::purge<AztecOO_StatusTest>();
}

}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for AztecOO_StatusTestCombo
//

#ifdef HAVE_CTRILINOS_AZTECOO

namespace CAztecOO {

/* get AztecOO_StatusTestCombo from non-const table using CT_AztecOO_StatusTestCombo_ID */
const Teuchos::RCP<AztecOO_StatusTestCombo>
getStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTestCombo>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id));
}

/* get AztecOO_StatusTestCombo from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestCombo>
getStatusTestCombo( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTestCombo>(id);
}

/* get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CT_AztecOO_StatusTestCombo_ID */
const Teuchos::RCP<const AztecOO_StatusTestCombo>
getConstStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTestCombo>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(id));
}

/* get const AztecOO_StatusTestCombo from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestCombo>
getConstStatusTestCombo( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTestCombo>(id);
}

/* store AztecOO_StatusTestCombo (owned) in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
storeNewStatusTestCombo( AztecOO_StatusTestCombo *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestCombo> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        tab->store(pobj, true));
}

/* store AztecOO_StatusTestCombo in non-const table */
CT_AztecOO_StatusTestCombo_ID_t
storeStatusTestCombo( AztecOO_StatusTestCombo *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestCombo> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        tab->store(pobj, false));
}

/* store const AztecOO_StatusTestCombo in const table */
CT_AztecOO_StatusTestCombo_ID_t
storeConstStatusTestCombo( const AztecOO_StatusTestCombo *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestCombo> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(
        tab->store(pobj, false));
}

/* remove AztecOO_StatusTestCombo from table using CT_AztecOO_StatusTestCombo_ID */
void
removeStatusTestCombo( CT_AztecOO_StatusTestCombo_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestCombo_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestCombo_ID_t>(aid);
}

/* purge AztecOO_StatusTestCombo table */
void
purgeStatusTestCombo(  )
{
    CTrilinos::TableRepos::purge<AztecOO_StatusTestCombo>();
}

}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for AztecOO_StatusTestMaxIters
//

#ifdef HAVE_CTRILINOS_AZTECOO

namespace CAztecOO {

/* get AztecOO_StatusTestMaxIters from non-const table using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
getStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get AztecOO_StatusTestMaxIters from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestMaxIters>
getStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTestMaxIters>(id);
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CT_AztecOO_StatusTestMaxIters_ID */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
getConstStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTestMaxIters>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(id));
}

/* get const AztecOO_StatusTestMaxIters from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestMaxIters>
getConstStatusTestMaxIters( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTestMaxIters>(id);
}

/* store AztecOO_StatusTestMaxIters (owned) in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
storeNewStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestMaxIters> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        tab->store(pobj, true));
}

/* store AztecOO_StatusTestMaxIters in non-const table */
CT_AztecOO_StatusTestMaxIters_ID_t
storeStatusTestMaxIters( AztecOO_StatusTestMaxIters *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestMaxIters> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        tab->store(pobj, false));
}

/* store const AztecOO_StatusTestMaxIters in const table */
CT_AztecOO_StatusTestMaxIters_ID_t
storeConstStatusTestMaxIters( const AztecOO_StatusTestMaxIters *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestMaxIters> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(
        tab->store(pobj, false));
}

/* remove AztecOO_StatusTestMaxIters from table using CT_AztecOO_StatusTestMaxIters_ID */
void
removeStatusTestMaxIters( CT_AztecOO_StatusTestMaxIters_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestMaxIters_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestMaxIters_ID_t>(aid);
}

/* purge AztecOO_StatusTestMaxIters table */
void
purgeStatusTestMaxIters(  )
{
    CTrilinos::TableRepos::purge<AztecOO_StatusTestMaxIters>();
}

}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for AztecOO_StatusTestResNorm
//

#ifdef HAVE_CTRILINOS_AZTECOO

namespace CAztecOO {

/* get AztecOO_StatusTestResNorm from non-const table using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
getStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
}

/* get AztecOO_StatusTestResNorm from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<AztecOO_StatusTestResNorm>
getStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<AztecOO_StatusTestResNorm>(id);
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CT_AztecOO_StatusTestResNorm_ID */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
getConstStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTestResNorm>(
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(id));
}

/* get const AztecOO_StatusTestResNorm from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const AztecOO_StatusTestResNorm>
getConstStatusTestResNorm( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<AztecOO_StatusTestResNorm>(id);
}

/* store AztecOO_StatusTestResNorm (owned) in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
storeNewStatusTestResNorm( AztecOO_StatusTestResNorm *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestResNorm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        tab->store(pobj, true));
}

/* store AztecOO_StatusTestResNorm in non-const table */
CT_AztecOO_StatusTestResNorm_ID_t
storeStatusTestResNorm( AztecOO_StatusTestResNorm *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestResNorm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        tab->store(pobj, false));
}

/* store const AztecOO_StatusTestResNorm in const table */
CT_AztecOO_StatusTestResNorm_ID_t
storeConstStatusTestResNorm( const AztecOO_StatusTestResNorm *pobj )
{
    CTrilinos::Table<AztecOO_StatusTestResNorm> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(
        tab->store(pobj, false));
}

/* remove AztecOO_StatusTestResNorm from table using CT_AztecOO_StatusTestResNorm_ID */
void
removeStatusTestResNorm( CT_AztecOO_StatusTestResNorm_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_AztecOO_StatusTestResNorm_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_AztecOO_StatusTestResNorm_ID_t>(aid);
}

/* purge AztecOO_StatusTestResNorm table */
void
purgeStatusTestResNorm(  )
{
    CTrilinos::TableRepos::purge<AztecOO_StatusTestResNorm>();
}

}

#endif /* HAVE_CTRILINOS_AZTECOO */


//
// Definitions for Ifpack_Preconditioner
//

#ifdef HAVE_CTRILINOS_IFPACK

namespace CIfpack {

/* get Ifpack_Preconditioner from non-const table using CT_Ifpack_Preconditioner_ID */
const Teuchos::RCP<Ifpack_Preconditioner>
getPreconditioner( CT_Ifpack_Preconditioner_ID_t id )
{
    return CTrilinos::TableRepos::get<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
}

/* get Ifpack_Preconditioner from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Ifpack_Preconditioner>
getPreconditioner( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Ifpack_Preconditioner>(id);
}

/* get const Ifpack_Preconditioner from either the const or non-const table
 * using CT_Ifpack_Preconditioner_ID */
const Teuchos::RCP<const Ifpack_Preconditioner>
getConstPreconditioner( CT_Ifpack_Preconditioner_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Ifpack_Preconditioner>(
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(id));
}

/* get const Ifpack_Preconditioner from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Ifpack_Preconditioner>
getConstPreconditioner( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Ifpack_Preconditioner>(id);
}

/* store Ifpack_Preconditioner (owned) in non-const table */
CT_Ifpack_Preconditioner_ID_t
storeNewPreconditioner( Ifpack_Preconditioner *pobj )
{
    CTrilinos::Table<Ifpack_Preconditioner> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        tab->store(pobj, true));
}

/* store Ifpack_Preconditioner in non-const table */
CT_Ifpack_Preconditioner_ID_t
storePreconditioner( Ifpack_Preconditioner *pobj )
{
    CTrilinos::Table<Ifpack_Preconditioner> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        tab->store(pobj, false));
}

/* store const Ifpack_Preconditioner in const table */
CT_Ifpack_Preconditioner_ID_t
storeConstPreconditioner( const Ifpack_Preconditioner *pobj )
{
    CTrilinos::Table<Ifpack_Preconditioner> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(
        tab->store(pobj, false));
}

/* remove Ifpack_Preconditioner from table using CT_Ifpack_Preconditioner_ID */
void
removePreconditioner( CT_Ifpack_Preconditioner_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Ifpack_Preconditioner_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Ifpack_Preconditioner_ID_t>(aid);
}

/* purge Ifpack_Preconditioner table */
void
purgePreconditioner(  )
{
    CTrilinos::TableRepos::purge<Ifpack_Preconditioner>();
}

}

#endif /* HAVE_CTRILINOS_IFPACK */


//
// Definitions for Epetra_SerialDenseVector
//

namespace CEpetra {

/* get Epetra_SerialDenseVector from non-const table using CT_Epetra_SerialDenseVector_ID */
const Teuchos::RCP<Epetra_SerialDenseVector>
getSerialDenseVector( CT_Epetra_SerialDenseVector_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SerialDenseVector>(
        CTrilinos::abstractType<CT_Epetra_SerialDenseVector_ID_t>(id));
}

/* get Epetra_SerialDenseVector from non-const table using CTrilinos_Universal_ID_t */
const Teuchos::RCP<Epetra_SerialDenseVector>
getSerialDenseVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::get<Epetra_SerialDenseVector>(id);
}

/* get const Epetra_SerialDenseVector from either the const or non-const table
 * using CT_Epetra_SerialDenseVector_ID */
const Teuchos::RCP<const Epetra_SerialDenseVector>
getConstSerialDenseVector( CT_Epetra_SerialDenseVector_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SerialDenseVector>(
        CTrilinos::abstractType<CT_Epetra_SerialDenseVector_ID_t>(id));
}

/* get const Epetra_SerialDenseVector from either the const or non-const table
 * using CTrilinos_Universal_ID_t */
const Teuchos::RCP<const Epetra_SerialDenseVector>
getConstSerialDenseVector( CTrilinos_Universal_ID_t id )
{
    return CTrilinos::TableRepos::getConst<Epetra_SerialDenseVector>(id);
}

/* store Epetra_SerialDenseVector (owned) in non-const table */
CT_Epetra_SerialDenseVector_ID_t
storeNewSerialDenseVector( Epetra_SerialDenseVector *pobj )
{
    CTrilinos::Table<Epetra_SerialDenseVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialDenseVector_ID_t>(
        tab->store(pobj, true));
}

/* store Epetra_SerialDenseVector in non-const table */
CT_Epetra_SerialDenseVector_ID_t
storeSerialDenseVector( Epetra_SerialDenseVector *pobj )
{
    CTrilinos::Table<Epetra_SerialDenseVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialDenseVector_ID_t>(
        tab->store(pobj, false));
}

/* store const Epetra_SerialDenseVector in const table */
CT_Epetra_SerialDenseVector_ID_t
storeConstSerialDenseVector( const Epetra_SerialDenseVector *pobj )
{
    CTrilinos::Table<Epetra_SerialDenseVector> *tab = 0;
    CTrilinos::TableRepos::getTable(tab);
    return CTrilinos::concreteType<CT_Epetra_SerialDenseVector_ID_t>(
        tab->store(pobj, false));
}

/* remove Epetra_SerialDenseVector from table using CT_Epetra_SerialDenseVector_ID */
void
removeSerialDenseVector( CT_Epetra_SerialDenseVector_ID_t *id )
{
    CTrilinos_Universal_ID_t aid = 
        CTrilinos::abstractType<CT_Epetra_SerialDenseVector_ID_t>(*id);
    CTrilinos::TableRepos::remove(&aid);
    *id = CTrilinos::concreteType<CT_Epetra_SerialDenseVector_ID_t>(aid);
}

/* purge Epetra_SerialDenseVector table */
void
purgeSerialDenseVector(  )
{
    CTrilinos::TableRepos::purge<Epetra_SerialDenseVector>();
}

}


