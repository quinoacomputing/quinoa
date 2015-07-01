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
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"
#include "Teuchos_RCP.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "CEpetra_SerialComm.h"
#include "CEpetra_Comm.h"
#include "CEpetra_Vector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"


#define JOIN_SET_0(A, B, C) A ## B ## C
#define JOIN_SET(A, B, C)   JOIN_SET_0(A, B, C)

#define BUILD_CALL(A, F) JOIN_SET( A , _ , F )
#define CLASS_TYPE(A)    JOIN_SET( CT_ , A , _ID_t )
#define CLASS_ENUM(A)    JOIN_SET( CT_ , A , _ID )
#define CLASS_ESTR(A)    XSTRFY(CLASS_ENUM(A))
#define STRFY(A)         #A
#define XSTRFY(A)        STRFY(A)
#define CONSTRUCTOR(A)   A


#define T Epetra_SerialComm
#define T1 Epetra_SerialComm
#define T2 Epetra_Comm
#define T3 Epetra_SerialComm
#define T4 Epetra_Vector


namespace {


using Teuchos::null;
using Teuchos::nonnull;
using Teuchos::is_null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::RangeError;
using Teuchos::NullReferenceError;
using Teuchos::m_bad_cast;
using CTrilinos::CTrilinosTypeMismatchError;
using CTrilinos::CTrilinosConstCastError;
using CTrilinos::CTrilinosWrongTableError;
using CTrilinos::Table;


/* Table::store() owned */

TEUCHOS_UNIT_TEST( Table, store )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id = table.store(new T, true));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY_CONST(id.is_const, FALSE);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(nonnull(table.get<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.get<T>(id)), false);
  TEST_EQUALITY_CONST(nonnull(table.getConst<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.getConst<T>(id)), false);
}

TEUCHOS_UNIT_TEST( Table, storeBase )
{
  ECHO(Table<T2> table(CONSTRUCTOR(CLASS_ENUM(T2))));
  ECHO(CTrilinos_Universal_ID_t id = table.store(new T1, true));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY_CONST(id.is_const, FALSE);
  TEST_EQUALITY(id.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(nonnull(table.get<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.get<T2>(id)), false);
  TEST_EQUALITY_CONST(nonnull(table.getConst<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.getConst<T2>(id)), false);
}

/* this should not compile -- ok */
/*
TEUCHOS_UNIT_TEST( Table, storeWrong )
{
  ECHO(Table<T4> table(CONSTRUCTOR(CLASS_ENUM(T4))));
  TEST_THROW(table.store(new T3, true), CTrilinosTypeMismatchError);
}
*/

TEUCHOS_UNIT_TEST( Table, storeNull )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(T* pobj = NULL);
  TEST_THROW(table.store(pobj, false), NullReferenceError); 
}


/* Table::store() non-owned */

TEUCHOS_UNIT_TEST( Table, storeShared )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(T *pobj = new T);
  ECHO(CTrilinos_Universal_ID_t id = table.store(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY(id.is_const, FALSE);
  TEST_EQUALITY_CONST(nonnull(table.get<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.get<T>(id)), false);
  TEST_EQUALITY_CONST(nonnull(table.getConst<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.getConst<T>(id)), false);
  ECHO(table.remove(&id));
  ECHO(delete pobj);
}

TEUCHOS_UNIT_TEST( Table, storeConstShared )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(const T *pobj = new T);
  ECHO(CTrilinos_Universal_ID_t id = table.store(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY(id.is_const, TRUE);
  TEST_EQUALITY_CONST(nonnull(table.getConst<T>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.getConst<T>(id)), false);
  TEST_THROW(nonnull(table.get<T>(id)), CTrilinosConstCastError);
  TEST_THROW(is_null(table.get<T>(id)), CTrilinosConstCastError);
  ECHO(table.remove(&id));
  ECHO(delete pobj);
}

TEUCHOS_UNIT_TEST( Table, storeSharedBase )
{
  ECHO(Table<T2> table(CONSTRUCTOR(CLASS_ENUM(T2))));
  ECHO(T1 *pobj = new T1);
  ECHO(CTrilinos_Universal_ID_t id = table.store(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T2));
  TEST_EQUALITY(id.is_const, FALSE);
  TEST_EQUALITY_CONST(nonnull(table.get<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.get<T2>(id)), false);
  TEST_EQUALITY_CONST(nonnull(table.getConst<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.getConst<T2>(id)), false);
  ECHO(table.remove(&id));
  ECHO(delete pobj);
}

TEUCHOS_UNIT_TEST( Table, storeConstSharedBase )
{
  ECHO(Table<T2> table(CONSTRUCTOR(CLASS_ENUM(T2))));
  ECHO(const T1 *pobj = new T1);
  ECHO(CTrilinos_Universal_ID_t id = table.store(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T2));
  TEST_EQUALITY(id.is_const, TRUE);
  TEST_EQUALITY_CONST(nonnull(table.getConst<T2>(id)), true);
  TEST_EQUALITY_CONST(is_null(table.getConst<T2>(id)), false);
  TEST_THROW(nonnull(table.get<T2>(id)), CTrilinosConstCastError);
  TEST_THROW(is_null(table.get<T2>(id)), CTrilinosConstCastError);
  ECHO(table.remove(&id));
  ECHO(delete pobj);
}

/* this should not compile -- ok */
/*
TEUCHOS_UNIT_TEST( Table, storeSharedWrong )
{
  ECHO(Table<T4> table(CONSTRUCTOR(CLASS_ENUM(T4))));
  ECHO(T3 *pobj = new T3);
  TEST_THROW(table.store(pobj, false), CTrilinosTypeMismatchError);
  ECHO(delete pobj);
}
*/

TEUCHOS_UNIT_TEST( Table, storeSharedNull )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(T* pobj = NULL);
  TEST_THROW(table.store(pobj, false), NullReferenceError); 
}

TEUCHOS_UNIT_TEST( Table, storeConstSharedNull )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(const T* pobj = NULL);
  TEST_THROW(table.store(pobj, false), NullReferenceError); 
}


/* Table::remove() */

TEUCHOS_UNIT_TEST( Table, remove )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id = table.store(new T, true));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(id.is_const, FALSE);
  ECHO(table.remove(&id));
  TEST_EQUALITY_CONST(id.index, -1);
  TEST_EQUALITY(id.table, CLASS_ENUM(Invalid));
}

TEUCHOS_UNIT_TEST( Table, removeConst )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(const T* pobj = new T);
  ECHO(CTrilinos_Universal_ID_t id = table.store(pobj, false));
  TEST_EQUALITY_CONST(id.index, 0);
  TEST_EQUALITY(id.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(id.is_const, TRUE);
  ECHO(table.remove(&id));
  TEST_EQUALITY_CONST(id.index, -1);
  TEST_EQUALITY(id.table, CLASS_ENUM(Invalid));
  ECHO(delete pobj);
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Table, removeInvalid )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = -1);
  ECHO(id.table = CLASS_ENUM(T));
  ECHO(id.is_const = FALSE);
  TEST_THROW(table.remove(&id), RangeError);
}

#endif /* TEUCHOS_DEBUG */

TEUCHOS_UNIT_TEST( Table, removeWrong )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id1 = table.store(new T, true));
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = id1.index);
  ECHO(id.table = CLASS_ENUM(T4));
  ECHO(id.is_const = FALSE);
  TEST_THROW(table.remove(&id), CTrilinosWrongTableError);
}


/* Table::get() */

TEUCHOS_UNIT_TEST( Table, get )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id = table.store(new T, true));
  ECHO(RCP<T> rcpT = table.get<T>(id));
  TEST_EQUALITY_CONST(nonnull(rcpT), true);
  TEST_EQUALITY_CONST(is_null(rcpT), false);
  ECHO(RCP<const T> rcpCT = table.getConst<T>(id));
  TEST_EQUALITY_CONST(nonnull(rcpCT), true);
  TEST_EQUALITY_CONST(is_null(rcpCT), false);
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Table, getInvalid )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = 0);
  ECHO(id.table = CLASS_ENUM(T));
  ECHO(id.is_const = FALSE);
  TEST_THROW(table.get<T>(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */

TEUCHOS_UNIT_TEST( Table, getWrong )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id1 = table.store(new T, true));
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = id1.index);
  ECHO(id.table = CLASS_ENUM(T4));
  ECHO(id.is_const = FALSE);
  TEST_THROW(table.get<T>(id), CTrilinosWrongTableError);
}


/* Table::alias() */

TEUCHOS_UNIT_TEST( Table, alias )
{
  ECHO(Table<T1> table1(CONSTRUCTOR(CLASS_ENUM(T1))));
  ECHO(Table<T2> table2(CONSTRUCTOR(CLASS_ENUM(T2))));

  ECHO(CTrilinos_Universal_ID_t id1 = table1.store(new T1, true));
  ECHO(CTrilinos_Universal_ID_t id2 = table2.alias(table1.get<T2>(id1)));

  TEST_EQUALITY(id2.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(id2.index, 0);
}

TEUCHOS_UNIT_TEST( Table, aliasConst )
{
  ECHO(Table<T1> table1(CONSTRUCTOR(CLASS_ENUM(T1))));
  ECHO(Table<T2> table2(CONSTRUCTOR(CLASS_ENUM(T2))));

  ECHO(CTrilinos_Universal_ID_t id1 = table1.store(new const T1, true));
  ECHO(CTrilinos_Universal_ID_t id2 = table2.alias(table1.getConst<T2>(id1)));

  TEST_EQUALITY(id2.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(id2.index, 0);
  TEST_EQUALITY_CONST(id2.is_const, TRUE);
}

TEUCHOS_UNIT_TEST( Table, aliasBad )
{
  ECHO(Table<T3> table3(CONSTRUCTOR(CLASS_ENUM(T3))));
  ECHO(Table<T4> table4(CONSTRUCTOR(CLASS_ENUM(T4))));

  ECHO(CTrilinos_Universal_ID_t id3 = table3.store(new T3, true));
  TEST_THROW(table4.alias(table3.get<T4>(id3)), m_bad_cast);
}


/* Table::purge() */

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Table, purge )
{
  ECHO(Table<T> table(CONSTRUCTOR(CLASS_ENUM(T))));
  ECHO(CTrilinos_Universal_ID_t id = table.store(new T, true));
  TEST_EQUALITY_CONST(nonnull(table.get<T>(id)), true);
  ECHO(table.purge());
  TEST_THROW(table.get<T>(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* Table::isType() */

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Table, isType )
{
  ECHO(Table<T3> table3(CONSTRUCTOR(CLASS_ENUM(T3))));
  ECHO(Table<T4> table4(CONSTRUCTOR(CLASS_ENUM(T4))));
  ECHO(CTrilinos_Universal_ID_t id = table3.store(new T3, true));
  TEST_EQUALITY_CONST(table3.isType(id.table), true);
  TEST_EQUALITY_CONST(table4.isType(id.table), false);
}

#endif /* TEUCHOS_DEBUG */


//
// Template Instantiations
//


#ifdef TEUCHOS_DEBUG

#  define DEBUG_UNIT_TEST_GROUP( TT ) \

#else

#  define DEBUG_UNIT_TEST_GROUP( TT )

#endif


#define UNIT_TEST_GROUP( TT ) \
  DEBUG_UNIT_TEST_GROUP( TT )


} // namespace
