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
#include "CTrilinos_utils_templ.hpp"
#include "CTrilinos_TableRepos.hpp"
#include "Teuchos_RCP.hpp"

#include "CTrilinos_UnitTestHelpers.hpp"
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
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::RangeError;
using Teuchos::NullReferenceError;
using Teuchos::m_bad_cast;
using CTrilinos::CTrilinosTypeMismatchError;
using CTrilinos::CTrilinosConstCastError;


/* Table::remove() */

TEUCHOS_UNIT_TEST( TableRepos, remove )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CLASS_TYPE(T) id = CEpetra::storeSerialComm(new T));
  ECHO(CTrilinos_Universal_ID_t aid = CTrilinos::abstractType(id));
  TEST_EQUALITY_CONST(aid.index, 0);
  TEST_EQUALITY(aid.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(aid.is_const, FALSE);
  ECHO(CTrilinos::TableRepos::remove(&aid));
  TEST_EQUALITY_CONST(aid.index, -1);
  TEST_EQUALITY(aid.table, CLASS_ENUM(Invalid));
}

TEUCHOS_UNIT_TEST( TableRepos, removeConst )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CLASS_TYPE(T) id = CEpetra::storeConstSerialComm(new T));
  ECHO(CTrilinos_Universal_ID_t aid = CTrilinos::abstractType(id));
  TEST_EQUALITY_CONST(aid.index, 0);
  TEST_EQUALITY(aid.table, CLASS_ENUM(T));
  TEST_EQUALITY_CONST(aid.is_const, TRUE);
  ECHO(CTrilinos::TableRepos::remove(&aid));
  TEST_EQUALITY_CONST(aid.index, -1);
  TEST_EQUALITY(aid.table, CLASS_ENUM(Invalid));
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( TableRepos, removeInvalid )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = -1);
  ECHO(id.table = CLASS_ENUM(T));
  ECHO(id.is_const = FALSE);
  TEST_THROW(CTrilinos::TableRepos::remove(&id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* Table::get() */

TEUCHOS_UNIT_TEST( TableRepos, get )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CLASS_TYPE(T1) id = CEpetra::storeSerialComm(new T1));
  ECHO(CTrilinos_Universal_ID_t aid = CTrilinos::abstractType(id));
  ECHO(RCP<T2> rcpT = CTrilinos::TableRepos::get<T2>(aid));
  TEST_EQUALITY_CONST(nonnull(rcpT), true);
  TEST_EQUALITY_CONST(is_null(rcpT), false);
  ECHO(RCP<const T2> rcpCT = CTrilinos::TableRepos::getConst<T2>(aid));
  TEST_EQUALITY_CONST(nonnull(rcpCT), true);
  TEST_EQUALITY_CONST(is_null(rcpCT), false);
}

#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( TableRepos, getInvalid )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CTrilinos_Universal_ID_t id);
  ECHO(id.index = 0);
  ECHO(id.table = CLASS_ENUM(T));
  ECHO(id.is_const = FALSE);
  TEST_THROW(CTrilinos::TableRepos::get<T>(id), RangeError);
}

#endif /* TEUCHOS_DEBUG */


/* Table::alias() */

TEUCHOS_UNIT_TEST( TableRepos, alias )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CLASS_TYPE(T1) id1 = CEpetra::storeSerialComm(new T1));
  ECHO(CTrilinos_Universal_ID_t aid1 = CTrilinos::abstractType(id1));
  ECHO(CTrilinos_Universal_ID_t aid2 = CTrilinos::TableRepos::alias(aid1, CLASS_ENUM(T2), true));
  TEST_EQUALITY(aid2.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(aid2.index, 0);
  TEST_EQUALITY_CONST(aid2.is_const, FALSE);
}

TEUCHOS_UNIT_TEST( TableRepos, aliasConst )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CLASS_TYPE(T1) id1 = CEpetra::storeConstSerialComm(new const T1));
  ECHO(CTrilinos_Universal_ID_t aid1 = CTrilinos::abstractType(id1));
  ECHO(CTrilinos_Universal_ID_t aid2 = CTrilinos::TableRepos::alias(aid1, CLASS_ENUM(T2), true));
  TEST_EQUALITY(aid2.table, CLASS_ENUM(T2));
  TEST_EQUALITY_CONST(aid2.index, 0);
  TEST_EQUALITY_CONST(aid2.is_const, TRUE);
}

TEUCHOS_UNIT_TEST( TableRepos, aliasBad )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CLASS_TYPE(T3) id = CEpetra::storeSerialComm(new T3));
  ECHO(CTrilinos_Universal_ID_t aid = CTrilinos::abstractType(id));
  TEST_THROW(CTrilinos::TableRepos::alias(aid, CLASS_ENUM(T4), true), m_bad_cast);
}


/* Table::typeCheck() */

TEUCHOS_UNIT_TEST( TableRepos, typeCheck )
{
  ECHO(CEpetra_Test_CleanSlate());
  ECHO(CLASS_TYPE(T1) id = CEpetra::storeSerialComm(new T1));
  ECHO(CTrilinos_Universal_ID_t aid = CTrilinos::abstractType(id));
  TEST_EQUALITY_CONST(CTrilinos::TableRepos::typeCheck(aid, CLASS_ENUM(T1)), true);
  TEST_EQUALITY_CONST(CTrilinos::TableRepos::typeCheck(aid, CLASS_ENUM(T2)), true);
  TEST_EQUALITY_CONST(CTrilinos::TableRepos::typeCheck(aid, CLASS_ENUM(T4)), false);
}

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
