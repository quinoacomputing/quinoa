
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


/*! @file CTrilinos_enums.h
 * @brief Defines structs and enums needed for CTrilinos. */


#ifndef CTRILINOS_ENUMS_H
#define CTRILINOS_ENUMS_H


#include "CTrilinos_config.h"
#include "CTrilinos_enum_wrappers.h"


#ifdef __cplusplus
extern "C" {
#endif


/*! C does not support the C++ bool type, so the enum below will
 * act as a custom boolean type for users using a C compiler.  #TRUE is
 * the equivalent of true and #FALSE is the equivalent of false. */

typedef int boolean;

/*! #boolean uses TRUE instead of true. */
#ifndef TRUE
#  define TRUE 1
#endif

/*! #boolean uses FALSE instead of false. */
#ifndef FALSE
#  define FALSE 0
#endif

/*! The enum below lists all the classes that CTrilinos supports.  Classes
 * that are derived from a base class listed below can also be used with some
 * restrictions: (1) the only methods that can be called on such classes are
 * the ones defined by the base class and (2) objects of the derived type can
 * only be passed to methods expecting an argument of the base-class type. */

typedef enum {
    CT_Invalid_ID,                       /*!< does not refer to a valid table entry */
    CT_Epetra_Distributor_ID,            /*!< refers to an Epetra_Distributor table entry */
    CT_Epetra_SerialComm_ID,             /*!< refers to an Epetra_SerialComm table entry */
    CT_Epetra_BLAS_ID,                   /*!< refers to an Epetra_BLAS table entry */
    CT_Epetra_Comm_ID,                   /*!< refers to an Epetra_Comm table entry */
    CT_Epetra_Operator_ID,               /*!< refers to an Epetra_Operator table entry */
    CT_Epetra_MultiVector_ID,            /*!< refers to an Epetra_MultiVector table entry */
    CT_Epetra_OffsetIndex_ID,            /*!< refers to an Epetra_OffsetIndex table entry */
    CT_Epetra_Object_ID,                 /*!< refers to an Epetra_Object table entry */
    CT_Epetra_RowMatrix_ID,              /*!< refers to an Epetra_RowMatrix table entry */
    CT_Epetra_CompObject_ID,             /*!< refers to an Epetra_CompObject table entry */
    CT_Epetra_Directory_ID,              /*!< refers to an Epetra_Directory table entry */
    CT_Epetra_Flops_ID,                  /*!< refers to an Epetra_Flops table entry */
    CT_Epetra_SrcDistObject_ID,          /*!< refers to an Epetra_SrcDistObject table entry */
    CT_Epetra_MpiComm_ID,                /*!< refers to an Epetra_MpiComm table entry */
    CT_Epetra_CrsMatrix_ID,              /*!< refers to an Epetra_CrsMatrix table entry */
    CT_Epetra_CrsGraph_ID,               /*!< refers to an Epetra_CrsGraph table entry */
    CT_Epetra_DistObject_ID,             /*!< refers to an Epetra_DistObject table entry */
    CT_Epetra_Vector_ID,                 /*!< refers to an Epetra_Vector table entry */
    CT_Epetra_Export_ID,                 /*!< refers to an Epetra_Export table entry */
    CT_Epetra_Map_ID,                    /*!< refers to an Epetra_Map table entry */
    CT_Epetra_BlockMap_ID,               /*!< refers to an Epetra_BlockMap table entry */
    CT_Epetra_Import_ID,                 /*!< refers to an Epetra_Import table entry */
    CT_Epetra_Time_ID,                   /*!< refers to an Epetra_Time table entry */
    CT_Epetra_JadMatrix_ID,              /*!< refers to an Epetra_JadMatrix table entry */
    CT_Epetra_LinearProblem_ID,          /*!< refers to an Epetra_LinearProblem table entry */
    CT_Epetra_LAPACK_ID,                 /*!< refers to an Epetra_LAPACK table entry */
    CT_Teuchos_CommandLineProcessor_ID,  /*!< refers to a Teuchos::CommandLineProcessor table entry */
    CT_Teuchos_ParameterList_ID,         /*!< refers to a Teuchos::ParameterList table entry */
    CT_Teuchos_ParameterEntry_ID,        /*!< refers to a Teuchos::ParameterEntry table entry */
    CT_Teuchos_any_ID,                   /*!< refers to a Teuchos::any table entry */
    CT_Amesos_BaseSolver_ID,             /*!< refers to an Amesos_BaseSolver table entry */
    CT_Amesos_ID,                        /*!< refers to an Amesos table entry */
    CT_Epetra_FECrsMatrix_ID,            /*!< refers to an Epetra_FECrsMatrix table entry */
    CT_Epetra_IntSerialDenseVector_ID,   /*!< refers to an Epetra_IntSerialDenseVector table entry */
    CT_Epetra_SerialDenseMatrix_ID,      /*!< refers to an Epetra_SerialDenseMatrix table entry */
    CT_AztecOO_ID,                       /*!< refers to an AztecOO table entry */
    CT_AztecOO_StatusTest_ID,            /*!< refers to an AztecOO_StatusTest table entry */
    CT_AztecOO_StatusTestCombo_ID,       /*!< refers to an AztecOO_StatusTestCombo table entry */
    CT_AztecOO_StatusTestMaxIters_ID,    /*!< refers to an AztecOO_StatusTestMaxIters table entry */
    CT_AztecOO_StatusTestResNorm_ID,     /*!< refers to an AztecOO_StatusTestResNorm table entry */
    CT_Ifpack_ID,                        /*!< refers to an Ifpack table entry */
    CT_Ifpack_Preconditioner_ID,         /*!< refers to an Ifpack_Preconditioner table entry */
    CT_Epetra_SerialDenseVector_ID,      /*!< refers to an Epetra_SerialDenseVector table entry */
    CT_Pliris_ID                         /*!< refers to a Pliris table entry */
} CTrilinos_Table_ID_t;

/*! The type in the struct below is actually used to identify the table in
 * which the object is stored and not the type of the underlying object.  The
 * end user should NEVER change the content of this struct as it will result
 * in segfaults and/or memory leaks.  When the user asks that the object
 * be removed from the table, the type will be changed to CT_Invalid_ID by the
 * destructor function. */

typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CTrilinos_Universal_ID_t;

/* All structs for specific types are identical but are declared separately
 * rather than as a typedef of the generic one so that they are not inter-
 * changeable.  This forces the compiler to perform type-checking on calls
 * to the wrapper functions and will make debugging easier for the user. */

/*! Struct used for referring to objects in the Epetra_Distributor table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Distributor.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Distributor_ID_t;

/*! Struct used for referring to objects in the Epetra_SerialComm table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_SerialComm.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_SerialComm_ID_t;

/*! Struct used for referring to objects in the Epetra_BLAS table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_BLAS.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_BLAS_ID_t;

/*! Struct used for referring to objects in the Epetra_Comm table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Comm.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Comm_ID_t;

/*! Struct used for referring to objects in the Epetra_Operator table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Operator.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Operator_ID_t;

/*! Struct used for referring to objects in the Epetra_MultiVector table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_MultiVector.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_MultiVector_ID_t;

/*! Struct used for referring to objects in the Epetra_OffsetIndex table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_OffsetIndex.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_OffsetIndex_ID_t;

/*! Struct used for referring to objects in the Epetra_Object table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Object.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Object_ID_t;

/*! Struct used for referring to objects in the Epetra_RowMatrix table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_RowMatrix.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_RowMatrix_ID_t;

/*! Struct used for referring to objects in the Epetra_CompObject table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_CompObject.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_CompObject_ID_t;

/*! Struct used for referring to objects in the Epetra_Directory table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Directory.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Directory_ID_t;

/*! Struct used for referring to objects in the Epetra_Flops table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Flops.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Flops_ID_t;

/*! Struct used for referring to objects in the Epetra_SrcDistObject table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_SrcDistObject.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_SrcDistObject_ID_t;

#ifdef HAVE_MPI
/*! Struct used for referring to objects in the Epetra_MpiComm table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_MpiComm.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_MpiComm_ID_t;
#endif /* HAVE_MPI */

/*! Struct used for referring to objects in the Epetra_CrsMatrix table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_CrsMatrix.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_CrsMatrix_ID_t;

/*! Struct used for referring to objects in the Epetra_CrsGraph table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_CrsGraph.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_CrsGraph_ID_t;

/*! Struct used for referring to objects in the Epetra_DistObject table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_DistObject.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_DistObject_ID_t;

/*! Struct used for referring to objects in the Epetra_Vector table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Vector.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Vector_ID_t;

/*! Struct used for referring to objects in the Epetra_Export table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Export.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Export_ID_t;

/*! Struct used for referring to objects in the Epetra_Map table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Map.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Map_ID_t;

/*! Struct used for referring to objects in the Epetra_BlockMap table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_BlockMap.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_BlockMap_ID_t;

/*! Struct used for referring to objects in the Epetra_Import table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Import.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Import_ID_t;

/*! Struct used for referring to objects in the Epetra_Time table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_Time.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_Time_ID_t;

/*! Struct used for referring to objects in the Epetra_JadMatrix table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_JadMatrix.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_JadMatrix_ID_t;

/*! Struct used for referring to objects in the Epetra_LinearProblem table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_LinearProblem.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_LinearProblem_ID_t;

/*! Struct used for referring to objects in the Epetra_LAPACK table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_LAPACK.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_LAPACK_ID_t;

/*! Struct used for referring to objects in the Teuchos::CommandLineProcessor table.  Methods
 * that can be invoked on the underlying objects are listed in CTeuchos_CommandLineProcessor.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Teuchos_CommandLineProcessor_ID_t;

/*! Struct used for referring to objects in the Teuchos::ParameterList table.  Methods
 * that can be invoked on the underlying objects are listed in CTeuchos_ParameterList.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Teuchos_ParameterList_ID_t;

/*! Struct used for referring to objects in the Teuchos::ParameterEntry table.  Methods
 * that can be invoked on the underlying objects are listed in CTeuchos_ParameterEntry.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Teuchos_ParameterEntry_ID_t;

/*! Struct used for referring to objects in the Teuchos::any table.  Methods
 * that can be invoked on the underlying objects are listed in CTeuchos_any.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Teuchos_any_ID_t;

#ifdef HAVE_CTRILINOS_AMESOS
/*! Struct used for referring to objects in the Amesos_BaseSolver table.  Methods
 * that can be invoked on the underlying objects are listed in CAmesos_BaseSolver.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Amesos_BaseSolver_ID_t;
#endif /* HAVE_CTRILINOS_AMESOS */

#ifdef HAVE_CTRILINOS_AMESOS
/*! Struct used for referring to objects in the Amesos table.  Methods
 * that can be invoked on the underlying objects are listed in CAmesos.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Amesos_ID_t;
#endif /* HAVE_CTRILINOS_AMESOS */

/*! Struct used for referring to objects in the Epetra_FECrsMatrix table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_FECrsMatrix.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_FECrsMatrix_ID_t;

/*! Struct used for referring to objects in the Epetra_IntSerialDenseVector table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_IntSerialDenseVector.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_IntSerialDenseVector_ID_t;

/*! Struct used for referring to objects in the Epetra_SerialDenseMatrix table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_SerialDenseMatrix.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_SerialDenseMatrix_ID_t;

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Struct used for referring to objects in the AztecOO table.  Methods
 * that can be invoked on the underlying objects are listed in CAztecOO.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_AztecOO_ID_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Struct used for referring to objects in the AztecOO_StatusTest table.  Methods
 * that can be invoked on the underlying objects are listed in CAztecOO_StatusTest.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_AztecOO_StatusTest_ID_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Struct used for referring to objects in the AztecOO_StatusTestCombo table.  Methods
 * that can be invoked on the underlying objects are listed in CAztecOO_StatusTestCombo.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_AztecOO_StatusTestCombo_ID_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Struct used for referring to objects in the AztecOO_StatusTestMaxIters table.  Methods
 * that can be invoked on the underlying objects are listed in CAztecOO_StatusTestMaxIters.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_AztecOO_StatusTestMaxIters_ID_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_AZTECOO
/*! Struct used for referring to objects in the AztecOO_StatusTestResNorm table.  Methods
 * that can be invoked on the underlying objects are listed in CAztecOO_StatusTestResNorm.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_AztecOO_StatusTestResNorm_ID_t;
#endif /* HAVE_CTRILINOS_AZTECOO */

#ifdef HAVE_CTRILINOS_IFPACK
/*! Struct used for referring to objects in the Ifpack table.  Methods
 * that can be invoked on the underlying objects are listed in CIfpack.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Ifpack_ID_t;
#endif /* HAVE_CTRILINOS_IFPACK */

#ifdef HAVE_CTRILINOS_IFPACK
/*! Struct used for referring to objects in the Ifpack_Preconditioner table.  Methods
 * that can be invoked on the underlying objects are listed in CIfpack_Preconditioner.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Ifpack_Preconditioner_ID_t;
#endif /* HAVE_CTRILINOS_IFPACK */

/*! Struct used for referring to objects in the Epetra_SerialDenseVector table.  Methods
 * that can be invoked on the underlying objects are listed in CEpetra_SerialDenseVector.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Epetra_SerialDenseVector_ID_t;

#ifdef HAVE_CTRILINOS_PLIRIS
#ifdef HAVE_MPI
/*! Struct used for referring to objects in the Pliris table.  Methods
 * that can be invoked on the underlying objects are listed in CPliris.h */
typedef struct {
    CTrilinos_Table_ID_t table;	/*!< Table holding reference to the object */
    int index;			/*!< Array index of the object */
    boolean is_const;		/*!< Whether or not object was declared const */
} CT_Pliris_ID_t;
#endif /* HAVE_MPI */
#endif /* HAVE_CTRILINOS_PLIRIS */


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
