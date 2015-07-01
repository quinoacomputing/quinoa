
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


/*! @file CTrilinos_TableRepos.hpp
 * @brief Central table repository for CTrilinos. */


#ifndef CTRILINOS_TABLEREPOS_HPP
#define CTRILINOS_TABLEREPOS_HPP


#include "CTrilinos_config.h"

#include <string>

#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#ifdef HAVE_MPI
#include "CEpetra_MpiComm_Cpp.hpp"
#endif /* HAVE_MPI */
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#ifdef HAVE_CTRILINOS_AMESOS
#include "CAmesos_BaseSolver_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AMESOS */
#include "CEpetra_FECrsMatrix_Cpp.hpp"
#include "CEpetra_IntSerialDenseVector_Cpp.hpp"
#include "CEpetra_SerialDenseMatrix_Cpp.hpp"
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTest_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTestCombo_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTestMaxIters_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
#include "CAztecOO_StatusTestResNorm_Cpp.hpp"
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
#include "CIfpack_Preconditioner_Cpp.hpp"
#endif /* HAVE_CTRILINOS_IFPACK */
#include "CEpetra_SerialDenseVector_Cpp.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_Table.hpp"
#include "Teuchos_RCP.hpp"


namespace CTrilinos {


namespace TableRepos
{
    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<T> get(CTrilinos_Universal_ID_t id);

    /*! retrieve the object */
    template <class T>
    const Teuchos::RCP<const T> getConst(CTrilinos_Universal_ID_t id);

    /*! remove an object from the table and invalidate the id struct */
    void remove(CTrilinos_Universal_ID_t * id);

    /*! create an alias for the object in another table */
    CTrilinos_Universal_ID_t alias(CTrilinos_Universal_ID_t id, CTrilinos_Table_ID_t tab, bool keepold);

    /*! see if the object is dynamic_cast'able */
    bool typeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type);

}

template <class T>
const Teuchos::RCP<T> TableRepos::get(CTrilinos_Universal_ID_t aid)
{
    switch (aid.table) {
    case CT_Epetra_Distributor_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Distributor>(CEpetra::getDistributor(aid), true);
    case CT_Epetra_SerialComm_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_SerialComm>(CEpetra::getSerialComm(aid), true);
    case CT_Epetra_BLAS_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_BLAS>(CEpetra::getBLAS(aid), true);
    case CT_Epetra_Comm_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Comm>(CEpetra::getComm(aid), true);
    case CT_Epetra_Operator_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Operator>(CEpetra::getOperator(aid), true);
    case CT_Epetra_MultiVector_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_MultiVector>(CEpetra::getMultiVector(aid), true);
    case CT_Epetra_OffsetIndex_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_OffsetIndex>(CEpetra::getOffsetIndex(aid), true);
    case CT_Epetra_Object_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Object>(CEpetra::getObject(aid), true);
    case CT_Epetra_RowMatrix_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_RowMatrix>(CEpetra::getRowMatrix(aid), true);
    case CT_Epetra_CompObject_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_CompObject>(CEpetra::getCompObject(aid), true);
    case CT_Epetra_Directory_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Directory>(CEpetra::getDirectory(aid), true);
    case CT_Epetra_Flops_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Flops>(CEpetra::getFlops(aid), true);
    case CT_Epetra_SrcDistObject_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_SrcDistObject>(CEpetra::getSrcDistObject(aid), true);
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_MpiComm>(CEpetra::getMpiComm(aid), true);
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_CrsMatrix>(CEpetra::getCrsMatrix(aid), true);
    case CT_Epetra_CrsGraph_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_CrsGraph>(CEpetra::getCrsGraph(aid), true);
    case CT_Epetra_DistObject_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_DistObject>(CEpetra::getDistObject(aid), true);
    case CT_Epetra_Vector_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Vector>(CEpetra::getVector(aid), true);
    case CT_Epetra_Export_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Export>(CEpetra::getExport(aid), true);
    case CT_Epetra_Map_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Map>(CEpetra::getMap(aid), true);
    case CT_Epetra_BlockMap_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_BlockMap>(CEpetra::getBlockMap(aid), true);
    case CT_Epetra_Import_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Import>(CEpetra::getImport(aid), true);
    case CT_Epetra_Time_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_Time>(CEpetra::getTime(aid), true);
    case CT_Epetra_JadMatrix_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_JadMatrix>(CEpetra::getJadMatrix(aid), true);
    case CT_Epetra_LinearProblem_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_LinearProblem>(CEpetra::getLinearProblem(aid), true);
    case CT_Epetra_LAPACK_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_LAPACK>(CEpetra::getLAPACK(aid), true);
    case CT_Teuchos_ParameterList_ID:
        return Teuchos::rcp_dynamic_cast<T,Teuchos::ParameterList>(CTeuchos::getParameterList(aid), true);
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        return Teuchos::rcp_dynamic_cast<T,Amesos_BaseSolver>(CAmesos::getBaseSolver(aid), true);
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_FECrsMatrix>(CEpetra::getFECrsMatrix(aid), true);
    case CT_Epetra_IntSerialDenseVector_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_IntSerialDenseVector>(CEpetra::getIntSerialDenseVector(aid), true);
    case CT_Epetra_SerialDenseMatrix_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_SerialDenseMatrix>(CEpetra::getSerialDenseMatrix(aid), true);
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        return Teuchos::rcp_dynamic_cast<T,AztecOO_StatusTest>(CAztecOO::getStatusTest(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        return Teuchos::rcp_dynamic_cast<T,AztecOO_StatusTestCombo>(CAztecOO::getStatusTestCombo(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        return Teuchos::rcp_dynamic_cast<T,AztecOO_StatusTestMaxIters>(CAztecOO::getStatusTestMaxIters(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        return Teuchos::rcp_dynamic_cast<T,AztecOO_StatusTestResNorm>(CAztecOO::getStatusTestResNorm(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        return Teuchos::rcp_dynamic_cast<T,Ifpack_Preconditioner>(CIfpack::getPreconditioner(aid), true);
#endif /* HAVE_CTRILINOS_IFPACK */
    case CT_Epetra_SerialDenseVector_ID:
        return Teuchos::rcp_dynamic_cast<T,Epetra_SerialDenseVector>(CEpetra::getSerialDenseVector(aid), true);
    default:
        throw CTrilinosInvalidTypeError("invalid table id");
    }

    return Teuchos::null;
}

template <class T>
const Teuchos::RCP<const T> TableRepos::getConst(CTrilinos_Universal_ID_t aid)
{
    switch (aid.table) {
    case CT_Epetra_Distributor_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Distributor>(CEpetra::getConstDistributor(aid), true);
    case CT_Epetra_SerialComm_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_SerialComm>(CEpetra::getConstSerialComm(aid), true);
    case CT_Epetra_BLAS_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_BLAS>(CEpetra::getConstBLAS(aid), true);
    case CT_Epetra_Comm_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Comm>(CEpetra::getConstComm(aid), true);
    case CT_Epetra_Operator_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Operator>(CEpetra::getConstOperator(aid), true);
    case CT_Epetra_MultiVector_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_MultiVector>(CEpetra::getConstMultiVector(aid), true);
    case CT_Epetra_OffsetIndex_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_OffsetIndex>(CEpetra::getConstOffsetIndex(aid), true);
    case CT_Epetra_Object_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Object>(CEpetra::getConstObject(aid), true);
    case CT_Epetra_RowMatrix_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_RowMatrix>(CEpetra::getConstRowMatrix(aid), true);
    case CT_Epetra_CompObject_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_CompObject>(CEpetra::getConstCompObject(aid), true);
    case CT_Epetra_Directory_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Directory>(CEpetra::getConstDirectory(aid), true);
    case CT_Epetra_Flops_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Flops>(CEpetra::getConstFlops(aid), true);
    case CT_Epetra_SrcDistObject_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_SrcDistObject>(CEpetra::getConstSrcDistObject(aid), true);
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_MpiComm>(CEpetra::getConstMpiComm(aid), true);
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_CrsMatrix>(CEpetra::getConstCrsMatrix(aid), true);
    case CT_Epetra_CrsGraph_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_CrsGraph>(CEpetra::getConstCrsGraph(aid), true);
    case CT_Epetra_DistObject_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_DistObject>(CEpetra::getConstDistObject(aid), true);
    case CT_Epetra_Vector_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Vector>(CEpetra::getConstVector(aid), true);
    case CT_Epetra_Export_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Export>(CEpetra::getConstExport(aid), true);
    case CT_Epetra_Map_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Map>(CEpetra::getConstMap(aid), true);
    case CT_Epetra_BlockMap_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_BlockMap>(CEpetra::getConstBlockMap(aid), true);
    case CT_Epetra_Import_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Import>(CEpetra::getConstImport(aid), true);
    case CT_Epetra_Time_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_Time>(CEpetra::getConstTime(aid), true);
    case CT_Epetra_JadMatrix_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_JadMatrix>(CEpetra::getConstJadMatrix(aid), true);
    case CT_Epetra_LinearProblem_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_LinearProblem>(CEpetra::getConstLinearProblem(aid), true);
    case CT_Epetra_LAPACK_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_LAPACK>(CEpetra::getConstLAPACK(aid), true);
    case CT_Teuchos_ParameterList_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Teuchos::ParameterList>(CTeuchos::getConstParameterList(aid), true);
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Amesos_BaseSolver>(CAmesos::getConstBaseSolver(aid), true);
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_FECrsMatrix>(CEpetra::getConstFECrsMatrix(aid), true);
    case CT_Epetra_IntSerialDenseVector_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_IntSerialDenseVector>(CEpetra::getConstIntSerialDenseVector(aid), true);
    case CT_Epetra_SerialDenseMatrix_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_SerialDenseMatrix>(CEpetra::getConstSerialDenseMatrix(aid), true);
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        return Teuchos::rcp_dynamic_cast<const T, const AztecOO_StatusTest>(CAztecOO::getConstStatusTest(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        return Teuchos::rcp_dynamic_cast<const T, const AztecOO_StatusTestCombo>(CAztecOO::getConstStatusTestCombo(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        return Teuchos::rcp_dynamic_cast<const T, const AztecOO_StatusTestMaxIters>(CAztecOO::getConstStatusTestMaxIters(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        return Teuchos::rcp_dynamic_cast<const T, const AztecOO_StatusTestResNorm>(CAztecOO::getConstStatusTestResNorm(aid), true);
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Ifpack_Preconditioner>(CIfpack::getConstPreconditioner(aid), true);
#endif /* HAVE_CTRILINOS_IFPACK */
    case CT_Epetra_SerialDenseVector_ID:
        return Teuchos::rcp_dynamic_cast<const T, const Epetra_SerialDenseVector>(CEpetra::getConstSerialDenseVector(aid), true);
    default:
        throw CTrilinosInvalidTypeError("invalid table id");
    }

    return Teuchos::null;
}


} // namespace CTrilinos


#endif

