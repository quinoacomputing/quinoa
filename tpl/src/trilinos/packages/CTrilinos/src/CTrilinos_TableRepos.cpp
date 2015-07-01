
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


/*! @file CTrilinos_TableRepos.cpp
 * @brief Central table repository for CTrilinos. */


#include "CTrilinos_config.h"
#include "CTrilinos_TableRepos.hpp"



namespace CTrilinos {


void TableRepos::remove(
    CTrilinos_Universal_ID_t * aid)
{
    switch (aid->table) {
    case CT_Epetra_Distributor_ID:
        CEpetra::removeDistributor(aid);
        break;
    case CT_Epetra_SerialComm_ID:
        CEpetra::removeSerialComm(aid);
        break;
    case CT_Epetra_BLAS_ID:
        CEpetra::removeBLAS(aid);
        break;
    case CT_Epetra_Comm_ID:
        CEpetra::removeComm(aid);
        break;
    case CT_Epetra_Operator_ID:
        CEpetra::removeOperator(aid);
        break;
    case CT_Epetra_MultiVector_ID:
        CEpetra::removeMultiVector(aid);
        break;
    case CT_Epetra_OffsetIndex_ID:
        CEpetra::removeOffsetIndex(aid);
        break;
    case CT_Epetra_Object_ID:
        CEpetra::removeObject(aid);
        break;
    case CT_Epetra_RowMatrix_ID:
        CEpetra::removeRowMatrix(aid);
        break;
    case CT_Epetra_CompObject_ID:
        CEpetra::removeCompObject(aid);
        break;
    case CT_Epetra_Directory_ID:
        CEpetra::removeDirectory(aid);
        break;
    case CT_Epetra_Flops_ID:
        CEpetra::removeFlops(aid);
        break;
    case CT_Epetra_SrcDistObject_ID:
        CEpetra::removeSrcDistObject(aid);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        CEpetra::removeMpiComm(aid);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        CEpetra::removeCrsMatrix(aid);
        break;
    case CT_Epetra_CrsGraph_ID:
        CEpetra::removeCrsGraph(aid);
        break;
    case CT_Epetra_DistObject_ID:
        CEpetra::removeDistObject(aid);
        break;
    case CT_Epetra_Vector_ID:
        CEpetra::removeVector(aid);
        break;
    case CT_Epetra_Export_ID:
        CEpetra::removeExport(aid);
        break;
    case CT_Epetra_Map_ID:
        CEpetra::removeMap(aid);
        break;
    case CT_Epetra_BlockMap_ID:
        CEpetra::removeBlockMap(aid);
        break;
    case CT_Epetra_Import_ID:
        CEpetra::removeImport(aid);
        break;
    case CT_Epetra_Time_ID:
        CEpetra::removeTime(aid);
        break;
    case CT_Epetra_JadMatrix_ID:
        CEpetra::removeJadMatrix(aid);
        break;
    case CT_Epetra_LinearProblem_ID:
        CEpetra::removeLinearProblem(aid);
        break;
    case CT_Epetra_LAPACK_ID:
        CEpetra::removeLAPACK(aid);
        break;
    case CT_Teuchos_ParameterList_ID:
        CTeuchos::removeParameterList(aid);
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        CAmesos::removeBaseSolver(aid);
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        CEpetra::removeFECrsMatrix(aid);
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        CEpetra::removeIntSerialDenseVector(aid);
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        CEpetra::removeSerialDenseMatrix(aid);
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        CAztecOO::removeStatusTest(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        CAztecOO::removeStatusTestCombo(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        CAztecOO::removeStatusTestMaxIters(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        CAztecOO::removeStatusTestResNorm(aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        CIfpack::removePreconditioner(aid);
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    case CT_Epetra_SerialDenseVector_ID:
        CEpetra::removeSerialDenseVector(aid);
        break;
    default:
        throw CTrilinosInvalidTypeError("invalid table id or non-polymorphic class");
    }
}

CTrilinos_Universal_ID_t TableRepos::alias(
    CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t tab, bool keepold)
{
    CTrilinos_Universal_ID_t newid;

    switch (tab) {
    case CT_Epetra_Distributor_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Distributor> robj = getConst<Epetra_Distributor>(aid);
            newid = CEpetra::aliasConstDistributor(robj);
        } else {
            Teuchos::RCP<Epetra_Distributor> robj = get<Epetra_Distributor>(aid);
            newid = CEpetra::aliasDistributor(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_SerialComm_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_SerialComm> robj = getConst<Epetra_SerialComm>(aid);
            newid = CEpetra::aliasConstSerialComm(robj);
        } else {
            Teuchos::RCP<Epetra_SerialComm> robj = get<Epetra_SerialComm>(aid);
            newid = CEpetra::aliasSerialComm(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_BLAS_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_BLAS> robj = getConst<Epetra_BLAS>(aid);
            newid = CEpetra::aliasConstBLAS(robj);
        } else {
            Teuchos::RCP<Epetra_BLAS> robj = get<Epetra_BLAS>(aid);
            newid = CEpetra::aliasBLAS(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Comm_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Comm> robj = getConst<Epetra_Comm>(aid);
            newid = CEpetra::aliasConstComm(robj);
        } else {
            Teuchos::RCP<Epetra_Comm> robj = get<Epetra_Comm>(aid);
            newid = CEpetra::aliasComm(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Operator_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Operator> robj = getConst<Epetra_Operator>(aid);
            newid = CEpetra::aliasConstOperator(robj);
        } else {
            Teuchos::RCP<Epetra_Operator> robj = get<Epetra_Operator>(aid);
            newid = CEpetra::aliasOperator(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_MultiVector_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_MultiVector> robj = getConst<Epetra_MultiVector>(aid);
            newid = CEpetra::aliasConstMultiVector(robj);
        } else {
            Teuchos::RCP<Epetra_MultiVector> robj = get<Epetra_MultiVector>(aid);
            newid = CEpetra::aliasMultiVector(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_OffsetIndex_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_OffsetIndex> robj = getConst<Epetra_OffsetIndex>(aid);
            newid = CEpetra::aliasConstOffsetIndex(robj);
        } else {
            Teuchos::RCP<Epetra_OffsetIndex> robj = get<Epetra_OffsetIndex>(aid);
            newid = CEpetra::aliasOffsetIndex(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Object_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Object> robj = getConst<Epetra_Object>(aid);
            newid = CEpetra::aliasConstObject(robj);
        } else {
            Teuchos::RCP<Epetra_Object> robj = get<Epetra_Object>(aid);
            newid = CEpetra::aliasObject(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_RowMatrix_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_RowMatrix> robj = getConst<Epetra_RowMatrix>(aid);
            newid = CEpetra::aliasConstRowMatrix(robj);
        } else {
            Teuchos::RCP<Epetra_RowMatrix> robj = get<Epetra_RowMatrix>(aid);
            newid = CEpetra::aliasRowMatrix(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_CompObject_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_CompObject> robj = getConst<Epetra_CompObject>(aid);
            newid = CEpetra::aliasConstCompObject(robj);
        } else {
            Teuchos::RCP<Epetra_CompObject> robj = get<Epetra_CompObject>(aid);
            newid = CEpetra::aliasCompObject(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Directory_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Directory> robj = getConst<Epetra_Directory>(aid);
            newid = CEpetra::aliasConstDirectory(robj);
        } else {
            Teuchos::RCP<Epetra_Directory> robj = get<Epetra_Directory>(aid);
            newid = CEpetra::aliasDirectory(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Flops_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Flops> robj = getConst<Epetra_Flops>(aid);
            newid = CEpetra::aliasConstFlops(robj);
        } else {
            Teuchos::RCP<Epetra_Flops> robj = get<Epetra_Flops>(aid);
            newid = CEpetra::aliasFlops(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_SrcDistObject_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_SrcDistObject> robj = getConst<Epetra_SrcDistObject>(aid);
            newid = CEpetra::aliasConstSrcDistObject(robj);
        } else {
            Teuchos::RCP<Epetra_SrcDistObject> robj = get<Epetra_SrcDistObject>(aid);
            newid = CEpetra::aliasSrcDistObject(robj);
        }
        if (!keepold) remove(&aid);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_MpiComm> robj = getConst<Epetra_MpiComm>(aid);
            newid = CEpetra::aliasConstMpiComm(robj);
        } else {
            Teuchos::RCP<Epetra_MpiComm> robj = get<Epetra_MpiComm>(aid);
            newid = CEpetra::aliasMpiComm(robj);
        }
        if (!keepold) remove(&aid);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_CrsMatrix> robj = getConst<Epetra_CrsMatrix>(aid);
            newid = CEpetra::aliasConstCrsMatrix(robj);
        } else {
            Teuchos::RCP<Epetra_CrsMatrix> robj = get<Epetra_CrsMatrix>(aid);
            newid = CEpetra::aliasCrsMatrix(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_CrsGraph_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_CrsGraph> robj = getConst<Epetra_CrsGraph>(aid);
            newid = CEpetra::aliasConstCrsGraph(robj);
        } else {
            Teuchos::RCP<Epetra_CrsGraph> robj = get<Epetra_CrsGraph>(aid);
            newid = CEpetra::aliasCrsGraph(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_DistObject_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_DistObject> robj = getConst<Epetra_DistObject>(aid);
            newid = CEpetra::aliasConstDistObject(robj);
        } else {
            Teuchos::RCP<Epetra_DistObject> robj = get<Epetra_DistObject>(aid);
            newid = CEpetra::aliasDistObject(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Vector_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Vector> robj = getConst<Epetra_Vector>(aid);
            newid = CEpetra::aliasConstVector(robj);
        } else {
            Teuchos::RCP<Epetra_Vector> robj = get<Epetra_Vector>(aid);
            newid = CEpetra::aliasVector(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Export_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Export> robj = getConst<Epetra_Export>(aid);
            newid = CEpetra::aliasConstExport(robj);
        } else {
            Teuchos::RCP<Epetra_Export> robj = get<Epetra_Export>(aid);
            newid = CEpetra::aliasExport(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Map_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Map> robj = getConst<Epetra_Map>(aid);
            newid = CEpetra::aliasConstMap(robj);
        } else {
            Teuchos::RCP<Epetra_Map> robj = get<Epetra_Map>(aid);
            newid = CEpetra::aliasMap(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_BlockMap_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_BlockMap> robj = getConst<Epetra_BlockMap>(aid);
            newid = CEpetra::aliasConstBlockMap(robj);
        } else {
            Teuchos::RCP<Epetra_BlockMap> robj = get<Epetra_BlockMap>(aid);
            newid = CEpetra::aliasBlockMap(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Import_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Import> robj = getConst<Epetra_Import>(aid);
            newid = CEpetra::aliasConstImport(robj);
        } else {
            Teuchos::RCP<Epetra_Import> robj = get<Epetra_Import>(aid);
            newid = CEpetra::aliasImport(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_Time_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_Time> robj = getConst<Epetra_Time>(aid);
            newid = CEpetra::aliasConstTime(robj);
        } else {
            Teuchos::RCP<Epetra_Time> robj = get<Epetra_Time>(aid);
            newid = CEpetra::aliasTime(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_JadMatrix_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_JadMatrix> robj = getConst<Epetra_JadMatrix>(aid);
            newid = CEpetra::aliasConstJadMatrix(robj);
        } else {
            Teuchos::RCP<Epetra_JadMatrix> robj = get<Epetra_JadMatrix>(aid);
            newid = CEpetra::aliasJadMatrix(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_LinearProblem_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_LinearProblem> robj = getConst<Epetra_LinearProblem>(aid);
            newid = CEpetra::aliasConstLinearProblem(robj);
        } else {
            Teuchos::RCP<Epetra_LinearProblem> robj = get<Epetra_LinearProblem>(aid);
            newid = CEpetra::aliasLinearProblem(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_LAPACK_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_LAPACK> robj = getConst<Epetra_LAPACK>(aid);
            newid = CEpetra::aliasConstLAPACK(robj);
        } else {
            Teuchos::RCP<Epetra_LAPACK> robj = get<Epetra_LAPACK>(aid);
            newid = CEpetra::aliasLAPACK(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Teuchos_ParameterList_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Teuchos::ParameterList> robj = getConst<Teuchos::ParameterList>(aid);
            newid = CTeuchos::aliasConstParameterList(robj);
        } else {
            Teuchos::RCP<Teuchos::ParameterList> robj = get<Teuchos::ParameterList>(aid);
            newid = CTeuchos::aliasParameterList(robj);
        }
        if (!keepold) remove(&aid);
        break;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_Amesos_BaseSolver_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Amesos_BaseSolver> robj = getConst<Amesos_BaseSolver>(aid);
            newid = CAmesos::aliasConstBaseSolver(robj);
        } else {
            Teuchos::RCP<Amesos_BaseSolver> robj = get<Amesos_BaseSolver>(aid);
            newid = CAmesos::aliasBaseSolver(robj);
        }
        if (!keepold) remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_Epetra_FECrsMatrix_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_FECrsMatrix> robj = getConst<Epetra_FECrsMatrix>(aid);
            newid = CEpetra::aliasConstFECrsMatrix(robj);
        } else {
            Teuchos::RCP<Epetra_FECrsMatrix> robj = get<Epetra_FECrsMatrix>(aid);
            newid = CEpetra::aliasFECrsMatrix(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_IntSerialDenseVector_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_IntSerialDenseVector> robj = getConst<Epetra_IntSerialDenseVector>(aid);
            newid = CEpetra::aliasConstIntSerialDenseVector(robj);
        } else {
            Teuchos::RCP<Epetra_IntSerialDenseVector> robj = get<Epetra_IntSerialDenseVector>(aid);
            newid = CEpetra::aliasIntSerialDenseVector(robj);
        }
        if (!keepold) remove(&aid);
        break;
    case CT_Epetra_SerialDenseMatrix_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_SerialDenseMatrix> robj = getConst<Epetra_SerialDenseMatrix>(aid);
            newid = CEpetra::aliasConstSerialDenseMatrix(robj);
        } else {
            Teuchos::RCP<Epetra_SerialDenseMatrix> robj = get<Epetra_SerialDenseMatrix>(aid);
            newid = CEpetra::aliasSerialDenseMatrix(robj);
        }
        if (!keepold) remove(&aid);
        break;
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTest_ID:
        if (aid.is_const) {
            Teuchos::RCP<const AztecOO_StatusTest> robj = getConst<AztecOO_StatusTest>(aid);
            newid = CAztecOO::aliasConstStatusTest(robj);
        } else {
            Teuchos::RCP<AztecOO_StatusTest> robj = get<AztecOO_StatusTest>(aid);
            newid = CAztecOO::aliasStatusTest(robj);
        }
        if (!keepold) remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestCombo_ID:
        if (aid.is_const) {
            Teuchos::RCP<const AztecOO_StatusTestCombo> robj = getConst<AztecOO_StatusTestCombo>(aid);
            newid = CAztecOO::aliasConstStatusTestCombo(robj);
        } else {
            Teuchos::RCP<AztecOO_StatusTestCombo> robj = get<AztecOO_StatusTestCombo>(aid);
            newid = CAztecOO::aliasStatusTestCombo(robj);
        }
        if (!keepold) remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestMaxIters_ID:
        if (aid.is_const) {
            Teuchos::RCP<const AztecOO_StatusTestMaxIters> robj = getConst<AztecOO_StatusTestMaxIters>(aid);
            newid = CAztecOO::aliasConstStatusTestMaxIters(robj);
        } else {
            Teuchos::RCP<AztecOO_StatusTestMaxIters> robj = get<AztecOO_StatusTestMaxIters>(aid);
            newid = CAztecOO::aliasStatusTestMaxIters(robj);
        }
        if (!keepold) remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    case CT_AztecOO_StatusTestResNorm_ID:
        if (aid.is_const) {
            Teuchos::RCP<const AztecOO_StatusTestResNorm> robj = getConst<AztecOO_StatusTestResNorm>(aid);
            newid = CAztecOO::aliasConstStatusTestResNorm(robj);
        } else {
            Teuchos::RCP<AztecOO_StatusTestResNorm> robj = get<AztecOO_StatusTestResNorm>(aid);
            newid = CAztecOO::aliasStatusTestResNorm(robj);
        }
        if (!keepold) remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    case CT_Ifpack_Preconditioner_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Ifpack_Preconditioner> robj = getConst<Ifpack_Preconditioner>(aid);
            newid = CIfpack::aliasConstPreconditioner(robj);
        } else {
            Teuchos::RCP<Ifpack_Preconditioner> robj = get<Ifpack_Preconditioner>(aid);
            newid = CIfpack::aliasPreconditioner(robj);
        }
        if (!keepold) remove(&aid);
        break;
#endif /* HAVE_CTRILINOS_IFPACK */
    case CT_Epetra_SerialDenseVector_ID:
        if (aid.is_const) {
            Teuchos::RCP<const Epetra_SerialDenseVector> robj = getConst<Epetra_SerialDenseVector>(aid);
            newid = CEpetra::aliasConstSerialDenseVector(robj);
        } else {
            Teuchos::RCP<Epetra_SerialDenseVector> robj = get<Epetra_SerialDenseVector>(aid);
            newid = CEpetra::aliasSerialDenseVector(robj);
        }
        if (!keepold) remove(&aid);
        break;
    default:
        throw CTrilinosInvalidTypeError("invalid table id or non-polymorphic class");
    }

    return newid;
}

bool TableRepos::typeCheck(CTrilinos_Universal_ID_t aid, CTrilinos_Table_ID_t type)
{
    if (aid.table == type)
        return true;

    try {
        if (aid.is_const) {
            switch (type) {
            case CT_Epetra_Distributor_ID:
                getConst<Epetra_Distributor>(aid);
                break;
            case CT_Epetra_SerialComm_ID:
                getConst<Epetra_SerialComm>(aid);
                break;
            case CT_Epetra_BLAS_ID:
                getConst<Epetra_BLAS>(aid);
                break;
            case CT_Epetra_Comm_ID:
                getConst<Epetra_Comm>(aid);
                break;
            case CT_Epetra_Operator_ID:
                getConst<Epetra_Operator>(aid);
                break;
            case CT_Epetra_MultiVector_ID:
                getConst<Epetra_MultiVector>(aid);
                break;
            case CT_Epetra_OffsetIndex_ID:
                getConst<Epetra_OffsetIndex>(aid);
                break;
            case CT_Epetra_Object_ID:
                getConst<Epetra_Object>(aid);
                break;
            case CT_Epetra_RowMatrix_ID:
                getConst<Epetra_RowMatrix>(aid);
                break;
            case CT_Epetra_CompObject_ID:
                getConst<Epetra_CompObject>(aid);
                break;
            case CT_Epetra_Directory_ID:
                getConst<Epetra_Directory>(aid);
                break;
            case CT_Epetra_Flops_ID:
                getConst<Epetra_Flops>(aid);
                break;
            case CT_Epetra_SrcDistObject_ID:
                getConst<Epetra_SrcDistObject>(aid);
                break;
#ifdef HAVE_MPI
            case CT_Epetra_MpiComm_ID:
                getConst<Epetra_MpiComm>(aid);
                break;
#endif /* HAVE_MPI */
            case CT_Epetra_CrsMatrix_ID:
                getConst<Epetra_CrsMatrix>(aid);
                break;
            case CT_Epetra_CrsGraph_ID:
                getConst<Epetra_CrsGraph>(aid);
                break;
            case CT_Epetra_DistObject_ID:
                getConst<Epetra_DistObject>(aid);
                break;
            case CT_Epetra_Vector_ID:
                getConst<Epetra_Vector>(aid);
                break;
            case CT_Epetra_Export_ID:
                getConst<Epetra_Export>(aid);
                break;
            case CT_Epetra_Map_ID:
                getConst<Epetra_Map>(aid);
                break;
            case CT_Epetra_BlockMap_ID:
                getConst<Epetra_BlockMap>(aid);
                break;
            case CT_Epetra_Import_ID:
                getConst<Epetra_Import>(aid);
                break;
            case CT_Epetra_Time_ID:
                getConst<Epetra_Time>(aid);
                break;
            case CT_Epetra_JadMatrix_ID:
                getConst<Epetra_JadMatrix>(aid);
                break;
            case CT_Epetra_LinearProblem_ID:
                getConst<Epetra_LinearProblem>(aid);
                break;
            case CT_Epetra_LAPACK_ID:
                getConst<Epetra_LAPACK>(aid);
                break;
            case CT_Teuchos_ParameterList_ID:
                getConst<Teuchos::ParameterList>(aid);
                break;
#ifdef HAVE_CTRILINOS_AMESOS
            case CT_Amesos_BaseSolver_ID:
                getConst<Amesos_BaseSolver>(aid);
                break;
#endif /* HAVE_CTRILINOS_AMESOS */
            case CT_Epetra_FECrsMatrix_ID:
                getConst<Epetra_FECrsMatrix>(aid);
                break;
            case CT_Epetra_IntSerialDenseVector_ID:
                getConst<Epetra_IntSerialDenseVector>(aid);
                break;
            case CT_Epetra_SerialDenseMatrix_ID:
                getConst<Epetra_SerialDenseMatrix>(aid);
                break;
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTest_ID:
                getConst<AztecOO_StatusTest>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestCombo_ID:
                getConst<AztecOO_StatusTestCombo>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestMaxIters_ID:
                getConst<AztecOO_StatusTestMaxIters>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestResNorm_ID:
                getConst<AztecOO_StatusTestResNorm>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
            case CT_Ifpack_Preconditioner_ID:
                getConst<Ifpack_Preconditioner>(aid);
                break;
#endif /* HAVE_CTRILINOS_IFPACK */
            case CT_Epetra_SerialDenseVector_ID:
                getConst<Epetra_SerialDenseVector>(aid);
                break;
            default:
                return false;
            }
        } else {
            switch (type) {
            case CT_Epetra_Distributor_ID:
                get<Epetra_Distributor>(aid);
                break;
            case CT_Epetra_SerialComm_ID:
                get<Epetra_SerialComm>(aid);
                break;
            case CT_Epetra_BLAS_ID:
                get<Epetra_BLAS>(aid);
                break;
            case CT_Epetra_Comm_ID:
                get<Epetra_Comm>(aid);
                break;
            case CT_Epetra_Operator_ID:
                get<Epetra_Operator>(aid);
                break;
            case CT_Epetra_MultiVector_ID:
                get<Epetra_MultiVector>(aid);
                break;
            case CT_Epetra_OffsetIndex_ID:
                get<Epetra_OffsetIndex>(aid);
                break;
            case CT_Epetra_Object_ID:
                get<Epetra_Object>(aid);
                break;
            case CT_Epetra_RowMatrix_ID:
                get<Epetra_RowMatrix>(aid);
                break;
            case CT_Epetra_CompObject_ID:
                get<Epetra_CompObject>(aid);
                break;
            case CT_Epetra_Directory_ID:
                get<Epetra_Directory>(aid);
                break;
            case CT_Epetra_Flops_ID:
                get<Epetra_Flops>(aid);
                break;
            case CT_Epetra_SrcDistObject_ID:
                get<Epetra_SrcDistObject>(aid);
                break;
#ifdef HAVE_MPI
            case CT_Epetra_MpiComm_ID:
                get<Epetra_MpiComm>(aid);
                break;
#endif /* HAVE_MPI */
            case CT_Epetra_CrsMatrix_ID:
                get<Epetra_CrsMatrix>(aid);
                break;
            case CT_Epetra_CrsGraph_ID:
                get<Epetra_CrsGraph>(aid);
                break;
            case CT_Epetra_DistObject_ID:
                get<Epetra_DistObject>(aid);
                break;
            case CT_Epetra_Vector_ID:
                get<Epetra_Vector>(aid);
                break;
            case CT_Epetra_Export_ID:
                get<Epetra_Export>(aid);
                break;
            case CT_Epetra_Map_ID:
                get<Epetra_Map>(aid);
                break;
            case CT_Epetra_BlockMap_ID:
                get<Epetra_BlockMap>(aid);
                break;
            case CT_Epetra_Import_ID:
                get<Epetra_Import>(aid);
                break;
            case CT_Epetra_Time_ID:
                get<Epetra_Time>(aid);
                break;
            case CT_Epetra_JadMatrix_ID:
                get<Epetra_JadMatrix>(aid);
                break;
            case CT_Epetra_LinearProblem_ID:
                get<Epetra_LinearProblem>(aid);
                break;
            case CT_Epetra_LAPACK_ID:
                get<Epetra_LAPACK>(aid);
                break;
            case CT_Teuchos_ParameterList_ID:
                get<Teuchos::ParameterList>(aid);
                break;
#ifdef HAVE_CTRILINOS_AMESOS
            case CT_Amesos_BaseSolver_ID:
                get<Amesos_BaseSolver>(aid);
                break;
#endif /* HAVE_CTRILINOS_AMESOS */
            case CT_Epetra_FECrsMatrix_ID:
                get<Epetra_FECrsMatrix>(aid);
                break;
            case CT_Epetra_IntSerialDenseVector_ID:
                get<Epetra_IntSerialDenseVector>(aid);
                break;
            case CT_Epetra_SerialDenseMatrix_ID:
                get<Epetra_SerialDenseMatrix>(aid);
                break;
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTest_ID:
                get<AztecOO_StatusTest>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestCombo_ID:
                get<AztecOO_StatusTestCombo>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestMaxIters_ID:
                get<AztecOO_StatusTestMaxIters>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
            case CT_AztecOO_StatusTestResNorm_ID:
                get<AztecOO_StatusTestResNorm>(aid);
                break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
            case CT_Ifpack_Preconditioner_ID:
                get<Ifpack_Preconditioner>(aid);
                break;
#endif /* HAVE_CTRILINOS_IFPACK */
            case CT_Epetra_SerialDenseVector_ID:
                get<Epetra_SerialDenseVector>(aid);
                break;
            default:
                return false;
            }
        }
    } catch (...) {
        return false;
    }

    return true;
}


} // namespace CTrilinos


