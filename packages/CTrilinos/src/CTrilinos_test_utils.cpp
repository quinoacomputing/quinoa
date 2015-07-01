
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


/*! @file CTrilinos_test_utils.cpp
 * @brief Utility functions for CTrilinos testing. */


#include "CTrilinos_config.h"
#include "CTrilinos_test_utils.hpp"
#include "CTrilinos_TableRepos.hpp"


namespace CTrilinos {


/* isSameObject(generic_id, generic_id) */
bool
isSameObject( CTrilinos_Universal_ID_t id1, CTrilinos_Universal_ID_t id2 )
{
    bool shares = false;

    if (id1.is_const) {
        switch (id1.table) {
        case CT_Epetra_Distributor_ID:
            shares = isSameObject(CEpetra::getConstDistributor(id1), id2);
            break;
        case CT_Epetra_SerialComm_ID:
            shares = isSameObject(CEpetra::getConstSerialComm(id1), id2);
            break;
        case CT_Epetra_BLAS_ID:
            shares = isSameObject(CEpetra::getConstBLAS(id1), id2);
            break;
        case CT_Epetra_Comm_ID:
            shares = isSameObject(CEpetra::getConstComm(id1), id2);
            break;
        case CT_Epetra_Operator_ID:
            shares = isSameObject(CEpetra::getConstOperator(id1), id2);
            break;
        case CT_Epetra_MultiVector_ID:
            shares = isSameObject(CEpetra::getConstMultiVector(id1), id2);
            break;
        case CT_Epetra_OffsetIndex_ID:
            shares = isSameObject(CEpetra::getConstOffsetIndex(id1), id2);
            break;
        case CT_Epetra_Object_ID:
            shares = isSameObject(CEpetra::getConstObject(id1), id2);
            break;
        case CT_Epetra_RowMatrix_ID:
            shares = isSameObject(CEpetra::getConstRowMatrix(id1), id2);
            break;
        case CT_Epetra_CompObject_ID:
            shares = isSameObject(CEpetra::getConstCompObject(id1), id2);
            break;
        case CT_Epetra_Directory_ID:
            shares = isSameObject(CEpetra::getConstDirectory(id1), id2);
            break;
        case CT_Epetra_Flops_ID:
            shares = isSameObject(CEpetra::getConstFlops(id1), id2);
            break;
        case CT_Epetra_SrcDistObject_ID:
            shares = isSameObject(CEpetra::getConstSrcDistObject(id1), id2);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            shares = isSameObject(CEpetra::getConstMpiComm(id1), id2);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            shares = isSameObject(CEpetra::getConstCrsMatrix(id1), id2);
            break;
        case CT_Epetra_CrsGraph_ID:
            shares = isSameObject(CEpetra::getConstCrsGraph(id1), id2);
            break;
        case CT_Epetra_DistObject_ID:
            shares = isSameObject(CEpetra::getConstDistObject(id1), id2);
            break;
        case CT_Epetra_Vector_ID:
            shares = isSameObject(CEpetra::getConstVector(id1), id2);
            break;
        case CT_Epetra_Export_ID:
            shares = isSameObject(CEpetra::getConstExport(id1), id2);
            break;
        case CT_Epetra_Map_ID:
            shares = isSameObject(CEpetra::getConstMap(id1), id2);
            break;
        case CT_Epetra_BlockMap_ID:
            shares = isSameObject(CEpetra::getConstBlockMap(id1), id2);
            break;
        case CT_Epetra_Import_ID:
            shares = isSameObject(CEpetra::getConstImport(id1), id2);
            break;
        case CT_Epetra_Time_ID:
            shares = isSameObject(CEpetra::getConstTime(id1), id2);
            break;
        case CT_Epetra_JadMatrix_ID:
            shares = isSameObject(CEpetra::getConstJadMatrix(id1), id2);
            break;
        case CT_Epetra_LinearProblem_ID:
            shares = isSameObject(CEpetra::getConstLinearProblem(id1), id2);
            break;
        case CT_Epetra_LAPACK_ID:
            shares = isSameObject(CEpetra::getConstLAPACK(id1), id2);
            break;
        case CT_Teuchos_CommandLineProcessor_ID:
            shares = isSameObject(CTeuchos::getConstCommandLineProcessor(id1), id2);
            break;
        case CT_Teuchos_ParameterList_ID:
            shares = isSameObject(CTeuchos::getConstParameterList(id1), id2);
            break;
        case CT_Teuchos_ParameterEntry_ID:
            shares = isSameObject(CTeuchos::getConstParameterEntry(id1), id2);
            break;
        case CT_Teuchos_any_ID:
            shares = isSameObject(CTeuchos::getConstany(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            shares = isSameObject(CAmesos::getConstBaseSolver(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_ID:
            shares = isSameObject(CAmesos::getConstAmesos(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            shares = isSameObject(CEpetra::getConstFECrsMatrix(id1), id2);
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            shares = isSameObject(CEpetra::getConstIntSerialDenseVector(id1), id2);
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            shares = isSameObject(CEpetra::getConstSerialDenseMatrix(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_ID:
            shares = isSameObject(CAztecOO::getConstAztecOO(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            shares = isSameObject(CAztecOO::getConstStatusTest(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            shares = isSameObject(CAztecOO::getConstStatusTestCombo(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            shares = isSameObject(CAztecOO::getConstStatusTestMaxIters(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            shares = isSameObject(CAztecOO::getConstStatusTestResNorm(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_ID:
            shares = isSameObject(CIfpack::getConstIfpack(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            shares = isSameObject(CIfpack::getConstPreconditioner(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        case CT_Epetra_SerialDenseVector_ID:
            shares = isSameObject(CEpetra::getConstSerialDenseVector(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_PLIRIS
#ifdef HAVE_MPI
        case CT_Pliris_ID:
            shares = isSameObject(CPliris::getConstPliris(id1), id2);
            break;
#endif /* HAVE_MPI */
#endif /* HAVE_CTRILINOS_PLIRIS */
        default:
            break;
        }
    } else {
        switch (id1.table) {
        case CT_Epetra_Distributor_ID:
            shares = isSameObject(CEpetra::getDistributor(id1), id2);
            break;
        case CT_Epetra_SerialComm_ID:
            shares = isSameObject(CEpetra::getSerialComm(id1), id2);
            break;
        case CT_Epetra_BLAS_ID:
            shares = isSameObject(CEpetra::getBLAS(id1), id2);
            break;
        case CT_Epetra_Comm_ID:
            shares = isSameObject(CEpetra::getComm(id1), id2);
            break;
        case CT_Epetra_Operator_ID:
            shares = isSameObject(CEpetra::getOperator(id1), id2);
            break;
        case CT_Epetra_MultiVector_ID:
            shares = isSameObject(CEpetra::getMultiVector(id1), id2);
            break;
        case CT_Epetra_OffsetIndex_ID:
            shares = isSameObject(CEpetra::getOffsetIndex(id1), id2);
            break;
        case CT_Epetra_Object_ID:
            shares = isSameObject(CEpetra::getObject(id1), id2);
            break;
        case CT_Epetra_RowMatrix_ID:
            shares = isSameObject(CEpetra::getRowMatrix(id1), id2);
            break;
        case CT_Epetra_CompObject_ID:
            shares = isSameObject(CEpetra::getCompObject(id1), id2);
            break;
        case CT_Epetra_Directory_ID:
            shares = isSameObject(CEpetra::getDirectory(id1), id2);
            break;
        case CT_Epetra_Flops_ID:
            shares = isSameObject(CEpetra::getFlops(id1), id2);
            break;
        case CT_Epetra_SrcDistObject_ID:
            shares = isSameObject(CEpetra::getSrcDistObject(id1), id2);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            shares = isSameObject(CEpetra::getMpiComm(id1), id2);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            shares = isSameObject(CEpetra::getCrsMatrix(id1), id2);
            break;
        case CT_Epetra_CrsGraph_ID:
            shares = isSameObject(CEpetra::getCrsGraph(id1), id2);
            break;
        case CT_Epetra_DistObject_ID:
            shares = isSameObject(CEpetra::getDistObject(id1), id2);
            break;
        case CT_Epetra_Vector_ID:
            shares = isSameObject(CEpetra::getVector(id1), id2);
            break;
        case CT_Epetra_Export_ID:
            shares = isSameObject(CEpetra::getExport(id1), id2);
            break;
        case CT_Epetra_Map_ID:
            shares = isSameObject(CEpetra::getMap(id1), id2);
            break;
        case CT_Epetra_BlockMap_ID:
            shares = isSameObject(CEpetra::getBlockMap(id1), id2);
            break;
        case CT_Epetra_Import_ID:
            shares = isSameObject(CEpetra::getImport(id1), id2);
            break;
        case CT_Epetra_Time_ID:
            shares = isSameObject(CEpetra::getTime(id1), id2);
            break;
        case CT_Epetra_JadMatrix_ID:
            shares = isSameObject(CEpetra::getJadMatrix(id1), id2);
            break;
        case CT_Epetra_LinearProblem_ID:
            shares = isSameObject(CEpetra::getLinearProblem(id1), id2);
            break;
        case CT_Epetra_LAPACK_ID:
            shares = isSameObject(CEpetra::getLAPACK(id1), id2);
            break;
        case CT_Teuchos_CommandLineProcessor_ID:
            shares = isSameObject(CTeuchos::getCommandLineProcessor(id1), id2);
            break;
        case CT_Teuchos_ParameterList_ID:
            shares = isSameObject(CTeuchos::getParameterList(id1), id2);
            break;
        case CT_Teuchos_ParameterEntry_ID:
            shares = isSameObject(CTeuchos::getParameterEntry(id1), id2);
            break;
        case CT_Teuchos_any_ID:
            shares = isSameObject(CTeuchos::getany(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            shares = isSameObject(CAmesos::getBaseSolver(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_ID:
            shares = isSameObject(CAmesos::getAmesos(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        case CT_Epetra_FECrsMatrix_ID:
            shares = isSameObject(CEpetra::getFECrsMatrix(id1), id2);
            break;
        case CT_Epetra_IntSerialDenseVector_ID:
            shares = isSameObject(CEpetra::getIntSerialDenseVector(id1), id2);
            break;
        case CT_Epetra_SerialDenseMatrix_ID:
            shares = isSameObject(CEpetra::getSerialDenseMatrix(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_ID:
            shares = isSameObject(CAztecOO::getAztecOO(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTest_ID:
            shares = isSameObject(CAztecOO::getStatusTest(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestCombo_ID:
            shares = isSameObject(CAztecOO::getStatusTestCombo(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestMaxIters_ID:
            shares = isSameObject(CAztecOO::getStatusTestMaxIters(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
        case CT_AztecOO_StatusTestResNorm_ID:
            shares = isSameObject(CAztecOO::getStatusTestResNorm(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_ID:
            shares = isSameObject(CIfpack::getIfpack(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
        case CT_Ifpack_Preconditioner_ID:
            shares = isSameObject(CIfpack::getPreconditioner(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_IFPACK */
        case CT_Epetra_SerialDenseVector_ID:
            shares = isSameObject(CEpetra::getSerialDenseVector(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_PLIRIS
#ifdef HAVE_MPI
        case CT_Pliris_ID:
            shares = isSameObject(CPliris::getPliris(id1), id2);
            break;
#endif /* HAVE_MPI */
#endif /* HAVE_CTRILINOS_PLIRIS */
        default:
            break;
        }
    }

    return shares;
}


void
purgeAllTables(  )
{
    CEpetra::purgeDistributor();
    CEpetra::purgeSerialComm();
    CEpetra::purgeBLAS();
    CEpetra::purgeComm();
    CEpetra::purgeOperator();
    CEpetra::purgeMultiVector();
    CEpetra::purgeOffsetIndex();
    CEpetra::purgeObject();
    CEpetra::purgeRowMatrix();
    CEpetra::purgeCompObject();
    CEpetra::purgeDirectory();
    CEpetra::purgeFlops();
    CEpetra::purgeSrcDistObject();
#ifdef HAVE_MPI
    CEpetra::purgeMpiComm();
#endif /* HAVE_MPI */
    CEpetra::purgeCrsMatrix();
    CEpetra::purgeCrsGraph();
    CEpetra::purgeDistObject();
    CEpetra::purgeVector();
    CEpetra::purgeExport();
    CEpetra::purgeMap();
    CEpetra::purgeBlockMap();
    CEpetra::purgeImport();
    CEpetra::purgeTime();
    CEpetra::purgeJadMatrix();
    CEpetra::purgeLinearProblem();
    CEpetra::purgeLAPACK();
    CTeuchos::purgeCommandLineProcessor();
    CTeuchos::purgeParameterList();
    CTeuchos::purgeParameterEntry();
    CTeuchos::purgeany();
#ifdef HAVE_CTRILINOS_AMESOS
    CAmesos::purgeBaseSolver();
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
    CAmesos::purgeAmesos();
#endif /* HAVE_CTRILINOS_AMESOS */
    CEpetra::purgeFECrsMatrix();
    CEpetra::purgeIntSerialDenseVector();
    CEpetra::purgeSerialDenseMatrix();
#ifdef HAVE_CTRILINOS_AZTECOO
    CAztecOO::purgeAztecOO();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    CAztecOO::purgeStatusTest();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    CAztecOO::purgeStatusTestCombo();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    CAztecOO::purgeStatusTestMaxIters();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_AZTECOO
    CAztecOO::purgeStatusTestResNorm();
#endif /* HAVE_CTRILINOS_AZTECOO */
#ifdef HAVE_CTRILINOS_IFPACK
    CIfpack::purgeIfpack();
#endif /* HAVE_CTRILINOS_IFPACK */
#ifdef HAVE_CTRILINOS_IFPACK
    CIfpack::purgePreconditioner();
#endif /* HAVE_CTRILINOS_IFPACK */
    CEpetra::purgeSerialDenseVector();
#ifdef HAVE_CTRILINOS_PLIRIS
#ifdef HAVE_MPI
    CPliris::purgePliris();
#endif /* HAVE_MPI */
#endif /* HAVE_CTRILINOS_PLIRIS */
}



} // namespace CTrilinos


