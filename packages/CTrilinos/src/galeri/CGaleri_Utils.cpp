
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


#ifdef HAVE_CTRILINOS_GALERI


#include "CTrilinos_enums.h"
#include "CGaleri_Utils.h"
#include "CGaleri_Utils_Cpp.hpp"
#include "Galeri_Utils.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"


//
// Definitions from CGaleri_Utils.h
//


extern "C" {


CT_Epetra_MultiVector_ID_t Galeri_Utils_CreateCartesianCoordinates ( 
  const char CoordType[], CT_Epetra_BlockMap_ID_t BlockMapID, 
  CT_Teuchos_ParameterList_ID_t ListID )
{
    const Teuchos::RCP<const Epetra_BlockMap> BlockMap = 
        CEpetra::getConstBlockMap(BlockMapID);
    const Teuchos::RCP<Teuchos::ParameterList> List = 
        CTeuchos::getParameterList(ListID);
    return CEpetra::storeMultiVector(Galeri::CreateCartesianCoordinates(
        std::string(CoordType), BlockMap.getRawPtr(), *List));
}

void Galeri_Utils_Solve_LinearProblem ( 
  CT_Epetra_LinearProblem_ID_t ProblemID )
{
    const Teuchos::RCP<const Epetra_LinearProblem> Problem = 
        CEpetra::getConstLinearProblem(ProblemID);
    Galeri::Solve(*Problem);
}

void Galeri_Utils_Solve_Matrix ( 
  CT_Epetra_RowMatrix_ID_t MatrixID, 
  CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID )
{
    const Teuchos::RCP<const Epetra_RowMatrix> Matrix = 
        CEpetra::getConstRowMatrix(MatrixID);
    const Teuchos::RCP<const Epetra_MultiVector> LHS = 
        CEpetra::getConstMultiVector(LHSID);
    const Teuchos::RCP<const Epetra_MultiVector> RHS = 
        CEpetra::getConstMultiVector(RHSID);
    Galeri::Solve(Matrix.getRawPtr(), LHS.getRawPtr(), RHS.getRawPtr());
}

double Galeri_Utils_ComputeNorm ( 
  CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID )
{
    const Teuchos::RCP<const Epetra_MultiVector> LHS = 
        CEpetra::getConstMultiVector(LHSID);
    const Teuchos::RCP<const Epetra_MultiVector> RHS = 
        CEpetra::getConstMultiVector(RHSID);
    return Galeri::ComputeNorm(LHS.getRawPtr(), RHS.getRawPtr());
}

double Galeri_Utils_ComputeNorm_Matrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID )
{
    const Teuchos::RCP<const Epetra_RowMatrix> A = CEpetra::getConstRowMatrix(
        AID);
    const Teuchos::RCP<const Epetra_MultiVector> LHS = 
        CEpetra::getConstMultiVector(LHSID);
    const Teuchos::RCP<const Epetra_MultiVector> RHS = 
        CEpetra::getConstMultiVector(RHSID);
    return Galeri::ComputeNorm(A.getRawPtr(), LHS.getRawPtr(), 
        RHS.getRawPtr());
}

const char * Galeri_Utils_toString_Int ( int x )
{
    return Galeri::toString(x).c_str();
}

const char * Galeri_Utils_toString_UInt ( unsigned int x )
{
    return Galeri::toString(x).c_str();
}

const char * Galeri_Utils_toString_Double ( double x )
{
    return Galeri::toString(x).c_str();
}

void Galeri_Utils_GetNeighboursCartesian2d ( 
  const int i, const int nx, const int ny, int * left, int * right, 
  int * lower, int * upper )
{
    Galeri::GetNeighboursCartesian2d(i, nx, ny, *left, *right, *lower, *upper);
}

void Galeri_Utils_GetNeighboursCartesian2d_Both ( 
  const int i, const int nx, const int ny, int * left, int * right, 
  int * lower, int * upper, int * left2, int * right2, int * lower2, 
  int * upper2 )
{
    Galeri::GetNeighboursCartesian2d(i, nx, ny, *left, *right, *lower, *upper, 
        *left2, *right2, *lower2, *upper2);
}

void Galeri_Utils_GetNeighboursCartesian3d ( 
  const int i, const int nx, const int ny, const int nz, int * left, 
  int * right, int * lower, int * upper, int * below, int * above )
{
    Galeri::GetNeighboursCartesian3d(i, nx, ny, nz, *left, *right, *lower, 
        *upper, *below, *above);
}

void Galeri_Utils_PrintStencil2D ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID, const int nx, const int ny, 
  int GID )
{
    const Teuchos::RCP<const Epetra_CrsMatrix> Matrix = 
        CEpetra::getConstCrsMatrix(MatrixID);
    Galeri::PrintStencil2D(Matrix.getRawPtr(), nx, ny, GID);
}


} // extern "C"


//
// Definitions from CGaleri_Utils_Cpp.hpp
//




#endif /* HAVE_CTRILINOS_GALERI */


