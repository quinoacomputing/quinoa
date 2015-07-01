#ifndef CGALERI_UTILS_H
#define CGALERI_UTILS_H

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



/*! @file CGaleri_Utils.h
 * @brief Wrappers for Galeri_Utils */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name Galeri_Utils static function wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_MultiVector* Galeri::CreateCartesianCoordinates(const string CoordType, const Epetra_BlockMap* BlockMap, Teuchos::ParameterList& List)
*/
CT_Epetra_MultiVector_ID_t Galeri_Utils_CreateCartesianCoordinates ( 
  const char CoordType[], CT_Epetra_BlockMap_ID_t BlockMapID, 
  CT_Teuchos_ParameterList_ID_t ListID );

/*! @brief Wrapper for 
   void Galeri::Solve(const Epetra_LinearProblem Problem)
*/
void Galeri_Utils_Solve_LinearProblem ( 
  CT_Epetra_LinearProblem_ID_t ProblemID );

/*! @brief Wrapper for 
   void Galeri::Solve(const Epetra_RowMatrix* Matrix, const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS)
*/
void Galeri_Utils_Solve_Matrix ( 
  CT_Epetra_RowMatrix_ID_t MatrixID, 
  CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID );

/*! @brief Wrapper for 
   double Galeri::ComputeNorm(const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS)
*/
double Galeri_Utils_ComputeNorm ( 
  CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID );

/*! @brief Wrapper for 
   double Galeri::ComputeNorm(const Epetra_RowMatrix* A, const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS)
*/
double Galeri_Utils_ComputeNorm_Matrix ( 
  CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t LHSID, 
  CT_Epetra_MultiVector_ID_t RHSID );

/*! @brief Wrapper for 
   string Galeri::toString(const int& x)
*/
const char * Galeri_Utils_toString_Int ( int x );

/*! @brief Wrapper for 
   string Galeri::toString(const unsigned int& x)
*/
const char * Galeri_Utils_toString_UInt ( unsigned int x );

/*! @brief Wrapper for 
   string Galeri::toString(const double& x)
*/
const char * Galeri_Utils_toString_Double ( double x );

/*! @brief Wrapper for 
   void Galeri::GetNeighboursCartesian2d(const int i, const int nx, const int ny, int & left, int & right, int & lower, int & upper)
*/
void Galeri_Utils_GetNeighboursCartesian2d ( 
  const int i, const int nx, const int ny, int * left, int * right, 
  int * lower, int * upper );

/*! @brief Wrapper for 
   void Galeri::GetNeighboursCartesian2d(const int i, const int nx, const int ny, int& left, int& right, int& lower, int& upper, int& left2, int& right2, int& lower2, int& upper2)
*/
void Galeri_Utils_GetNeighboursCartesian2d_Both ( 
  const int i, const int nx, const int ny, int * left, int * right, 
  int * lower, int * upper, int * left2, int * right2, int * lower2, 
  int * upper2 );

/*! @brief Wrapper for 
   void Galeri::GetNeighboursCartesian3d(const int i, const int nx, const int ny, const int nz, int& left, int& right, int& lower, int& upper, int& below, int& above)
*/
void Galeri_Utils_GetNeighboursCartesian3d ( 
  const int i, const int nx, const int ny, const int nz, int * left, 
  int * right, int * lower, int * upper, int * below, int * above );

/*! @brief Wrapper for 
   void Galeri::PrintStencil2D(const Epetra_CrsMatrix* Matrix, const int nx, const int ny, int GID = -1)
*/
void Galeri_Utils_PrintStencil2D ( 
  CT_Epetra_CrsMatrix_ID_t MatrixID, const int nx, const int ny, 
  int GID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_CTRILINOS_GALERI */
#endif /* CGALERI_UTILS_H */

