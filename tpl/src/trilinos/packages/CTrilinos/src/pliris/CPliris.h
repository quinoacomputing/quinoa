
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



/*! @file CPliris.h
 * @brief Wrappers for Pliris */

/* True C header file! */


#ifndef CPLIRIS_H
#define CPLIRIS_H

#ifdef HAVE_CTRILINOS_PLIRIS

#ifdef HAVE_MPI


#include "CTrilinos_config.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name Pliris constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Pliris::Pliris(Epetra_Vector * A, Epetra_MultiVector * X, Epetra_MultiVector * B)
*/
CT_Pliris_ID_t Pliris_Create ( 
  CT_Epetra_Vector_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
  CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   Pliris::Pliris()
*/
CT_Pliris_ID_t Pliris_Create_Default (  );

/*@}*/

/*! @name Pliris destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Pliris::~Pliris(void)
*/
void Pliris_Destroy ( CT_Pliris_ID_t * selfID );

/*@}*/

/*! @name Pliris member wrappers */
/*@{*/

/*! @brief Wrapper for 
   int Pliris::SetLHS(Epetra_MultiVector * X)
*/
int Pliris_SetLHS ( 
  CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t XID );

/*! @brief Wrapper for 
   int Pliris::SetRHS(Epetra_MultiVector * B)
*/
int Pliris_SetRHS ( 
  CT_Pliris_ID_t selfID, CT_Epetra_MultiVector_ID_t BID );

/*! @brief Wrapper for 
   int Pliris::SetMatrix(Epetra_Vector * A)
*/
int Pliris_SetMatrix ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID );

/*! @brief Wrapper for 
   int Pliris::SetMatrix(Epetra_SerialDenseVector * A)
*/
int Pliris_SetMatrix_Serial ( 
  CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AID );

/*! @brief Wrapper for 
   int Pliris::GetDistribution( int * nprocs_row, int * number_of_unknowns, int * nrhs, int * my_rows, int * my_cols, int * my_first_row, int * my_first_col, int * my_rhs, int * my_row, int * my_col )
*/
int Pliris_GetDistribution ( 
  CT_Pliris_ID_t selfID, int * nprocs_row, int * number_of_unknowns, 
  int * nrhs, int * my_rows, int * my_cols, int * my_first_row, 
  int * my_first_col, int * my_rhs, int * my_row, int * my_col );

/*! @brief Wrapper for 
   int Pliris::FactorSolve( Epetra_Vector * A, int my_rows, int my_cols, int* matrix_size, int* num_procsr, int* num_rhs, double* secs)
*/
int Pliris_FactorSolve ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, int my_rows, 
  int my_cols, int * matrix_size, int * num_procsr, int * num_rhs, 
  double * secs );

/*! @brief Wrapper for 
   int Pliris::FactorSolve( Epetra_SerialDenseVector * AA, int my_rows, int my_cols, int* matrix_size, int* num_procsr, int* num_rhs, double* secs)
*/
int Pliris_FactorSolve_Serial ( 
  CT_Pliris_ID_t selfID, CT_Epetra_SerialDenseVector_ID_t AAID, 
  int my_rows, int my_cols, int * matrix_size, int * num_procsr, 
  int * num_rhs, double * secs );

/*! @brief Wrapper for 
   int Pliris::Factor( Epetra_Vector* A, int* matrix_size, int* num_procsr, int* permute, double* secs)
*/
int Pliris_Factor ( 
  CT_Pliris_ID_t selfID, CT_Epetra_Vector_ID_t AID, 
  int * matrix_size, int * num_procsr, int * permute, 
  double * secs );

/*! @brief Wrapper for 
   int Pliris::Solve(int* permute, int* num_rhs)
*/
int Pliris_Solve ( 
  CT_Pliris_ID_t selfID, int * permute, int * num_rhs );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_MPI */

#endif /* HAVE_CTRILINOS_PLIRIS */

#endif /* CPLIRIS_H */

