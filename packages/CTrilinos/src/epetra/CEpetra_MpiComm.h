#ifndef CEPETRA_MPICOMM_H
#define CEPETRA_MPICOMM_H

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


#ifdef HAVE_MPI



/*! @file CEpetra_MpiComm.h
 * @brief Wrappers for Epetra_MpiComm */

/* True C header file! */

#include "mpi.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_MpiComm_Generalize ( 
  CT_Epetra_MpiComm_ID_t id );

/*@}*/

/*! @name Epetra_MpiComm constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_MpiComm::Epetra_MpiComm(MPI_Comm comm)
*/
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create ( MPI_Comm comm );

/*! @brief Wrapper for 
   Epetra_MpiComm::Epetra_MpiComm(const Epetra_MpiComm & Comm)
*/
CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( 
  CT_Epetra_MpiComm_ID_t CommID );

/*@}*/

/*! @name Epetra_MpiComm destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_MpiComm::~Epetra_MpiComm()
*/
void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );

/*@}*/

/*! @name Epetra_MpiComm member wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Comm * Epetra_MpiComm::Clone() const
*/
CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( 
  CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_MpiComm::Barrier() const
*/
void Epetra_MpiComm_Barrier ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(double * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, int Count, 
  int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(int * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(long * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::Broadcast(char * MyVals, int Count, int Root) const
*/
int Epetra_MpiComm_Broadcast_Char ( 
  CT_Epetra_MpiComm_ID_t selfID, char * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GatherAll(double * MyVals, double * AllVals, int Count) const
*/
int Epetra_MpiComm_GatherAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GatherAll(int * MyVals, int * AllVals, int Count) const
*/
int Epetra_MpiComm_GatherAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GatherAll(long * MyVals, long * AllVals, int Count) const
*/
int Epetra_MpiComm_GatherAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::SumAll(double * PartialSums, double * GlobalSums, int Count) const
*/
int Epetra_MpiComm_SumAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::SumAll(int * PartialSums, int * GlobalSums, int Count) const
*/
int Epetra_MpiComm_SumAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::SumAll(long * PartialSums, long * GlobalSums, int Count) const
*/
int Epetra_MpiComm_SumAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const
*/
int Epetra_MpiComm_MaxAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const
*/
int Epetra_MpiComm_MaxAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const
*/
int Epetra_MpiComm_MaxAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MinAll(double * PartialMins, double * GlobalMins, int Count) const
*/
int Epetra_MpiComm_MinAll_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MinAll(int * PartialMins, int * GlobalMins, int Count) const
*/
int Epetra_MpiComm_MinAll_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MinAll(long * PartialMins, long * GlobalMins, int Count) const
*/
int Epetra_MpiComm_MinAll_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::ScanSum(double * MyVals, double * ScanSums, int Count) const
*/
int Epetra_MpiComm_ScanSum_Double ( 
  CT_Epetra_MpiComm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::ScanSum(int * MyVals, int * ScanSums, int Count) const
*/
int Epetra_MpiComm_ScanSum_Int ( 
  CT_Epetra_MpiComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_MpiComm::ScanSum(long * MyVals, long * ScanSums, int Count) const
*/
int Epetra_MpiComm_ScanSum_Long ( 
  CT_Epetra_MpiComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   MPI_Comm Epetra_MpiComm::Comm() const
*/
MPI_Comm Epetra_MpiComm_Comm ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::MyPID() const
*/
int Epetra_MpiComm_MyPID ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::NumProc() const
*/
int Epetra_MpiComm_NumProc ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Distributor * Epetra_MpiComm::CreateDistributor() const
*/
CT_Epetra_Distributor_ID_t Epetra_MpiComm_CreateDistributor ( 
  CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Directory * Epetra_MpiComm::CreateDirectory(const Epetra_BlockMap & Map) const
*/
CT_Epetra_Directory_ID_t Epetra_MpiComm_CreateDirectory ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/*! @brief Wrapper for 
   int Epetra_MpiComm::GetMpiTag() const
*/
int Epetra_MpiComm_GetMpiTag ( CT_Epetra_MpiComm_ID_t selfID );

/*! @brief Wrapper for 
   MPI_Comm Epetra_MpiComm::GetMpiComm() const
*/
MPI_Comm Epetra_MpiComm_GetMpiComm ( CT_Epetra_MpiComm_ID_t selfID );

/*@}*/

/*! @name Epetra_MpiComm operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_MpiComm & Epetra_MpiComm::operator=(const Epetra_MpiComm & Comm)
*/
void Epetra_MpiComm_Assign ( 
  CT_Epetra_MpiComm_ID_t selfID, CT_Epetra_MpiComm_ID_t CommID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* HAVE_MPI */
#endif /* CEPETRA_MPICOMM_H */

