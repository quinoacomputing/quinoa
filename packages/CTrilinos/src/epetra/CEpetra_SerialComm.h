#ifndef CEPETRA_SERIALCOMM_H
#define CEPETRA_SERIALCOMM_H

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


/*! @file CEpetra_SerialComm.h
 * @brief Wrappers for Epetra_SerialComm */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_SerialComm_Generalize ( 
  CT_Epetra_SerialComm_ID_t id );

/*@}*/

/*! @name Epetra_SerialComm constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_SerialComm::Epetra_SerialComm()
*/
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  );

/*! @brief Wrapper for 
   Epetra_SerialComm::Epetra_SerialComm(const Epetra_SerialComm& Comm)
*/
CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( 
  CT_Epetra_SerialComm_ID_t CommID );

/*@}*/

/*! @name Epetra_SerialComm destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_SerialComm::~Epetra_SerialComm()
*/
void Epetra_SerialComm_Destroy ( CT_Epetra_SerialComm_ID_t * selfID );

/*@}*/

/*! @name Epetra_SerialComm member wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_Comm * Epetra_SerialComm::Clone() const
*/
CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( 
  CT_Epetra_SerialComm_ID_t selfID );

/*! @brief Wrapper for 
   void Epetra_SerialComm::Barrier() const
*/
void Epetra_SerialComm_Barrier ( CT_Epetra_SerialComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_SerialComm::Broadcast(double * MyVals, int Count, int Root) const
*/
int Epetra_SerialComm_Broadcast_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, int Count, 
  int Root );

/*! @brief Wrapper for 
   int Epetra_SerialComm::Broadcast(int * MyVals, int Count, int Root) const
*/
int Epetra_SerialComm_Broadcast_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int Count, 
  int Root );

/*! @brief Wrapper for 
   int Epetra_SerialComm::Broadcast(long * MyVals, int Count, int Root) const
*/
int Epetra_SerialComm_Broadcast_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, int Count, 
  int Root );

/*! @brief Wrapper for 
   int Epetra_SerialComm::Broadcast(char * MyVals, int Count, int Root) const
*/
int Epetra_SerialComm_Broadcast_Char ( 
  CT_Epetra_SerialComm_ID_t selfID, char * MyVals, int Count, 
  int Root );

/*! @brief Wrapper for 
   int Epetra_SerialComm::GatherAll(double * MyVals, double * AllVals, int Count) const
*/
int Epetra_SerialComm_GatherAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * AllVals, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::GatherAll(int * MyVals, int * AllVals, int Count) const
*/
int Epetra_SerialComm_GatherAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::GatherAll(long * MyVals, long * AllVals, int Count) const
*/
int Epetra_SerialComm_GatherAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::SumAll(double * PartialSums, double * GlobalSums, int Count) const
*/
int Epetra_SerialComm_SumAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::SumAll(int * PartialSums, int * GlobalSums, int Count) const
*/
int Epetra_SerialComm_SumAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialSums, 
  int * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::SumAll(long * PartialSums, long * GlobalSums, int Count) const
*/
int Epetra_SerialComm_SumAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialSums, 
  long * GlobalSums, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const
*/
int Epetra_SerialComm_MaxAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const
*/
int Epetra_SerialComm_MaxAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMaxs, 
  int * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const
*/
int Epetra_SerialComm_MaxAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMaxs, 
  long * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::MinAll(double * PartialMins, double * GlobalMins, int Count) const
*/
int Epetra_SerialComm_MinAll_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::MinAll(int * PartialMins, int * GlobalMins, int Count) const
*/
int Epetra_SerialComm_MinAll_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * PartialMins, 
  int * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::MinAll(long * PartialMins, long * GlobalMins, int Count) const
*/
int Epetra_SerialComm_MinAll_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * PartialMins, 
  long * GlobalMins, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::ScanSum(double * MyVals, double * ScanSums, int Count) const
*/
int Epetra_SerialComm_ScanSum_Double ( 
  CT_Epetra_SerialComm_ID_t selfID, double * MyVals, 
  double * ScanSums, int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::ScanSum(int * MyVals, int * ScanSums, int Count) const
*/
int Epetra_SerialComm_ScanSum_Int ( 
  CT_Epetra_SerialComm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::ScanSum(long * MyVals, long * ScanSums, int Count) const
*/
int Epetra_SerialComm_ScanSum_Long ( 
  CT_Epetra_SerialComm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   int Epetra_SerialComm::MyPID() const
*/
int Epetra_SerialComm_MyPID ( CT_Epetra_SerialComm_ID_t selfID );

/*! @brief Wrapper for 
   int Epetra_SerialComm::NumProc() const
*/
int Epetra_SerialComm_NumProc ( CT_Epetra_SerialComm_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Distributor * Epetra_SerialComm::CreateDistributor() const
*/
CT_Epetra_Distributor_ID_t Epetra_SerialComm_CreateDistributor ( 
  CT_Epetra_SerialComm_ID_t selfID );

/*! @brief Wrapper for 
   Epetra_Directory * Epetra_SerialComm::CreateDirectory(const Epetra_BlockMap & Map) const
*/
CT_Epetra_Directory_ID_t Epetra_SerialComm_CreateDirectory ( 
  CT_Epetra_SerialComm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/*@}*/

/*! @name Epetra_SerialComm operator wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_SerialComm & Epetra_SerialComm::operator=(const Epetra_SerialComm & Comm)
*/
void Epetra_SerialComm_Assign ( 
  CT_Epetra_SerialComm_ID_t selfID, 
  CT_Epetra_SerialComm_ID_t CommID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_SERIALCOMM_H */

