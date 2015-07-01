#ifndef CEPETRA_COMM_H
#define CEPETRA_COMM_H

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


/*! @file CEpetra_Comm.h
 * @brief Wrappers for Epetra_Comm */

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
CT_Epetra_Comm_ID_t Epetra_Comm_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_Comm_Generalize ( 
  CT_Epetra_Comm_ID_t id );

/*@}*/

/*! @name Epetra_Comm destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Comm::~Epetra_Comm()
*/
void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );

/*@}*/

/*! @name Epetra_Comm member wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_Comm * Epetra_Comm::Clone() const = 0
*/
CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );

/*! @brief Wrapper for 
   virtual void Epetra_Comm::Barrier() const = 0
*/
void Epetra_Comm_Barrier ( CT_Epetra_Comm_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::Broadcast(double * MyVals, int Count, int Root) const = 0
*/
int Epetra_Comm_Broadcast_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::Broadcast(int * MyVals, int Count, int Root) const = 0
*/
int Epetra_Comm_Broadcast_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::Broadcast(long * MyVals, int Count, int Root) const = 0
*/
int Epetra_Comm_Broadcast_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::Broadcast(char * MyVals, int Count, int Root) const = 0
*/
int Epetra_Comm_Broadcast_Char ( 
  CT_Epetra_Comm_ID_t selfID, char * MyVals, int Count, int Root );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::GatherAll(double * MyVals, double * AllVals, int Count) const = 0
*/
int Epetra_Comm_GatherAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * AllVals, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::GatherAll(int * MyVals, int * AllVals, int Count) const = 0
*/
int Epetra_Comm_GatherAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * AllVals, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::GatherAll(long * MyVals, long * AllVals, int Count) const = 0
*/
int Epetra_Comm_GatherAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * AllVals, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::SumAll(double * PartialSums, double * GlobalSums, int Count) const = 0
*/
int Epetra_Comm_SumAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialSums, 
  double * GlobalSums, int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::SumAll(int * PartialSums, int * GlobalSums, int Count) const = 0
*/
int Epetra_Comm_SumAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialSums, int * GlobalSums, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::SumAll(long * PartialSums, long * GlobalSums, int Count) const = 0
*/
int Epetra_Comm_SumAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialSums, long * GlobalSums, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const = 0
*/
int Epetra_Comm_MaxAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMaxs, 
  double * GlobalMaxs, int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const = 0
*/
int Epetra_Comm_MaxAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMaxs, int * GlobalMaxs, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::MaxAll(long * PartialMaxs, long * GlobalMaxs, int Count) const = 0
*/
int Epetra_Comm_MaxAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMaxs, long * GlobalMaxs, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::MinAll(double * PartialMins, double * GlobalMins, int Count) const = 0
*/
int Epetra_Comm_MinAll_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * PartialMins, 
  double * GlobalMins, int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::MinAll(int * PartialMins, int * GlobalMins, int Count) const = 0
*/
int Epetra_Comm_MinAll_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * PartialMins, int * GlobalMins, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::MinAll(long * PartialMins, long * GlobalMins, int Count) const = 0
*/
int Epetra_Comm_MinAll_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * PartialMins, long * GlobalMins, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::ScanSum(double * MyVals, double * ScanSums, int Count) const = 0
*/
int Epetra_Comm_ScanSum_Double ( 
  CT_Epetra_Comm_ID_t selfID, double * MyVals, double * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::ScanSum(int * MyVals, int * ScanSums, int Count) const = 0
*/
int Epetra_Comm_ScanSum_Int ( 
  CT_Epetra_Comm_ID_t selfID, int * MyVals, int * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::ScanSum(long * MyVals, long * ScanSums, int Count) const = 0
*/
int Epetra_Comm_ScanSum_Long ( 
  CT_Epetra_Comm_ID_t selfID, long * MyVals, long * ScanSums, 
  int Count );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::MyPID() const = 0
*/
int Epetra_Comm_MyPID ( CT_Epetra_Comm_ID_t selfID );

/*! @brief Wrapper for 
   virtual int Epetra_Comm::NumProc() const = 0
*/
int Epetra_Comm_NumProc ( CT_Epetra_Comm_ID_t selfID );

/*! @brief Wrapper for 
   virtual Epetra_Distributor * Epetra_Comm::CreateDistributor() const = 0
*/
CT_Epetra_Distributor_ID_t Epetra_Comm_CreateDistributor ( 
  CT_Epetra_Comm_ID_t selfID );

/*! @brief Wrapper for 
   virtual Epetra_Directory * Epetra_Comm::CreateDirectory(const Epetra_BlockMap & Map) const = 0
*/
CT_Epetra_Directory_ID_t Epetra_Comm_CreateDirectory ( 
  CT_Epetra_Comm_ID_t selfID, CT_Epetra_BlockMap_ID_t MapID );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_COMM_H */

