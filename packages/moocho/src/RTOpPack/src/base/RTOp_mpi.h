/*
// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/* */
/* MPI declarations used by RTOp example program. */
/* These where taken from mpich for Windows NT. */
/* */

#ifndef RTOP_MPI_H
#define RTOP_MPI_H

#include "RTOp_ConfigDefs.hpp" /* This C++ file has a proper C mode */

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////// */
/* MPI declarations */

#define MPI_Aint int
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_FLOAT          ((MPI_Datatype)10)
#define MPI_DOUBLE         ((MPI_Datatype)11)
typedef int MPI_Comm;
#define MPI_COMM_WORLD 91
#define MPI_COMM_NULL      ((MPI_Comm)0)
typedef int MPI_Op;
#define MPI_OP_NULL        ((MPI_Op)0)
#define MPI_MAX            (MPI_Op)(100)
#define MPI_MIN            (MPI_Op)(101)
#define MPI_SUM            (MPI_Op)(102)
#define MPI_DATATYPE_NULL  ((MPI_Datatype)0)
typedef struct { int MPI_SOURCE; int MPI_TAG; int MPI_ERROR; } RTOP_MPI_Status;
typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * ); 

/* ////////////////////////// */
/* MPI functions */

#define EXPORT_MPI_API
EXPORT_MPI_API int MPI_Init(int *, char ***);
EXPORT_MPI_API int MPI_Finalize(void);
EXPORT_MPI_API int MPI_Comm_size(MPI_Comm, int *);
EXPORT_MPI_API int MPI_Comm_rank(MPI_Comm, int *);
EXPORT_MPI_API int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
EXPORT_MPI_API int MPI_Type_commit(MPI_Datatype *);
EXPORT_MPI_API int MPI_Type_free(MPI_Datatype *);
EXPORT_MPI_API int MPI_Op_create(MPI_User_function *, int, MPI_Op *);
EXPORT_MPI_API int MPI_Op_free( MPI_Op *);
EXPORT_MPI_API int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
EXPORT_MPI_API int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, RTOP_MPI_Status*);
EXPORT_MPI_API int MPI_Sendrecv_replace(void*, int, MPI_Datatype, int, int, int, int, MPI_Comm, RTOP_MPI_Status*);
EXPORT_MPI_API int MPI_Reduce(void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
EXPORT_MPI_API int MPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);
EXPORT_MPI_API int MPI_Barrier(MPI_Comm);
EXPORT_MPI_API int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );
EXPORT_MPI_API int MPI_Gather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm); 

#ifdef __cplusplus
}
#endif

#endif /* RTOP_MPI_H */
