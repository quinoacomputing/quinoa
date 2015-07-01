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

/*
 * This file includes selected hollow MPI function definitions for a
 * sinlge process implementation.
 */

#include <assert.h>

/*
 * RAB: 2004/01/22: This file is included because it includes
 * Thyra_Config.h which then defines RTOp_USE_MPI or not.  If
 * RTOp_USE_MPI is defined then this header file will also include
 * RTOp_mpi.h for these delcarations.
 */
#include "RTOp_MPI_config.h"

#ifndef RTOp_USE_MPI

int MPI_Init(int *argc, char ***argv)
{
  return 0;
}

int MPI_Finalize(void)
{
  return 0;
}

int MPI_Comm_size(MPI_Comm comm, int *size)
{
  *size = 1;
  return 0;
}

int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  *rank = 0;
  return 0;
}

int MPI_Type_struct(int count , int *array_of_blocklengths, MPI_Aint *array_of_displacements
  , MPI_Datatype *array_of_types, MPI_Datatype *data_type)
{
  /* Make the mpi datatype just the extent (needed latter!) */
  int len = 0, extent = 0, k = 0;
  for( k = 0; k < count; ++k ) {
    switch( array_of_types[k] ) {
      case MPI_CHAR:
        len = sizeof(char);
        break;
      case MPI_INT:
        len = sizeof(int);
        break;
      case MPI_FLOAT:
        len = sizeof(float);
        break;
      case MPI_DOUBLE:
        len = sizeof(double);
        break;
      default:
        assert(0);
    }
    len = array_of_displacements[k] + array_of_blocklengths[k] * len;
    if( len > extent )
      extent = len;
  }
  *data_type = extent;
  return 0;
}

int MPI_Type_commit(MPI_Datatype *datatype)
{
  return 0;
}

int MPI_Type_free(MPI_Datatype *op)
{
  *op = MPI_DATATYPE_NULL;
  return 0;
}

int MPI_Op_create(MPI_User_function *func, int communitive, MPI_Op *op)
{
  *op = (MPI_Op)*func;
  return 0;
}
int MPI_Op_free( MPI_Op *op)
{
  *op = MPI_OP_NULL;
  return 0;
}

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
  assert(0); /* Should never be called in serial mode */
  return 0;
}
int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, RTOP_MPI_Status* status)
{
  assert(0); /* Should never be called in serial mode */
  return 0;
}

int MPI_Sendrecv_replace(void* buff, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, RTOP_MPI_Status* status)
{
  assert(0); /* Should never be called in serial mode */
  return 0;
}

int MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op
  , int root, MPI_Comm comm)
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  for( k = 0; k < count * datatype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype
  , MPI_Op op, MPI_Comm comm)
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  for( k = 0; k < count * datatype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

int MPI_Barrier(MPI_Comm comm)
{
  return 0;
}

int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
{
  return 0;
}

int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype
         , void* recvbuf, int recvcount, MPI_Datatype recvtype, int root , MPI_Comm comm )
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  assert(sendtype == recvtype);
  assert(sendcount == recvcount);
  for( k = 0; k < sendcount * sendtype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

#endif /* RTOp_USE_MPI */
