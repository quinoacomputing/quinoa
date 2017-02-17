/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: rrdrake $
 *    Date: 2009/07/17 15:14:49 $
 *    Revision: 1.3 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ***********************  PMPI_Type_size.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_size ( MPI_Datatype datatype, int *size )
{
  int index;
  *size = 0;
  if ( _MPI_BasicType(datatype) == MPI_SUCCESS ) {
    *size = _MPI_getSize(datatype);
  }
  else
  {
    index = _MPI_FindType (datatype);
    if (index == _MPI_NOT_OK)
    {
      _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI_TYPE_SIZE: datatype error");
      MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE); 
    }
    *size = _MPI_TYPE_LIST[index].size;
  }
  return MPI_SUCCESS;
}

