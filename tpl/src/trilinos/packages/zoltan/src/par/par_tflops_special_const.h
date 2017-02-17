/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifndef __PAR_TFLOPS_SPECIAL_H
#define __PAR_TFLOPS_SPECIAL_H

#include <mpi.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* prototypes for TFLOPS_SPECIAL */
extern void Zoltan_RB_scan_double(double *, double *, int, MPI_Comm, int, int, int);
extern void Zoltan_RB_sum_double(double *, int, int, int, int, MPI_Comm);
extern void Zoltan_RB_max_double(double *, int, int, int, int, MPI_Comm);

extern void Zoltan_RB_bcast_doubles(double *, int, int, int, int, int, MPI_Comm comm);
extern void Zoltan_RB_gather_double( double, double *, int, int, int, int, MPI_Comm comm);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
