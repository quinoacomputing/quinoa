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

#include "dr_const.h"
#include "dr_eval_const.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Utility functions for evaluating a partition.
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void driver_eval(MESH_INFO_PTR mesh)
{
/*
 * Function to evaluate a partition.  Largely duplicates functionality
 * of Zoltan_LB_Eval, but provides sanity checking.
 *
 * Currently uses only the first cpu weight.
 */
int i;
int proc;
double load = 0.;
ZOLTAN_ID_TYPE cuts = 0;
ZOLTAN_ID_TYPE gsumcuts, gmaxcuts, gmincuts, elemcount;
ZOLTAN_ID_TYPE gsumelems, gmaxelems, gminelems;
double gsumload, gmaxload, gminload;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  for (i = 0; i < mesh->necmap; i++) {
    cuts += mesh->ecmap_cnt[i];
  }
  
  for (i = 0; i < mesh->num_elems; i++) {
    if (mesh->blank_count && (mesh->blank[i] == 1)) continue;
    load += mesh->elements[i].cpu_wgt[0];
  }

  MPI_Allreduce(&cuts, &gsumcuts, 1, ZOLTAN_ID_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&cuts, &gmaxcuts, 1, ZOLTAN_ID_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&cuts, &gmincuts, 1, ZOLTAN_ID_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD);

  elemcount = mesh->num_elems - mesh->blank_count;

  MPI_Allreduce(&elemcount, &gsumelems, 1, ZOLTAN_ID_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&elemcount, &gmaxelems, 1, ZOLTAN_ID_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&elemcount, &gminelems, 1, ZOLTAN_ID_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD);

  MPI_Allreduce(&load, &gsumload, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&load, &gmaxload, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&load, &gminload, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  if (proc == 0) {
    printf("DRIVER EVAL:  load:  max %f  min %f  sum %f\n", 
           gmaxload, gminload, gsumload);
    printf("DRIVER EVAL:  objs:  max " ZOLTAN_ID_SPEC "  min " ZOLTAN_ID_SPEC "  sum " ZOLTAN_ID_SPEC "\n", 
           gmaxelems, gminelems, gsumelems);
    printf("DRIVER EVAL:  cuts:  max " ZOLTAN_ID_SPEC "  min " ZOLTAN_ID_SPEC "  sum " ZOLTAN_ID_SPEC "\n",
           gmaxcuts, gmincuts, gsumcuts);
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
