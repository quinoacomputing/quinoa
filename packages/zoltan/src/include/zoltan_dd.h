/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef ZOLTAN_DD_DDIRECTORY_H
#define ZOLTAN_DD_DDIRECTORY_H

#include "zoltan_types.h"
#include <mpi.h>

/*
** Must define this function prototype before #ifdef __cplusplus
** to avoid warning when compiling with C++ on solaris
*/
typedef unsigned int ZOLTAN_HASH_FN(ZOLTAN_ID_PTR, int, unsigned int);

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

struct Zoltan_DD_Struct;
typedef struct Zoltan_DD_Struct Zoltan_DD_Directory;


/***********  Distributed Directory Function Prototypes ************/

int Zoltan_DD_Create (Zoltan_DD_Directory **dd, MPI_Comm comm, 
 int num_gid,  /* number of ZOLTAN_ID_TYPEs in a GID */
 int num_lid,  /* number of ZOLTAN_ID_TYPEs in an LID (optional) */
 int user_length,  /* number of chars in the user data (optional) */
 int table_length, int debug_level) ;

int Zoltan_DD_Copy_To(Zoltan_DD_Directory **toptr, Zoltan_DD_Directory *from);
Zoltan_DD_Directory *Zoltan_DD_Copy(Zoltan_DD_Directory *from);

void Zoltan_DD_Destroy (Zoltan_DD_Directory **dd) ;

int Zoltan_DD_Update (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 ZOLTAN_ID_PTR lid, char *user, int *partition, int count) ;

int Zoltan_DD_Find (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 ZOLTAN_ID_PTR lid, char *data, int *partition, int count,
 int *owner) ;

int Zoltan_DD_GetLocalKeys(Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR* gid, int* size);

int Zoltan_DD_Remove (Zoltan_DD_Directory *dd, ZOLTAN_ID_PTR gid,
 int count) ;

int Zoltan_DD_Set_Hash_Fn (Zoltan_DD_Directory *dd, ZOLTAN_HASH_FN *hash);

void Zoltan_DD_Stats (Zoltan_DD_Directory *dd) ;

int Zoltan_DD_Set_Neighbor_Hash_Fn1 (Zoltan_DD_Directory *dd, int size) ;

int Zoltan_DD_Set_Neighbor_Hash_Fn2 (Zoltan_DD_Directory *dd, int *proc,
 int *low, int *high, int count) ;

int Zoltan_DD_Set_Neighbor_Hash_Fn3 (Zoltan_DD_Directory *dd, int total) ;

int Zoltan_DD_Print (Zoltan_DD_Directory *dd) ;

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
