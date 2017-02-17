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
#include "dr_util_const.h"
#include "dr_err_const.h"
#include "dr_dd.h"

#include "zoltan_dd_cpp.h"

static Zoltan_DD *dd = NULL;

/****************************************************************************/
int build_elem_dd(MESH_INFO_PTR mesh) 
{
  destroy_elem_dd();

  dd = new Zoltan_DD(MPI_COMM_WORLD, 1, 0, 0, 0, 0);

  return update_elem_dd(mesh);
}
/****************************************************************************/
void destroy_elem_dd()
{
  if (dd)
    {
    delete dd;
    dd = NULL;
    }
}
/****************************************************************************/
int update_elem_dd(MESH_INFO_PTR mesh)
{
/*  Update a distributed directory of the elements with processor
 *  assignments initially and after migration.
 */

  ZOLTAN_ID_PTR gids = new ZOLTAN_ID_TYPE [mesh->num_elems];
  int *parts = new int [mesh->num_elems];

  int i, j;

  for (j = 0, i = 0; i < mesh->elem_array_len; i++) {
    ELEM_INFO_PTR current_elem = &(mesh->elements[i]);
    if (current_elem->globalID != ZOLTAN_ID_INVALID) {
      gids[j] = current_elem->globalID;
      parts[j] = current_elem->my_part;
      j++;
    }
  }

  int rc = 1;

  if (dd->Update(gids, NULL,NULL, parts, mesh->num_elems)!=0) {
    Gen_Error(0, "fatal:  NULL returned from Zoltan_DD::Update()\n");
    rc = 0;
  }

  delete [] gids;
  delete [] parts;
  
  return rc;
}

/****************************************************************************/
int update_hvertex_proc(MESH_INFO_PTR mesh)
{
  if (dd->Find(mesh->hvertex, NULL, NULL, NULL, 
                     mesh->hindex[mesh->nhedges], mesh->hvertex_proc) != 0) {
    Gen_Error(0, "fatal:  NULL returned from Zoltan_DD::Find()\n");
    return 0;
  }

  return 1;
}
