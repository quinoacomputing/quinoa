/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
#include "zz_const.h"
#include "rib.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RIB_Build_Structure(ZZ *zz, int *num_obj, int *max_obj, int wgtflag,
                               double overalloc, int use_ids, int gen_tree)
{
/* Function to build the geometry-based data structures for RIB method. */
char           *yo = "Zoltan_RIB_Build_Structure";
RIB_STRUCT     *rib;                  /* Data structure for RIB.             */
struct rib_tree *treeptr;
int            i, ierr = 0;

  /* Allocate an RIB data structure for this Zoltan structure.
     If the previous data structure is still there, free the Dots and IDs first;
     the other fields can be reused. */

  if (zz->LB.Data_Structure == NULL) {
    rib = (RIB_STRUCT *) ZOLTAN_MALLOC(sizeof(RIB_STRUCT));
    if (rib == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    zz->LB.Data_Structure = (void *) rib;
    rib->Tree_Ptr = NULL;
    rib->Global_IDs = NULL;
    rib->Local_IDs = NULL;

    Zoltan_Initialize_Transformation(&(rib->Tran));

    if (gen_tree) {
      rib->Tree_Ptr = (struct rib_tree *)
                ZOLTAN_CALLOC(zz->LB.Num_Global_Parts, sizeof(struct rib_tree));
      if (rib->Tree_Ptr == NULL) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        Zoltan_RIB_Free_Structure(zz);
        return(ZOLTAN_MEMERR);
      }
      /* initialize Tree_Ptr */
      for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
        treeptr = &(rib->Tree_Ptr[i]);
        treeptr->cm[0] = treeptr->cm[1] = treeptr->cm[2] = 0.0;
        treeptr->ev[0] = treeptr->ev[1] = treeptr->ev[2] = 0.0;
        treeptr->cut = 0.0;
        treeptr->parent = treeptr->left_leaf = treeptr->right_leaf = 0;
      }
    }
  }
  else {
    rib = (RIB_STRUCT *) zz->LB.Data_Structure;
    ZOLTAN_FREE(&(rib->Global_IDs));
    ZOLTAN_FREE(&(rib->Local_IDs));
    ZOLTAN_FREE(&(rib->Dots));
  }

  ierr = Zoltan_RB_Build_Structure(zz, &(rib->Global_IDs), &(rib->Local_IDs),
                               &(rib->Dots), num_obj, max_obj, &(rib->Num_Geom),
                               wgtflag, overalloc, use_ids, 0);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Zoltan_RB_Build_Structure.");
    Zoltan_RIB_Free_Structure(zz);
    return(ierr);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_RIB_Free_Structure(ZZ *zz)
{
/* Deallocate the persistent RIB data structures in zz->Structure.  */
RIB_STRUCT    *rib;                   /* Data structure for RIB. */

  rib = (RIB_STRUCT *) (zz->LB.Data_Structure);

  if (rib != NULL) {
    ZOLTAN_FREE(&(rib->Tree_Ptr));
    ZOLTAN_FREE(&(rib->Global_IDs));
    ZOLTAN_FREE(&(rib->Local_IDs));
    Zoltan_Free_And_Reset_Dot_Structure(&rib->Dots);
    ZOLTAN_FREE(&(zz->LB.Data_Structure));
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#define COPY_BUFFER(buf, type, num) \
  if (from->buf) { \
    to->buf = (type *)ZOLTAN_MALLOC((num) * sizeof(type)); \
    if (!to->buf) { \
      Zoltan_RIB_Free_Structure(toZZ); \
      ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory."); \
      return(ZOLTAN_MEMERR); \
    } \
    memcpy(to->buf, from->buf, (num) * sizeof(type)); \
  } \
  else { \
    to->buf = NULL; \
  }

int Zoltan_RIB_Copy_Structure(ZZ *toZZ, ZZ const *fromZZ)
{
  char *yo = "Zoltan_RIB_Copy_Structure";
  RIB_STRUCT *to;
  RIB_STRUCT const *from;

  from = (RIB_STRUCT const *)fromZZ->LB.Data_Structure;
  Zoltan_RIB_Free_Structure(toZZ);

  if (!from){
    return ZOLTAN_OK;
  }

  to = (RIB_STRUCT *)ZOLTAN_MALLOC(sizeof(RIB_STRUCT));
  if (to == NULL) {
    ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  toZZ->LB.Data_Structure = (void *)to;

  *to = *from;

  COPY_BUFFER(Tree_Ptr, struct rib_tree, fromZZ->LB.Num_Global_Parts);

  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
