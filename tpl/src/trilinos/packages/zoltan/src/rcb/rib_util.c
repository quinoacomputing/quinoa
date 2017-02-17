/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


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
void Zoltan_RIB_Print_Structure(ZZ *zz, int howMany)
{
  int num_obj, i, len;
  RIB_STRUCT *rib;
  struct rib_tree r;
  int printed = 0;

  rib = (RIB_STRUCT *)zz->LB.Data_Structure;
  num_obj = Zoltan_Print_Obj_List(zz, rib->Global_IDs, rib->Local_IDs,
    0, NULL, NULL, howMany);

#if 0
  struct Dot_Struct dot;
  for (i=0; rib->Dots && (i<num_obj); i++){
    dot = rib->Dots[i];
    printf("(Dots %d) (%6.4f %6.4f %6.4f) (%6.4f %6.4f %6.4f %6.4f) proc %d, part %d, new part %dn",
     i, dot.X[0], dot.X[1], dot.X[2],
     dot.Weight[0], dot.Weight[1], dot.Weight[2], dot.Weight[3],
     dot.Proc, dot.Input_Part, dot.Part);
    printed = 1;
  }
  if (!printed){
    printf("Dots: NULL\n");
  }
#endif

  len = zz->LB.Num_Global_Parts;
  printed = 0;

  for (i=0; rib->Tree_Ptr && (i<len); i++){
    r = rib->Tree_Ptr[i];
    printf("(Tree %d) cm %6.4f %6.4f %6.4f, ev %6.4f %6.4f %6.4f, cut: %6.4f, up %d, left %d, right %d\n",
      i, 
      r.cm[0], r.cm[1], r.cm[2],
      r.ev[0], r.ev[1], r.ev[2],
      r.cut, r.parent, r.left_leaf, r.right_leaf);
    printed=1;
  }
  if (!printed){
    printf("Tree: NULL\n");
  }

  printf("Num_Geom: %d\n", rib->Num_Geom);

  Zoltan_Print_Transformation(&(rib->Tran));
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
