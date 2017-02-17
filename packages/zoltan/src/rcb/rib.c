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


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <float.h>
#include "zz_const.h"
#include "rib.h"
#include "params_const.h"
#include "ha_const.h"
#include "par_median_const.h"
#include "rib_params.h"

/* Inertial recursive bisection (RIB) load balancing routine operates on
   "dots" as defined in shared_const.h */

/* Notes:
   dots are balanced across procs by weight (if used)
   on return, proc owns dotnum "dots" in dense array of max-length dotmax
   input weights (if used) are real numbers > 0.0
   can extend "Dot_Struct" data structure in calling program, see shared_const.h
   returned tree only contains one cut on each proc,
   need to do MPI_Allgather if wish to collect it on all procs */

/*  RIB_OUTPUT_LEVEL = 0  No statistics logging */
/*  RIB_OUTPUT_LEVEL = 1  Log times and counts, print summary */
/*  RIB_OUTPUT_LEVEL = 2  Log times and counts, print for each proc */
#define RIB_DEFAULT_OUTPUT_LEVEL 0
#define RIB_DEFAULT_OVERALLOC 1.2


/*---------------------------------------------------------------------------*/
static int rib_fn(ZZ *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **,
                  int *, ZOLTAN_ID_PTR *,
                  double, int, int, int, int, int, float *);
static void print_rib_tree(ZZ *, int, int, struct rib_tree *);
static int compute_rib_direction(ZZ *, int, int, double *, double *,
  struct Dot_Struct *, int *, int, int, double *, double *, double *,
  MPI_Comm, int, int, int);
static int serial_rib(ZZ *, struct Dot_Struct *, int *, int *, int, int,
  int, double, int, int, int *, int *, int, int, int, int, int, int, int,
  struct rib_tree *, double *, double *, float *, double *);

/*---------------------------------------------------------------------------*/

int Zoltan_RIB_Set_Param(
  char *name,                 /* name of variable */
  char *val                   /* value of variable */
)
{
  int status;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */

  status = Zoltan_Check_Param(name, val, RIB_params, &result, &index);

  return(status);
}

/*---------------------------------------------------------------------------*/

int Zoltan_RIB(
  ZZ *zz,                       /* The Zoltan structure with info for
                                the RIB balancer.                       */
  float *part_sizes,            /* Input:  Array of size 
                                zz->Num_Global_Parts * max(zz->Obj_Weight_Dim, 1)
                                containing the percentage of work to be
                                assigned to each part.               */
  int *num_import,              /* Returned value: Number of non-local 
                                objects assigned to
                                this processor in the new decomposition.*/
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value: array of global IDs for
                                non-local objects in this processor's new
                                decomposition.                          */
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value: array of local IDs for
                                non-local objects in this processor's new
                                decomposition.                          */
  int **import_procs,           /* Returned value: array of processor IDs for
                                processors owning the non-local objects in
                                this processor's new decomposition.     */
  int **import_to_part,         /* Returned value:  array of parts to
                                 which imported objects should be assigned. */
  int *num_export,              /* Not computed, set to -1 */
  ZOLTAN_ID_PTR *export_global_ids, /* Not computed. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Not computed. */
  int **export_procs,           /* Not computed. */
  int **export_to_part          /* Not computed. */
)
{
  /* Wrapper routine to set parameter values and call the real rib. */
  double overalloc;           /* amount to overallocate by when realloc
                              of dot array must be done.
                              1.0 = no extra; 1.5 = 50% extra; etc. */
  int wgtflag;                /* No. of weights per dot. */
  int check_geom;             /* Check input & output for consistency? */
  int stats;                  /* Print timing & count summary? */
  int gen_tree;               /* (0) don't (1) generate whole treept to use
                              later for point and box drop. */
  int average_cuts;           /* (0) don't (1) compute the cut to be the
                              average of the closest dots. */
  int idummy;
  int final_output;
  double ddummy;
  int ierr;

  Zoltan_Bind_Param(RIB_params, "RIB_OVERALLOC", (void *) &overalloc);
  Zoltan_Bind_Param(RIB_params, "CHECK_GEOM", (void *) &check_geom);
  Zoltan_Bind_Param(RIB_params, "RIB_OUTPUT_LEVEL", (void *) &stats);
  Zoltan_Bind_Param(RIB_params, "AVERAGE_CUTS", (void *) &average_cuts);
  Zoltan_Bind_Param(RIB_params, "KEEP_CUTS", (void *) &gen_tree);
  Zoltan_Bind_Param(RIB_params, "REDUCE_DIMENSIONS", (void *) &idummy);
  Zoltan_Bind_Param(RIB_params, "DEGENERATE_RATIO", (void *) &ddummy);
  Zoltan_Bind_Param(RIB_params, "FINAL_OUTPUT", (void *) &final_output);

  overalloc = RIB_DEFAULT_OVERALLOC;
  check_geom = DEFAULT_CHECK_GEOM;
  stats = RIB_DEFAULT_OUTPUT_LEVEL;
  gen_tree = 0;
  final_output = 0;
  average_cuts = 0;
  wgtflag = zz->Obj_Weight_Dim;
  idummy = 0;
  ddummy = 0.0;

  Zoltan_Assign_Param_Vals(zz->Params, RIB_params, zz->Debug_Level, zz->Proc,
                    zz->Debug_Proc);

  /* Initializations in case of early exit. */
  *num_import = -1;
  *num_export = -1;  /* We don't compute the export map. */

  if (final_output && (stats < 1)){
    /* FINAL_OUTPUT is a graph/phg param, corresponds to our OUTPUT_LEVEL 1 */
    stats = 1; 
  }

  ierr = rib_fn(zz, num_import, import_global_ids, import_local_ids,
                import_procs, import_to_part, num_export, export_global_ids,
                overalloc, wgtflag, check_geom, stats, gen_tree, average_cuts,
                part_sizes);

  return(ierr);

}

/*---------------------------------------------------------------------------*/

static int rib_fn(
  ZZ *zz,                       /* The Zoltan structure with info for
                                the RIB balancer. */
  int *num_import,              /* Number of non-local objects assigned to
                                this processor in the new decomposition.      
                                When LB.Return_Lists==CANDIDATE_LISTS,
                                num_import returns the number of input
                                objects as given by ZOLTAN_NUM_OBJ_FN. */
  ZOLTAN_ID_PTR *import_global_ids, /* Returned value:  array of global IDs for
                                non-local objects in this processor's new
                                decomposition.  
                                When LB.Return_Lists==CANDIDATE_LISTS,
                                this array contains GIDs for all input 
                                objs as given by ZOLTAN_OBJ_LIST_FN.*/
  ZOLTAN_ID_PTR *import_local_ids,  /* Returned value:  array of local IDs for
                                non-local objects in this processor's new
                                decomposition.                            
                                When LB.Return_Lists==CANDIDATE_LISTS,
                                this array contains LIDs for all input 
                                objs as given by ZOLTAN_OBJ_LIST_FN.*/
  int **import_procs,           /* Returned value: array of processor IDs for
                                processors owning the non-local objects in
                                this processor's new decomposition.       
                                When LB.Return_Lists==CANDIDATE_LISTS,
                                the returned array is NULL. */
  int **import_to_part,         /* Returned value: array of parts to
                                which objects are imported.      
                                When LB.Return_Lists==CANDIDATE_LISTS,
                                the returned array is NULL.  */
  int *num_export,              /* Returned value only when 
                                LB.Return_Lists==CANDIDATE_LISTS; number of
                                input objs as given by ZOLTAN_NUM_OBJ_FN */
  ZOLTAN_ID_PTR *export_global_ids, /* Returned value only when
                                LB.Return_Lists==CANDIDATE_LISTS; for each
                                input obj (from ZOLTAN_OBJ_LIST_FN), 
                                return a candidate obj from the part to which
                                the obj is assigned; used in PHG matching */
  double overalloc,             /* amount to overallocate by when realloc
                                of dot array must be done.
                                  1.0 = no extra; 1.5 = 50% extra; etc. */
  int wgtflag,                  /* No. of weights per dot supplied by user. */
  int check_geom,               /* Check input & output for consistency? */
  int stats,                    /* Print timing & count summary? */
  int gen_tree,                 /* (0) do not (1) do generate full treept */
  int average_cuts,             /* (0) don't (1) compute the cut to be the
                                average of the closest dots. */
  float *part_sizes            /* Input:  Array of size
                                zz->Num_Global_Parts * max(zz->Obj_Weight_Dim, 1)
                                containing the percentage of work to be
                                assigned to each part.               */
)
{
  char    yo[] = "rib_fn";
  int     proc,nprocs;        /* my proc id, total # of procs */
  struct Dot_Struct *dotpt;   /* temporary pointer to local dot arrays */
  int     pdotnum;            /* # of dots - decomposition changes it */
  int    *dotmark = NULL;     /* which side of median for each dot */
  int     dotnum;             /* number of dots */
  int     dotmax = 0;         /* max # of dots arrays can hold */
  int     dottop;             /* dots >= this index are new */
  int     proclower;          /* 1st proc in lower set */
  int     procmid;            /* 1st proc in upper set */
  int     partlower;          /* 1st part in lower set */
  int     partmid;            /* 1st part in upper set */
  int     set;                /* which set processor is in = 0/1 */
  int     old_set;            /* set processor was in last cut = 0/1 */
  int     root;               /* part that stores last cut */
  int     num_procs;          /* number of procs in current set */
  int     num_parts;          /* number of parts in current set */
  int     ierr = ZOLTAN_OK;   /* error flag. */
  double *value = NULL;       /* temp array for median_find */
  double *wgts = NULL;        /* temp array for serial_rib */
  double  valuehalf;          /* median cut position */
  double  cm[3];              /* Center of mass of objects */
  double  evec[3];            /* Eigenvector defining direction */
  int     first_guess = 0;    /* flag if first guess for median search */
  int     allocflag;          /* have to re-allocate space */
  double  time1=0,time2=0;    /* timers */
  double  time3=0,time4=0;    /* timers */
  double  timestart=0,timestop=0; /* timers */
  double  timers[4]={0.,0.,0.,0.}; 
                              /* diagnostic timers
                                 0 = start-up time before recursion
                                 1 = time before median iterations
                                 2 = time in median iterations
                                 3 = communication time */
  ZOLTAN_GNO_TYPE counters[7];  /* diagnostic counts
                                 0 = unused
                                 1 = # of dots sent
                                 2 = # of dots received
                                 3 = most dots this proc ever owns
                                 4 = most dot memory this proc ever allocs
                                 5 = # of times a previous cut is re-used
                                 6 = # of reallocs of dot array */
  int     i, j;               /* local variables */
  int     use_ids;            /* When true, global and local IDs will be
                                 stored along with dots in the RCB_STRUCT.
                                 When false, storage, manipulation, and
                                 communication of IDs is avoided.     
                                 Set by call to Zoltan_RB_Use_IDs().         */

  RIB_STRUCT *rib = NULL;     /* Pointer to data structures for RIB */
  struct rib_tree *treept = NULL; /* tree of cuts - single cut on exit*/

  double start_time, end_time;
  double lb_time[2]={0,0};
  int tfs[2], tmp_tfs[2];     /* added for Tflops_Special; max number
                                 of procs and parts over all processors
                                 in each iteration (while loop) of
                                 parallel partitioning.  */
  int old_nprocs;             /* added for Tflops_Special */
  int old_nparts;             /* added for Tflops_Special */
  double valuelo;             /* smallest value of value[i] */
  double valuehi;             /* largest value of value[i] */
  double weight[RB_MAX_WGTS]; /* weight for current set */
  double weightlo[RB_MAX_WGTS]; /* weight of lower side of cut */
  double weighthi[RB_MAX_WGTS]; /* weight of upper side of cut */
  double fractionlo[RB_MAX_WGTS]; /* desired wt in lower half */
  int *dotlist = NULL;        /* list of dots for find_median.
                                 allocated above find_median for
                                 better efficiency (don't necessarily
                                 have to realloc for each find_median).*/
  int rectilinear_blocks = 0; /* parameter for find_median (not used by rib) */
  int fp=0;                     /* first part assigned to this proc. */
  int np=0;                     /* number of parts assigned to this proc. */
  int wgtdim;                   /* max(wgtflag,1) */
  int *dindx = NULL, *tmpdindx = NULL;

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;
  int free_comm = FALSE;            /* Flag indicating whether MPI_Comm_free
                                       should be called on local_comm at end. */

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    MPI_Barrier(zz->Communicator);
    timestart = time1 = Zoltan_Time(zz->Timer);
  }

  /* setup for parallel */

  proc = zz->Proc;
  nprocs = zz->Num_Proc;
  num_parts = zz->LB.Num_Global_Parts;

  /*
   * Determine whether to store, manipulate, and communicate global and
   * local IDs.
   */
  use_ids = Zoltan_RB_Use_IDs(zz);

  /*
   *  Build the RIB Data structure and
   *  set pointers to information in it.
   */

  start_time = Zoltan_Time(zz->Timer);
  ierr = Zoltan_RIB_Build_Structure(zz, &pdotnum, &dotmax, wgtflag, overalloc,
                                    use_ids, gen_tree);
  if (ierr < 0) {
    ZOLTAN_PRINT_ERROR(proc, yo, 
      "Error returned from Zoltan_RIB_Build_Structure.");
    goto End;
  }

  rib = (RIB_STRUCT *) (zz->LB.Data_Structure);

  treept = rib->Tree_Ptr;
  end_time = Zoltan_Time(zz->Timer);
  lb_time[0] = end_time - start_time;
  start_time = end_time;

  /* local copies of calling parameters */

  dottop = dotnum = pdotnum;

  /* initialize timers and counters */

  counters[0] = 0;
  counters[1] = 0;
  counters[2] = 0;
  counters[3] = dotnum;
  counters[4] = dotmax;
  counters[5] = 0;
  counters[6] = 0;

  /* Ensure there are dots */
    
  MPI_Allreduce(&dotnum, &i, 1, MPI_INT, MPI_MAX, zz->Communicator);
  
  if (i == 0){
    if (proc == 0){
      ZOLTAN_PRINT_WARN(proc, yo, "RIB partitioning called with no objects");
    }
    timestart = timestop = 0;
    goto EndReporting;
  }

  /* If using RIB for matching, need to generate candidate lists.
   * Candidate lists include input GIDs, LIDs as provided by the application.
   * We need to capture that input here before we move any dots!
   * We return it in the import lists.
   * Candidates will be computed after partitioning and returned in the
   * export lists.
   */
  if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
    ierr = Zoltan_RB_Candidates_Copy_Input(zz, dotnum,
                                           rib->Global_IDs, rib->Local_IDs,
                                           &rib->Dots,
                                           num_import,
                                           import_global_ids, import_local_ids,
                                           import_procs, import_to_part);
    if (ierr < 0) {
       ZOLTAN_PRINT_ERROR(proc,yo,
                        "Error returned from Zoltan_RB_Candidates_Copy_Input.");
       goto End;
    }
  }

  /* create mark and list arrays for dots */

  allocflag = 0;
  if (dotmax > 0) {
    if (!(dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))
     || !(value = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double)))
     || !(dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))) {
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }
  else {
    dotmark = NULL;
    value = NULL;
    dotlist = NULL;
  }

  /* set dot weights = 1.0 if user didn't and determine total weight */

  dotpt = &rib->Dots;
  if (dotpt->nWeights == 0) {
    weightlo[0] = (double) dotnum;
    dotpt->uniformWeight = 1.0;
    wgtdim = 1;
  }
  else {
    double *wgt;
    for (j=0; j<dotpt->nWeights; j++){
      weightlo[j] = 0.0;
      wgt = dotpt->Weight + j;
      for (i=0; i < dotnum; i++){
        weightlo[j] += *wgt;
        wgt += dotpt->nWeights;
      }
    }
    wgtdim = dotpt->nWeights;
  }

  MPI_Allreduce(weightlo, weight, wgtdim, MPI_DOUBLE, MPI_SUM, zz->Communicator);

  if (check_geom) {
    ierr = Zoltan_RB_check_geom_input(zz, dotpt, dotnum);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo,
        "Error returned from Zoltan_RB_check_geom_input");
      goto End;
    }
  }

  /* create local communicator for use in recursion */

  if (zz->Tflops_Special)
     local_comm = zz->Communicator;
  else {
     MPI_Comm_dup(zz->Communicator,&local_comm);
     free_comm = TRUE;
  }

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    time2 = Zoltan_Time(zz->Timer);
    timers[0] = time2 - time1;
  }

  /* recursively halve until just one part or proc in set */

  old_nprocs = num_procs = nprocs;
  old_nparts = num_parts;
  partlower = 0;
  root = 0;
  old_set = 1;
  ierr = Zoltan_LB_Proc_To_Part(zz, proc, &np, &fp);
  if (treept) {
    for (i = fp; i < (fp + np); i++) {
      treept[i].parent = 0;
      treept[i].left_leaf = 0;
    }
  }
  if (zz->Tflops_Special) {
    proclower = 0;
    tfs[0] = nprocs;
    tfs[1] = num_parts;
  }

  while ((num_parts > 1 && num_procs > 1) || 
         (zz->Tflops_Special && tfs[0] > 1 && tfs[1] > 1)) {

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time1 = Zoltan_Time(zz->Timer);

    ierr = Zoltan_Divide_Machine(zz, zz->Obj_Weight_Dim, part_sizes, 
                                 proc, local_comm, &set, 
                                 &proclower, &procmid, &num_procs, 
                                 &partlower, &partmid, &num_parts, 
                                 fractionlo);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error in Zoltan_Divide_Machine.");
      goto End;
    }

    /* tfs[0] is max number of processors in all sets over all processors -
     * tfs[1] is max number of parts in all sets over all processors -
     * force all processors to go through all levels of parallel rib */
    if (zz->Tflops_Special) {
      tmp_tfs[0] = num_procs;
      tmp_tfs[1] = num_parts;
      MPI_Allreduce(tmp_tfs, tfs, 2, MPI_INT, MPI_MAX, local_comm);
    }

    /* create mark array and active list for dots */

    if (allocflag) {
      allocflag = 0;
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_FREE(&dotlist);
      if (!(dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))
       || !(value = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double)))
       || !(dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))) {
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }

    dotpt = &rib->Dots;
    
    if (old_nparts > 1 && old_nprocs > 1) { /* test added for Tflops_Special;
                                               compute values only if looping
                                               to decompose, not if looping to
                                               keep Tflops_Special happy.  */
      ierr = compute_rib_direction(zz, zz->Tflops_Special, rib->Num_Geom, 
                                   &valuelo, &valuehi, dotpt, NULL, dotnum, 
                                   wgtflag, cm, evec, value,
                                   local_comm, proc, old_nprocs, proclower);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(proc, yo, 
          "Error returned from compute_rib_direction");
        goto End;
      }
    }
    else {  /* For Tflops_Special: initialize value when looping only 
                                   for Tflops_Special */
      for (i = 0; i < dotmax; i++)
        value[i] = 0.0;
      valuelo = valuehi = 0.0;
    }

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time2 = Zoltan_Time(zz->Timer);

    if (!Zoltan_RB_find_median(
                   zz->Tflops_Special, value, dotpt->Weight, dotpt->uniformWeight, dotmark, dotnum, proc, 
                   fractionlo, local_comm, &valuehalf, first_guess,
                   nprocs, old_nprocs, proclower, old_nparts,
                   wgtflag, valuelo, valuehi, weight[0], weightlo,
                   weighthi, dotlist, rectilinear_blocks, average_cuts)) {
      ZOLTAN_PRINT_ERROR(proc, yo,
        "Error returned from Zoltan_RB_find_median.");
      ierr = ZOLTAN_FATAL;
      goto End;
    }
  
    if (set)    /* set weight for current part */
      for (j=0; j<wgtdim; j++) weight[j] = weighthi[j];
    else
      for (j=0; j<wgtdim; j++) weight[j] = weightlo[j];

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      time3 = Zoltan_Time(zz->Timer);

    /* store cut info in tree only if proc "owns" partmid */
    /* test of partmid > 0 prevents treept[0] being set when this cut is 
       only removing low-numbered processors (proclower to procmid-1) that
       have no parts in them from the processors remaining to 
       be partitioned. */

    if (treept && partmid > 0 && partmid == fp) {
      treept[partmid].cm[0] = cm[0];
      treept[partmid].cm[1] = cm[1];
      treept[partmid].cm[2] = cm[2];
      treept[partmid].ev[0] = evec[0];
      treept[partmid].ev[1] = evec[1];
      treept[partmid].ev[2] = evec[2];
      treept[partmid].cut = valuehalf;
      treept[partmid].parent = old_set ? -(root+1) : root+1;
      /* The following two will get overwritten when the information
         is assembled if this is not a terminal cut */
      treept[partmid].left_leaf = -partlower;
      treept[partmid].right_leaf = -partmid;
    }

    if (old_nprocs > 1 && partmid > 0 && partmid != partlower + old_nparts) {  
      /* old_nprocs > 1 test: Don't reset these values if proc is in loop only 
       * because of other procs for Tflops_Special.
       * partmid > 0 test:  Don't reset these values if low-numbered processors
       * (proclower to procmid-1) have zero parts and this cut is removing
       * them from the processors remaining to be partitioned. 
       * partmid != partlower + old_nparts test:  Don't reset these values if
       * cut is removing high-numbered processors with zero parts from
       * the processors remaining to be partitioned.
       */
      old_set = set;
      root = partmid;
    }

    ierr = Zoltan_RB_Send_Outgoing(zz, &(rib->Global_IDs), &(rib->Local_IDs), 
                               &(rib->Dots), &dotmark,
                               &dottop, &dotnum, &dotmax,
                               set, &allocflag, overalloc,
                               stats, counters, use_ids,  
                                local_comm, proclower,
                               old_nprocs, partlower, partmid);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo,
        "Error returned from Zoltan_RB_Send_Outgoing.");
      goto End;
    }
    
    /* create new communicators */

    if (zz->Tflops_Special) {
      if (set) {
        proclower = procmid;
        partlower = partmid;
      }
      old_nprocs = num_procs;
      old_nparts = num_parts;
    }
    else {
      if (set) partlower = partmid;
      MPI_Comm_split(local_comm,set,proc,&tmp_comm);
      MPI_Comm_free(&local_comm);
      local_comm = tmp_comm;
      old_nprocs = num_procs;
      old_nparts = num_parts;
    }

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
      time4 = Zoltan_Time(zz->Timer);
      timers[1] += time2 - time1;
      timers[2] += time3 - time2;
      timers[3] += time4 - time3;
    }
  }

  /* have recursed all the way to a single processor sub-domain */

  /* Send dots to correct processors for their parts.  This is needed
     most notably when a processor has zero parts on it, but still has
     some dots after the parallel partitioning. */

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
    time1 = Zoltan_Time(zz->Timer);

  ierr = Zoltan_RB_Send_To_Part(zz, &(rib->Global_IDs), &(rib->Local_IDs),
                               &(rib->Dots), &dotmark, &dottop,
                               &dotnum, &dotmax, &allocflag, overalloc,
                               stats, counters, use_ids);

  if (ierr < 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Zoltan_RB_Send_To_Part");
    goto End;
  }

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
    timers[3] += (Zoltan_Time(zz->Timer) - time1);
  

  /* All dots are now on the processors they will end up on; now generate
   * more parts if needed. */

  if (num_parts > 1) {
    if (dotpt->nWeights) 
      wgts = (double *) ZOLTAN_MALLOC(dotpt->nWeights * dotnum * sizeof(double));
    dindx = (int *) ZOLTAN_MALLOC(dotnum * 2 * sizeof(int));
    tmpdindx = dindx + dotnum;
    if (allocflag) {
      ZOLTAN_FREE(&dotmark);
      ZOLTAN_FREE(&value);
      ZOLTAN_FREE(&dotlist);
      if (!(dotmark = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))
       || !(value = (double *) ZOLTAN_MALLOC(dotmax*sizeof(double)))
       || !(dotlist = (int *) ZOLTAN_MALLOC(dotmax*sizeof(int)))) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }
    for (i = 0; i < dotnum; i++)
      dindx[i] = i;

    ierr = serial_rib(zz, &rib->Dots, dotmark, dotlist, old_set, root,
                      rib->Num_Geom, weight[0], dotnum, num_parts,
                      &(dindx[0]), &(tmpdindx[0]), partlower,
                      proc, wgtflag, stats, gen_tree,
                      rectilinear_blocks, average_cuts,
                      treept, value, wgts, part_sizes, timers);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo, "Error returned from serial_rib");
      goto End;
    }
    ZOLTAN_FREE(&wgts);
  }

  end_time = Zoltan_Time(zz->Timer);
  lb_time[1] = end_time - start_time;

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
    MPI_Barrier(zz->Communicator);
    timestop = time1 = Zoltan_Time(zz->Timer);
  }

  /* error checking and statistics */

  if (check_geom) {
    ierr = Zoltan_RB_check_geom_output(zz, &rib->Dots, part_sizes, np, fp,
                                       dotnum, pdotnum, NULL);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo,
                         "Error returned from Zoltan_RB_check_geom_output");
      goto End;
    }
  }

EndReporting:

  /* update calling routine parameters */

  start_time = Zoltan_Time(zz->Timer);

  pdotnum = dotnum;

  /* Perform remapping (if requested) */

  if (zz->LB.Remap_Flag) {
    ierr = Zoltan_RB_Remap(zz, &(rib->Global_IDs), &(rib->Local_IDs),
                               &(rib->Dots), &dotnum, &dotmax,
                  &allocflag, overalloc, stats, counters, use_ids);
    /* Note:  dottop is no longer valid after remapping.  Remapping might
       destroy the nice local-followed-by-non-local ordering of the 
       dots array.  Do not use dottop after remapping. */
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_RB_Remap.");
      goto End;
    }
  }

  /*  build return arguments */

  if (zz->LB.Return_Lists != ZOLTAN_LB_NO_LISTS &&
      zz->LB.Return_Lists != ZOLTAN_LB_CANDIDATE_LISTS) {
    /* zz->LB.Return_Lists is true ==> use_ids is true */
    ierr = Zoltan_RB_Return_Arguments(zz, rib->Global_IDs, rib->Local_IDs, 
                                      &rib->Dots, num_import,
                                      import_global_ids, import_local_ids,
                                      import_procs, import_to_part, 
                                      dotnum);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo,
                     "Error returned from Zoltan_RB_Return_Arguments.");
      goto End;
    }
  }
  else if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
    /* Select a candidate for each part and return it in the export_GIDs. */
    ierr = Zoltan_RB_Candidates_Output(zz, dotnum, dindx,
                                       rib->Global_IDs, rib->Local_IDs,
                                       &rib->Dots,
                                       *num_import, *import_global_ids,
                                       num_export, export_global_ids);
    if (ierr < 0) {
       ZOLTAN_PRINT_ERROR(proc,yo,
                          "Error returned from Zoltan_RB_Return_Candidates.");
       goto End;
    }
  }
  ZOLTAN_FREE(&dindx);


  if (gen_tree) {
    int *displ, *recvcount;
    int sendcount;
    struct rib_tree *treetmp = NULL; /* temporary tree of cuts; used to keep
					valgrind from reporting overlapped
					memory in MPI_Allgatherv */

    treetmp = (struct rib_tree *)
	       ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(struct rib_tree));
    displ = (int *) ZOLTAN_MALLOC(2 * zz->Num_Proc * sizeof(int));
    if (!displ || !treetmp) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      ZOLTAN_FREE(&displ);
      ZOLTAN_FREE(&treetmp);
      goto End;
    }
    recvcount = displ + zz->Num_Proc;

    ierr = Zoltan_RB_Tree_Gatherv(zz, sizeof(struct rib_tree), &sendcount,
                                  recvcount, displ);
    /* 
     * Create copy of treept so that MPI_Allgatherv doesn't use same
     * memory for sending and receiving; removes valgrind warning.
     */
    for (i = 0; i < zz->LB.Num_Global_Parts; i++)
      treetmp[i] = treept[i];

    MPI_Allgatherv(&treetmp[fp], sendcount, MPI_BYTE, treept, recvcount, displ,
                   MPI_BYTE, zz->Communicator);
    for (i = 1; i < zz->LB.Num_Global_Parts; i++){
      if (treept[i].parent > 0)
        treept[treept[i].parent - 1].left_leaf = i;
      else if (treept[i].parent < 0)
        treept[-treept[i].parent - 1].right_leaf = i;
    }

    ZOLTAN_FREE(&displ);
    ZOLTAN_FREE(&treetmp);
  }

  if (treept && zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    print_rib_tree(zz, np, fp, &(treept[fp]));

  end_time = Zoltan_Time(zz->Timer);
  lb_time[0] += (end_time - start_time);

  if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME))
    Zoltan_RB_stats(zz, timestop-timestart, &rib->Dots, dotnum, 
                part_sizes, timers, counters, stats, NULL, NULL, FALSE);

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
    if (zz->Proc == zz->Debug_Proc)
      printf("ZOLTAN RIB Times:  \n");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, lb_time[0], 
                   "ZOLTAN       Build:       ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, lb_time[1], 
                   "ZOLTAN         RIB:         ");
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    /* zz->Debug_Level >= ZOLTAN_DEBUG_ALL ==> use_ids is true */
    Zoltan_RB_Print_All(zz, rib->Global_IDs, &rib->Dots, 
                    dotnum, *num_import, 
                    *import_global_ids, *import_procs);
  }

End:

  /* Free memory allocated by the algorithm.  */

  if (free_comm) MPI_Comm_free(&local_comm);
  ZOLTAN_FREE(&dotmark);
  ZOLTAN_FREE(&value);
  ZOLTAN_FREE(&dotlist);

  if (!gen_tree &&                         /* don't need parts */
      rib && (rib->Tran.Target_Dim < 0)) { /* don't need transformation */
    /* Free all memory used. */
    Zoltan_RIB_Free_Structure(zz);
  }
  else if (rib != NULL) {
    /* Free only Dots and IDs; keep other structures. */
    ZOLTAN_FREE(&(rib->Global_IDs));
    ZOLTAN_FREE(&(rib->Local_IDs));
    Zoltan_Free_And_Reset_Dot_Structure(&(rib->Dots));
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

/*****************************************************************************/

static void print_rib_tree(ZZ *zz, int np, int fp, struct rib_tree *treept_arr)
{
int i;
struct rib_tree *treept = treept_arr;

  Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
  for (i = fp; i < fp+np; i++) {
    printf("Proc %d Part %d:  Tree Struct:\n", zz->Proc, i);
    printf("          cm         = (%e,%e,%e)\n", 
           treept->cm[0], treept->cm[1], treept->cm[2]);
    printf("          ev         = (%e,%e,%e)\n", 
           treept->ev[0], treept->ev[1], treept->ev[2]);
    printf("          cut        = %e\n", treept->cut);
    printf("          parent     = %d\n", treept->parent);
    printf("          left_leaf  = %d\n", treept->left_leaf);
    printf("          right_leaf = %d\n", treept->right_leaf);
    treept++;
  }
  Zoltan_Print_Sync_End(zz->Communicator, TRUE);
}

/*****************************************************************************/

void Zoltan_RIB_min_max(
   double   *min,             /* minimum value */
   double   *max,             /* maximum value */
   int      proclower,        /* smallest processor in partition */
   int      proc,             /* rank of processor in global partition */
   int      nprocs,           /* number of processors in partition */
   MPI_Comm comm
)
{
   double   tmp[2], tmp1[2];  /* temporaries for min/max */
   int      tag = 32100;      /* message tag */
   int      rank;             /* rank of processor in partition */
   int      partner;          /* message partner in binary exchange */
   int      to;               /* message partner not in binary exchange */
   int      mask;             /* mask to determine communication partner */
   int      nprocs_small;     /* largest power of 2 contained in nprocs */
   int      hbit;             /* 2^hbit = nproc_small */
   MPI_Status status;

   /* This routine finds the global min of min and the global max of max */

   rank = proc - proclower;
   tmp[0] = *min;
   tmp[1] = *max;

   /* Find next lower power of 2. */
   for (hbit = 0; (nprocs >> hbit) != 1; hbit++);
 
   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nprocs) {
      nprocs_small *= 2;
      hbit++;
   }
 
   to = proclower + (rank ^ nprocs_small);
   if (rank & nprocs_small) {  /* processors greater than largest power of 2 */
      MPI_Send(tmp, 2, MPI_DOUBLE, to, tag, comm);
      tag += hbit + 1;
      MPI_Recv(tmp, 2, MPI_DOUBLE, to, tag, comm, &status);
   }
   else {   /* processors within greatest power of 2 */
      if (rank + nprocs_small < nprocs) {
         MPI_Recv(tmp1, 2, MPI_DOUBLE, to, tag, comm, &status);
         if (tmp1[0] < tmp[0]) tmp[0] = tmp1[0];
         if (tmp1[1] > tmp[1]) tmp[1] = tmp1[1];
      }  
      for (mask = nprocs_small >> 1; mask; mask >>= 1) { /* binary exchange */
         tag++;
         partner = proclower + (rank ^ mask);
         MPI_Send(tmp, 2, MPI_DOUBLE, partner, tag, comm);
         MPI_Recv(tmp1, 2, MPI_DOUBLE, partner, tag, comm, &status);
         if (tmp1[0] < tmp[0]) tmp[0] = tmp1[0];
         if (tmp1[1] > tmp[1]) tmp[1] = tmp1[1];
      }  
      tag++;
      if (rank + nprocs_small < nprocs)
         MPI_Send(tmp, 2, MPI_DOUBLE, to, tag, comm);
   }

   *min = tmp[0];
   *max = tmp[1];
}

/*****************************************************************************/

static int compute_rib_direction(
  ZZ *zz, 
  int Tflops_Special,         /* Tflops_Special flag for special processing.
                                 Should be 0 when called by serial_rib. */
  int num_geom,               /* number of dimensions */
  double *valuelo,            /* smallest value of value[i] */
  double *valuehi,            /* largest value of value[i] */
  struct Dot_Struct *dotpt,   /* local dot array */
  int *dindx,                 /* index array into dotpt; if NULL, access dotpt
                                 directly */
  int dotnum,                 /* number of dots */
  int wgtflag,                /* (0) do not (1) do use weights.
                                 Multidimensional weights not supported */
  double *cm,                 /* Center of mass of objects */
  double *evec,               /* Eigenvector defining direction */
  double *value,              /* temp array for median_find; rotated coords */
  MPI_Comm local_comm,        /* MPI communicator for set */
  int proc,                   /* Current processor; needed for Tflops_Special */
  int nprocs,                 /* Number of procs in operation; needed for
                                 Tflops_Special */
  int proclower               /* Lowest numbered proc; needed for 
                                 Tflops_Special */
)
{
int i, ierr = 0;
double tmp;
RIB_STRUCT *rib;

  rib = (RIB_STRUCT *)zz->LB.Data_Structure;

  if (rib->Tran.Target_Dim > 0){ 
    num_geom = rib->Tran.Target_Dim; /* degenerate geometry */
  }

  switch (num_geom) {
  case 3:
    ierr = Zoltan_RIB_inertial3d(Tflops_Special, dotpt, dindx, dotnum, wgtflag, 
                                 cm, evec, value,
                                 local_comm, proc, nprocs, proclower);
    break;
  case 2:
    ierr = Zoltan_RIB_inertial2d(Tflops_Special, dotpt, dindx, dotnum, wgtflag, 
                                 cm, evec, value,
                                 local_comm, proc, nprocs, proclower);
    break;
  case 1:
    ierr = Zoltan_RIB_inertial1d(dotpt, dindx, dotnum, wgtflag,
                                 cm, evec, value);
    break;
  }
  *valuelo = DBL_MAX;
  *valuehi = -DBL_MAX;
  for (i = 0; i < dotnum; i++) {
    if (value[i] < *valuelo) *valuelo = value[i];
    if (value[i] > *valuehi) *valuehi = value[i];
  }
  if (Tflops_Special)
    Zoltan_RIB_min_max(valuelo, valuehi, proclower, proc, nprocs,
                       local_comm);
  else {
    tmp = *valuehi;
    MPI_Allreduce(&tmp, valuehi, 1, MPI_DOUBLE, MPI_MAX, local_comm);
    tmp = *valuelo;
    MPI_Allreduce(&tmp, valuelo, 1, MPI_DOUBLE, MPI_MIN, local_comm);
  }

  return ierr;
}

/*****************************************************************************/

static int serial_rib(
  ZZ *zz, 
  struct Dot_Struct *dotpt,  /* local dot array */
  int *dotmark,              /* which side of median for each dot */
  int *dotlist,              /* list of dots used only in find_median;
                                allocated above find_median for
                                better efficiency (don't necessarily
                                have to realloc for each find_median).*/
  int old_set,               /* Set the objects to be partitioned were in
                                for last cut */
  int root,                  /* part for which last cut was stored. */
  int num_geom,              /* number of dimensions */
  double weight,             /* Weight for current set */
  int dotnum,                /* number of dots */
  int num_parts,             /* number of parts to create. */
  int *dindx,                /* Index into dotpt for dotnum dots to be
                                partitioned; reordered in serial_rib so set0
                                dots are followed by set1 dots. */
  int *tmpdindx,             /* Temporary memory used in reordering dindx. */
  int partlower,             /* smallest part number to be created. */
  int proc,                  /* processor number. */
  int wgtflag,               /* No. of weights per dot provided by user. */
  int stats,                 /* Print timing & count summary?             */
  int gen_tree,              /* (0) do not (1) do generate full treept    */
  int rectilinear_blocks,    /* parameter for find_median (not used by rib) */
  int average_cuts,          /* (0) don't (1) compute the cut to be the
                                average of the closest dots. */
  struct rib_tree *treept,   /* tree of RCB cuts */
  double *value,             /* temp array for median_find */
  double *wgts,              /* temp array for serial_rib */
  float *part_sizes,         /* Array of size zz->LB.Num_Global_Parts
                                containing the percentage of work to be
                                assigned to each part.               */
  double timers[]            /* as in rib_fn */
)
{
char *yo = "serial_rib";
int ierr = ZOLTAN_OK;
double valuelo;            /* smallest value of value[i] */
double valuehi;            /* largest value of value[i] */
int partmid;
int new_nparts;
double fractionlo[RB_MAX_WGTS];
double valuehalf;
double weightlo, weighthi;
double cm[3];                 /* Center of mass of objects */
double evec[3];               /* Eigenvector defining direction */
int set0, set1;
int i;
double start_time=0., end_time;

  if (num_parts == 1) {
    for (i = 0; i < dotnum; i++)
      dotpt->Part[dindx[i]] = partlower;
  }
  else {

    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) 
      start_time = Zoltan_Time(zz->Timer);

    ierr = Zoltan_Divide_Parts(zz, zz->Obj_Weight_Dim, part_sizes, num_parts,
                               &partlower, &partmid, fractionlo);

    ierr = compute_rib_direction(zz, 0, num_geom, &valuelo, &valuehi, 
                                 dotpt, dindx, dotnum, wgtflag, cm, evec, value,
                                 MPI_COMM_SELF, proc, 1, proc);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(proc, yo, 
        "Error returned from compute_rib_direction");
      goto End;
    }

    if (dotpt->nWeights){
      double *w = wgts;
      int j;
      for (i = 0; i < dotnum; i++) {
        for (j=0; j<dotpt->nWeights; j++){
          *w++ = dotpt->Weight[dindx[i] * dotpt->nWeights + j];
        }
      }
    }
    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
      end_time = Zoltan_Time(zz->Timer);
      timers[1] += end_time - start_time;
      start_time = end_time;
    }

    if (!Zoltan_RB_find_median(0, value, wgts, dotpt->uniformWeight, dotmark, dotnum, proc, 
                               fractionlo, MPI_COMM_SELF, &valuehalf, 
                               0, zz->Num_Proc, 1, proc, num_parts,
                               wgtflag, valuelo, valuehi, weight, &weightlo,
                               &weighthi, dotlist, rectilinear_blocks, 
                               average_cuts)) {
      ZOLTAN_PRINT_ERROR(proc, yo, 
        "Error returned from Zoltan_RB_find_median.");
      ierr = ZOLTAN_FATAL;
      goto End;
    }
    if (treept) {
      treept[partmid].cm[0] = cm[0];
      treept[partmid].cm[1] = cm[1];
      treept[partmid].cm[2] = cm[2];
      treept[partmid].ev[0] = evec[0];
      treept[partmid].ev[1] = evec[1];
      treept[partmid].ev[2] = evec[2];
      treept[partmid].cut = valuehalf;
      treept[partmid].parent = old_set ? -(root+1) : root+1;
      /* The following two will get overwritten when the information
         is assembled if this is not a terminal cut */
      treept[partmid].left_leaf = -partlower;
      treept[partmid].right_leaf = -partmid;
    }

    root = partmid;

      /* Create new dindx, grouping set 0 and set 1 dots together */
    for (set0 = 0, set1 = dotnum, i = 0; i < dotnum; i++) {
      if (dotmark[i] == 0)
        tmpdindx[set0++] = dindx[i];
      else
        tmpdindx[--set1] = dindx[i];
    }
    memcpy(dindx, tmpdindx, dotnum * sizeof(int));
    if (stats || (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME)) {
      end_time = Zoltan_Time(zz->Timer);
      timers[2] += (end_time - start_time);
    }

    /* If set 0 has at least one part and at least one dot,
     * call serial_rib for set 0 */
    new_nparts = partmid - partlower;
    if (new_nparts > 0 && set1 != 0) {
      ierr = serial_rib(zz, dotpt, dotmark, dotlist, 0, root, num_geom,
                        weightlo, set0, new_nparts,
                        &(dindx[0]), &(tmpdindx[0]), partlower,
                        proc, wgtflag, stats, gen_tree, 
                        rectilinear_blocks, average_cuts,
                        treept, value, wgts, part_sizes, timers);
      if (ierr < 0) {
        goto End;
      }
    }

    /* If set 1 has at least one part and at least one dot,
     * call serial_rcb for set 1 */
    new_nparts = partlower + num_parts - partmid;
    if (new_nparts > 0 && set0 != dotnum) {
      ierr = serial_rib(zz, dotpt, dotmark, dotlist, 1, root, num_geom,
                        weighthi, dotnum-set0, new_nparts,
                        &(dindx[set1]), &(tmpdindx[set1]), partmid,
                        proc, wgtflag, stats, gen_tree,
                        rectilinear_blocks, average_cuts,
                        treept, value, wgts, part_sizes, timers);
      if (ierr < 0) {
        goto End;
      }
    }
  }
    
End:
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
