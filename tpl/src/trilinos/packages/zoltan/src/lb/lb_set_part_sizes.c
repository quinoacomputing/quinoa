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


#include "zz_const.h"
#include "zz_util_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines to set the part sizes.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Set_Part_Sizes(ZZ *zz, int global_num,
    int len, int *part_ids, int *wgt_idx, float *part_sizes)
{
/*
 *  Function to set the desired part sizes. This function
 *  only sets values locally. Later, Zoltan_LB_Get_Part_Sizes
 *  collects all the information across processors.
 *
 *  Input:
 *    zz            --  The Zoltan structure to which this method
 *                      applies.
 *    global_num    --  Global part numbers? (0 for local numbers)
 *    len           --  Length of arrays wgt_idx, part_idx, part_sizes
 *    part_ids      --  Array of part ids (local or global)
 *    wgt_idx       --  Array of indices between 0 and Obj_Wgt_Dim-1
 *    part_sizes    --  Array of floats that gives the desired part 
 *                      size for each weight and each part, i.e., 
 *                      part_sizes[i] corresponds to wgt_idx[i] and part_id[i]
 *
 *  Output:
 *    zz->LB.*      --  Appropriate fields set to designated values.
 *    Return value  --  Error code.
 */

  char *yo = "Zoltan_LB_Set_Part_Sizes";
  int i, j, maxlen=0;
  int error = ZOLTAN_OK;
  const int INIT_NUM_PART = 64; /* Initial allocation for Part_Info array. */

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* len = -1 will nullify all part sizes set on this proc */
  if (len == -1){
    zz->LB.Part_Info_Len = 0;
    goto End;
  }

  /* Verify input. */
  if ((part_ids==NULL) || (part_sizes==NULL)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid input argument NULL.");
    error = ZOLTAN_FATAL;
    goto End;
  }

  /* Do we need more space? */
  if ((!zz->LB.Part_Info) || (zz->LB.Part_Info_Max_Len==0)){
    maxlen = INIT_NUM_PART;          /* Start with space for INIT_NUM_PART */
    zz->LB.Part_Info = (struct Zoltan_part_info *) ZOLTAN_MALLOC(maxlen *
      sizeof(struct Zoltan_part_info));
  }
  if (zz->LB.Part_Info_Len + len > zz->LB.Part_Info_Max_Len){
    maxlen = 2*(zz->LB.Part_Info_Len + len);  /* Double the length */
    zz->LB.Part_Info = (struct Zoltan_part_info *) ZOLTAN_REALLOC( zz->LB.Part_Info, 
                         maxlen * sizeof(struct Zoltan_part_info));
  }

  if (zz->LB.Part_Info == NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      error = ZOLTAN_MEMERR;
      goto End;
  }

  /* Add new data to part info array. */
  for (i=0,j=zz->LB.Part_Info_Len; i<len; i++,j++){
    zz->LB.Part_Info[j].Size = part_sizes[i];
    zz->LB.Part_Info[j].Part_id = part_ids[i]; 
    zz->LB.Part_Info[j].Idx = (wgt_idx ? wgt_idx[i] : 0); 
    zz->LB.Part_Info[j].Global_num = global_num;
  }

  /* Update values in LB. */
  zz->LB.Part_Info_Len += len;
  if (maxlen > zz->LB.Part_Info_Max_Len)
    zz->LB.Part_Info_Max_Len = maxlen;

End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return error;
}

/***************************************************************************/

int Zoltan_LB_Get_Part_Sizes(ZZ *zz, int part_dim, float *part_sizes)
{
/*
 *  Function to get the scaled part sizes.
 *
 *  Input:
 *    zz            --  The Zoltan structure to which this method
 *                      applies.
 *    part_dim      --  The number of object weights per part.
 *                      (This usually equals lb->Obj_Wgt_Dim.)
 *
 *  Output:
 *    part_sizes    --  Array of floats that gives the set part 
 *                      sizes, scaled such that they sum to one.
 */
  int i, j, nparts, fpart;
  int num_global_parts = zz->LB.Num_Global_Parts;
  float *temp_part_sizes=NULL, *sum=NULL;
  int error = ZOLTAN_OK;
  char msg[128];
  static char *yo = "Zoltan_LB_Get_Part_Sizes";

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("[%1d] Debug: num_global_parts = %d\n", zz->Proc, num_global_parts);

  /* Barrier to make sure all procs have finished Zoltan_LB_Set_Part_Sizes */
  MPI_Barrier(zz->Communicator);

  /* For convenience, if no weights are used, set part_dim to 1 */
  if (part_dim==0) part_dim = 1;

  if (part_sizes == NULL){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input argument part_sizes is NULL.");
    error = ZOLTAN_FATAL;
    goto End;
  }

  /* Find max Part_Info_Len over all procs to see if they are all zero. */
  MPI_Allreduce((void*) &(zz->LB.Part_Info_Len), (void*) &j, 
      1, MPI_INT, MPI_MAX, zz->Communicator);

  if (j == 0){
    /* Uniform part sizes. */
    zz->LB.Uniform_Parts = 1;
    for (i = 0; i < num_global_parts*part_dim; i++)
      part_sizes[i] = 1.0 / (float)num_global_parts;
  }
  else {
   /* Get the part sizes set by the user (application).
    * Each processor puts its data in a part_dim * num_global_parts
    * array. Then we gather all the data across processors.
    * Out-of-range part size data is ignored.
    */
    zz->LB.Uniform_Parts = 0;

    /* Pack LB.Part_Info into temp array */
    temp_part_sizes = (float *)ZOLTAN_MALLOC(num_global_parts*part_dim
                                             *sizeof(float));
    sum = (float *)ZOLTAN_MALLOC(part_dim*sizeof(float));
    if ((!temp_part_sizes) || (!sum)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      error = ZOLTAN_MEMERR;
      goto End;
    }
    for (i = 0; i < num_global_parts*part_dim; i++){
      temp_part_sizes[i] = -1.0;
    }
    for (i = 0; i < zz->LB.Part_Info_Len; i++){
      /* Only assemble part sizes for parts and weights
         in the requested range. */
      if (zz->LB.Part_Info[i].Idx < part_dim){
        j = zz->LB.Part_Info[i].Part_id;
        if (zz->LB.Part_Info[i].Global_num == 0) {
          Zoltan_LB_Proc_To_Part(zz, zz->Proc, &nparts, &fpart);
          j += fpart;
        }
        if (j >= num_global_parts){
          sprintf(msg, "Part number %d is >= num_global_parts %d.",
            j, num_global_parts);
          ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
          error = ZOLTAN_WARN;
        }
        else
          temp_part_sizes[j*part_dim + zz->LB.Part_Info[i].Idx] 
            = zz->LB.Part_Info[i].Size;
      }
    }

    /* Reduce over all procs */
    MPI_Allreduce((void*) temp_part_sizes, (void*) part_sizes, 
      num_global_parts*part_dim, MPI_FLOAT, MPI_MAX, zz->Communicator);
  
    /* Check for errors. Scale the sizes so they sum to one for each weight. */
    for (j = 0; j < part_dim; j++) 
      sum[j] = 0.0;

    for (i = 0; i < num_global_parts; i++){
      for (j = 0; j < part_dim; j++){
        if (part_sizes[i*part_dim+j]<0)
          part_sizes[i*part_dim+j] = 1.0; /* default value if not set */
        sum[j] += part_sizes[i*part_dim+j];
      }

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] In %s: Part size %1d (before scaling) = ",  
            zz->Proc, yo, i);
        for (j = 0; j < part_dim; j++)
          printf("%f, ",  part_sizes[i*part_dim+j]);
        printf("\n");
      }
    }

    /* Check for sum[j] == 0 (error). */
    for (j = 0; j < part_dim; j++) {
      if (sum[j] == 0.0) {
        sprintf(msg, "Sum of weights (component %1d) is zero.", j);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
        error = ZOLTAN_FATAL;
        goto End;
      }
    }

    /* Normalize part sizes */
    for (i = 0; i < num_global_parts; i++)
      for (j = 0; j < part_dim; j++)
        part_sizes[i*part_dim+j] /= sum[j];

  }
 
End:
  if (temp_part_sizes) ZOLTAN_FREE(&temp_part_sizes);
  if (sum)             ZOLTAN_FREE(&sum);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return error;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
