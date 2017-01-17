/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function(s) contained in this file:
 *
 *      ex_put_partial_node_set_df()
 *
 *****************************************************************************
 *
 *  Variable Index:
 *
 *      exoid               - The NetCDF ID of an already open NemesisI file.
 *      node_set_id        - ID of node set to write.
 *      start_num          - The starting index of the df's to be written.
 *      num_df_to_get      - The number of distribution factors to write out.
 *      node_set_dist_fact - List of node distribution factors in node set.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include "exodusII.h"     // for ex_err, exerrval, etc
#include "exodusII_int.h" // for EX_FATAL, EX_WARN, etc
#include "netcdf.h"       // for NC_NOERR, nc_inq_dimid, etc
#include <inttypes.h>     // for PRId64
#include <stddef.h>       // for size_t
#include <stdio.h>
#include <sys/types.h> // for int64_t

/*
 * writes the node set distribution factors for a single node set
 */

int ex_put_partial_node_set_df(int exoid, ex_entity_id node_set_id, int64_t start_num,
                               int64_t num_df_to_get, void *node_set_dist_fact)
{
  int    status;
  int    dimid, dist_id, node_set_id_ndx;
  size_t num_nodes_in_set, start[1], count[1];
  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* first check if any node sets are specified */

  if ((status = nc_inq_dimid(exoid, DIM_NUM_NS, &dimid)) < 0) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: no node sets specified in file id %d", exoid);
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Lookup index of node set id in VAR_NS_IDS array */
  if ((node_set_id_ndx = ex_id_lkup(exoid, EX_NODE_SET, node_set_id)) < 0) {
    if (exerrval == EX_NULLENTITY) {
      snprintf(errmsg, MAX_ERR_LENGTH,
               "Warning: no data allowed for NULL node set %" PRId64 " in file id %d", node_set_id,
               exoid);
      ex_err("ex_put_partial_node_set_df", errmsg, EX_NULLENTITY);
      return (EX_WARN);
    }
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate node set id %" PRId64 " in VAR_NS_IDS array in file id %d",
             node_set_id, exoid);
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* inquire id's of previously defined dimensions  */
  if ((status = nc_inq_dimid(exoid, DIM_NUM_NOD_NS(node_set_id_ndx), &dimid)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate number of nodes in node set %" PRId64 " in file id %d",
             node_set_id, exoid);
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  if ((status = nc_inq_dimlen(exoid, dimid, &num_nodes_in_set)) != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to get number of nodes in set %" PRId64 " in file id %d", node_set_id,
             exoid);
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* Check input parameters for a valid range of numbers */
  if (start_num < 0 || start_num > num_nodes_in_set) {
    exerrval = EX_BADPARAM;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Invalid input");
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  if (num_df_to_get < 0) {
    exerrval = EX_BADPARAM;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: Invalid number of nodes in nodes set!");
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* start_num now starts at 1, not 0 */
  if ((start_num + num_df_to_get - 1) > num_nodes_in_set) {
    exerrval = EX_BADPARAM;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: request larger than number of nodes in set!");
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* find id of distribution factors variable */
  if ((status = nc_inq_varid(exoid, VAR_FACT_NS(node_set_id_ndx), &dist_id)) != NC_NOERR) {
    if (status == NC_ENOTVAR) {
      exerrval = EX_BADPARAM;
      snprintf(errmsg, MAX_ERR_LENGTH,
               "Warning: no dist factors defined for node set %" PRId64 " in file id %d",
               node_set_id, exoid);
      ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
      return (EX_WARN);
    }
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to locate node set %" PRId64 " dist factors in file id %d", node_set_id,
             exoid);
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }

  /* write out the distribution factors array */
  start[0] = --start_num;
  count[0] = num_df_to_get;
  if (count[0] == 0) {
    start[0] = 0;
  }

  if (ex_comp_ws(exoid) == 4) {
    status = nc_put_vara_float(exoid, dist_id, start, count, node_set_dist_fact);
  }
  else {
    status = nc_put_vara_double(exoid, dist_id, start, count, node_set_dist_fact);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH,
             "ERROR: failed to store node set %" PRId64 " dist factors in file id %d", node_set_id,
             exoid);
    ex_err("ex_put_partial_node_set_df", errmsg, exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
