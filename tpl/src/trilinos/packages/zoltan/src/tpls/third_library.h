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


#ifndef __THIRD_LIBRARY_H
#define __THIRD_LIBRARY_H

#include <limits.h>
#include "zoltan_comm.h"
#include "third_library_const.h"
#include "graph.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#define INT_EPSILON (1e-5) /* How close we need to be to say it's an integer */

/* Macro for error handling */
#define ZOLTAN_PARMETIS_ERROR(error,str) {ierr = error ; \
 ZOLTAN_PRINT_ERROR(zz->Proc, yo, str) ; goto End ;}

#define ZOLTAN_THIRD_ERROR(error,str) { \
 ZOLTAN_PRINT_ERROR(zz->Proc, yo , str) ; return(error) ;}


/* Structure that defines a graph for third party libraries like ParMetis. */
typedef struct ZOLTAN_Third_Graph_ {
  int graph_type;                       /* Type of the graph */
  int check_graph;                      /* We have to check graph consistency */
  int final_output;                     /* Do some output computations */
  int showMoveVol;                      /* How works the final output */
  int scatter;                          /* Graph has been scattered */
  int scatter_min;                      /* Minimum level of scatter */
  int get_data;                         /* Construct edge datas */
  int obj_wgt_dim;                      /* Number of weights by vertex */
  int edge_wgt_dim;                     /* Number of weights by edge */
  int num_obj;                          /* Local number of vertices */
  int num_obj_orig;                     /* Local number of vertices in original graph */
  int num_edges;                        /* Local number of edges */
  indextype *vtxdist;                   /* How vertices are distributed */
  indextype * xadj;                     /* Indexes on adjency array */
  indextype *adjncy;                    /* adjency array (CSR) */
  weighttype * vwgt;                    /* Array of vertex weights */
  weighttype * ewgts;                   /* Array of edge weights */
  float * float_ewgts;
  int * adjproc;                        /* Array of procs ? */
  ZOLTAN_COMM_OBJ *comm_plan;           /* Communication plan used by scattering process */
  ZG graph;
} ZOLTAN_Third_Graph;

/* Structure that defines a geometry for third party libraries like ParMetis. */
typedef struct ZOLTAN_Third_Geom_ {
  int ndims;                            /* Number of dimensions */
  realtype *xyz;                        /* Coordinates */
} ZOLTAN_Third_Geom;

/* Structure that defines a partition for third party libraries like ParMetis. */
typedef struct ZOLTAN_Third_Part_ {
  indextype *part;                      /* Current partition */
  indextype *input_part;
  indextype *part_orig;
  realtype  *part_sizes;
  realtype  *input_part_sizes;
} ZOLTAN_Third_Part;

typedef struct ZOLTAN_Third_Vsize_ {
  int vsize_malloc;
  indextype *vsize;
  indextype *vsizeBACKUP;
} ZOLTAN_Third_Vsize;

/* Structure that defines an ordering output for third party libraries like ParMetis. */
typedef struct ZOLTAN_Output_Order_ {
  int num_part;
  indextype start_index;
  indextype *rank;      /* rank[i] is the rank of gids[i] */
  indextype *iperm;     /* inverse permutation of rank */
  ZOOS *order_opt;	/* ordering options */
  ZTPL_OS *order_info;	/* ordering info */
  indextype *sep_sizes;
} ZOLTAN_Output_Order;

/* Structure that describes a partitioning output in Zoltan lingo
   for third party libraries like ParMetis. */
typedef struct ZOLTAN_Output_Part_ {
  int  compute_only_part_changes;

  int  num_imp;
  ZOLTAN_ID_PTR* imp_gids;
  ZOLTAN_ID_PTR* imp_lids;
  int **imp_procs;
  int **imp_part;

  int  num_exp;
  ZOLTAN_ID_PTR *exp_gids;
  ZOLTAN_ID_PTR *exp_lids;
  int **exp_procs;
  int **exp_part;
} ZOLTAN_Output_Part;

/*****************************
 *  Functions which can be used to interface with third party libraries
 ****************************************************************************/

/* Construct graph and associated datas in order to call third library */
int Zoltan_Preprocess_Graph(
  ZZ *zz,                               /* Zoltan structure */
  ZOLTAN_ID_PTR *global_ids,
  ZOLTAN_ID_PTR *local_ids,
  ZOLTAN_Third_Graph *gr,              /* Graph for third part libs */
  ZOLTAN_Third_Geom  *geo,
  ZOLTAN_Third_Part  *prt,
  ZOLTAN_Third_Vsize *vsp
  );

/* Activate timers if requested */
int Zoltan_Preprocess_Timer(ZZ *zz, int *use_timer);

/* Display timing informations */
void Zoltan_Third_DisplayTime(ZZ* zz, double* times);

/* Initialize Zoltan internal structures */
int Zoltan_Third_Init(ZOLTAN_Third_Graph *gr, ZOLTAN_Third_Part  *prt, ZOLTAN_Third_Vsize *vsp, ZOLTAN_Output_Part *part,
		      ZOLTAN_ID_PTR *imp_gids, ZOLTAN_ID_PTR *imp_lids, int **imp_procs, int **imp_to_part,
		      ZOLTAN_ID_PTR *exp_gids, ZOLTAN_ID_PTR *exp_lids, int **exp_procs, int **exp_to_part);


/* export to user variables */
int Zoltan_Third_Export_User(ZOLTAN_Output_Part *part,
		     int *num_imp, ZOLTAN_ID_PTR *imp_gids, ZOLTAN_ID_PTR *imp_lids, int **imp_procs, int **imp_to_part,
		     int *num_exp, ZOLTAN_ID_PTR *exp_gids, ZOLTAN_ID_PTR *exp_lids, int **exp_procs, int **exp_to_part);



/* Free temporary datas */
void Zoltan_Third_Exit(ZOLTAN_Third_Graph *gr, ZOLTAN_Third_Geom *geo,
		       ZOLTAN_Third_Part *prt, ZOLTAN_Third_Vsize *vsp,
		       ZOLTAN_Output_Part *part, ZOLTAN_Output_Order *ord);

/* Do some postprocesing in order to exploit third party results in Zoltan */
int Zoltan_Postprocess_Graph(
  ZZ *zz,                               /* Zoltan structure */
  ZOLTAN_ID_PTR      global_ids,
  ZOLTAN_ID_PTR      local_ids,
  ZOLTAN_Third_Graph *gr,              /* Graph for third part libs */
  ZOLTAN_Third_Geom  *geo,
  ZOLTAN_Third_Part  *prt,
  ZOLTAN_Third_Vsize *vsp,
  ZOLTAN_Output_Order *ord,
  ZOLTAN_Output_Part  *part);

int
Zoltan_Postprocess_FinalOutput (ZZ* zz, ZOLTAN_Third_Graph *gr,
				ZOLTAN_Third_Part *prt, ZOLTAN_Third_Vsize *vsp,
				int use_timers, realtype itr);


int Zoltan_matrix_Print(Zoltan_matrix *m, char *s);
int Zoltan_Third_Graph_Print(ZZ *zz, ZOLTAN_Third_Graph *gr, char *s);
int Zoltan_ZG_Print(ZZ *zz, ZG *gr, char *s);
int Zoltan_Check_TPL_Data_Sizes(ZZ *, int);


#ifdef __cplusplus
}
#endif

#endif
