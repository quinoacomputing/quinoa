/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"defs.h"
#include	"structs.h"


/* Set up data structures for refine_part. */
void 
make_maps_ref (
    struct vtx_data **graph,	/* graph data structure */
    struct bilist *set_list,	/* lists of vertices in each set */
    struct bilist *vtx_elems,	/* start of storage for vertices */
    int *assignment,		/* set assignments for graph */
    int *sub_assign,		/* assignment file for subgraph */
    int set1,
    int set2,		/* set value denoting subgraph */
    int *glob2loc,		/* graph -> subgraph numbering map */
    int *loc2glob,		/* subgraph -> graph numbering map */
    int *psub_nvtxs,		/* number of vtxs in subgraph */
    int *pvwgt_max,		/* returned largest vwgt */
    int *pvwgt_sum1,
    int *pvwgt_sum2	/* returned set sizes */
)
{
    struct bilist *ptr;		/* loops through set lists */
    int       vwgt_max;		/* largest vertex weight in subgraph */
    int       vwgt_sum1, vwgt_sum2;	/* sum of vertex weights in sets */
    int       vtx;		/* vertex in subgraph */
    int       size;		/* array spacing */
    int       j;		/* loop counter */

    size = (int) (&(vtx_elems[1]) - &(vtx_elems[0]));
    j = 1;
    vwgt_max = vwgt_sum1 = vwgt_sum2 = 0;
    for (ptr = set_list[set1].next; ptr != NULL; ptr = ptr->next) {
	vtx = ((int) (ptr - vtx_elems)) / size;
	sub_assign[j] = 0;
	glob2loc[vtx] = j;
	loc2glob[j] = vtx;
	if (graph[vtx]->vwgt > vwgt_max) {
	    vwgt_max = graph[vtx]->vwgt;
	}
	vwgt_sum1 += graph[vtx]->vwgt;
	j++;
    }

    for (ptr = set_list[set2].next; ptr != NULL; ptr = ptr->next) {
	vtx = ((int) (ptr - vtx_elems)) / size;
	sub_assign[j] = 1;
	glob2loc[vtx] = j;
	loc2glob[j] = vtx;
	if (graph[vtx]->vwgt > vwgt_max) {
	    vwgt_max = graph[vtx]->vwgt;
	}
	vwgt_sum2 += graph[vtx]->vwgt;
	assignment[vtx] = (int) set1;
	j++;
    }
    *pvwgt_sum1 = vwgt_sum1;
    *pvwgt_sum2 = vwgt_sum2;
    *pvwgt_max = vwgt_max;
    *psub_nvtxs = j - 1;
}
