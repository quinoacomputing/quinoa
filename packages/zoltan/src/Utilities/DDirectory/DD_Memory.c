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


#include <stdio.h>
#include <stdlib.h>

#include "zoltan_dd_const.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*
 * Idea:  Store hash table as array of bucket heads, pointing into
 * an array-based linked list of nodes (instead of to individually 
 * allocated nodes).
 * The array-based memory would have a free-memory list running
 * through it, allowing items to be added to or deleted from the 
 * directory.
 */

/******************************************************************************/
int DD_Memory_Alloc_Nodelist(
  Zoltan_DD_Directory *dd,  /* directory state information    */
  DD_NodeIdx count,         /* Number of GIDs in update list  */
  float overalloc           /* Percentage to extra nodes to
                               allocate (for future dynamic 
                               additions).                    */
)
{
/* Allocate node memory and initialize it to be all free. */
/* Return error code if memory alloc fails.               */
  int ierr = ZOLTAN_OK;
  DD_Node *ptr;
  char *dataptr;
  DD_NodeIdx nodeidx;
  DD_NodeIdx len = count * (1. + overalloc);

  dd->nodelistlen = len;

  if (len > 0) {

    dd->nodelist = (DD_Node *) ZOLTAN_MALLOC(sizeof(DD_Node) * len);
    dd->nodedata = (char *) ZOLTAN_MALLOC(dd->nodedata_size * len);
    /* TODO ADD ERROR CHECK */

    /* Initialize the freenode list; all nodes are initially free. */
    /* Also initialize gid pointers to nodedata.                   */
    dd->nextfreenode = 0;
    for (ptr = dd->nodelist, dataptr = dd->nodedata, nodeidx = 0; 
         nodeidx < len;
         nodeidx++, ptr++, dataptr += dd->nodedata_size) {
      ptr->next = nodeidx + 1;
      ptr->gid = (ZOLTAN_ID_PTR) dataptr;
      ptr->free = 1;
    }
    dd->nodelist[len-1].next = -1;  /* NULL value at end of list */
  }
  dd->nodecnt = 0;

  return ierr;
}

/******************************************************************************/

DD_NodeIdx DD_Memory_Alloc_Node(Zoltan_DD_Directory *dd)
{
/* "allocate" a node by removing it from the free list and returning
 * its NodeIdx.
 */
DD_NodeIdx returnnode;

  if (dd->nextfreenode == -1) {
    /* No more room for new nodes in the node list.
     * Realloc memory and set up a new freelist.
     */
    DD_NodeIdx newlen = dd->nodelistlen * 2;
    DD_NodeIdx nodeidx;
    DD_Node *ptr;
    char *dataptr;

    dd->nodelist = (DD_Node *) ZOLTAN_REALLOC(dd->nodelist,
                                              sizeof(DD_Node) * newlen);
    dd->nodedata = (char *) ZOLTAN_REALLOC(dd->nodedata,
                                           dd->nodedata_size * newlen);
    /* TODO ADD ERROR CHECK */

    /* Reinitialize the gid pointers in the realloc'ed nodelist. */
    for (ptr = dd->nodelist, dataptr = dd->nodedata, nodeidx = 0;
         nodeidx < newlen;
         ptr++, dataptr += dd->nodedata_size, nodeidx++) 
      ptr->gid = (ZOLTAN_ID_PTR) dataptr;

    /* Initialize free list in the newly extended memory */
    dd->nextfreenode = dd->nodelistlen;
    for (nodeidx = dd->nodelistlen; nodeidx < newlen-1; nodeidx++) {
      dd->nodelist[nodeidx].next = nodeidx+1;
      dd->nodelist[nodeidx].free = 1;
    }
    dd->nodelist[newlen-1].next = -1;
    dd->nodelist[newlen-1].free = 1;
    dd->nodelistlen = newlen;
  }

  returnnode = dd->nextfreenode;
  dd->nextfreenode = dd->nodelist[returnnode].next;
  dd->nodelist[returnnode].next = -1;
  dd->nodelist[returnnode].free = 0;
  dd->nodecnt++;

  return returnnode;
}

/******************************************************************************/

void DD_Memory_Free_Node(
  Zoltan_DD_Directory *dd,
  DD_NodeIdx freenode
)
{
/* "free" a node by returning it to the free list.  */
  
  /* TODO Error check:  freenode should be < nodelistlen */

  dd->nodelist[freenode].next = dd->nextfreenode;
  dd->nodelist[freenode].free = 1;
  dd->nextfreenode = freenode;
  dd->nodecnt--;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
