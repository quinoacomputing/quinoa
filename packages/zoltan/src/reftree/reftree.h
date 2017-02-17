/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
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

#ifndef __REFTREE_H
#define __REFTREE_H

#include "reftree_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Some constants */

/* Maximum number of vertices per element */
/* Used for dimensioning space for a query function to return vertices */
#define MAXVERT 8

/* Default parameter values */
#define DEFAULT_HASH_TABLE_SIZE 16384
#define DEFAULT_INITPATH "REFTREE_DEFAULT"

/* Data structures for refinement tree */

/* The main refinement tree structure */

struct Zoltan_Reftree_Struct {
   ZOLTAN_ID_PTR global_id;  /* global ID of the corresponding element */
   ZOLTAN_ID_PTR local_id;   /* local ID of the corresponding element */
   struct Zoltan_Reftree_Struct *children; /* array of the children in the tree */
   int num_child;        /* number of children */
   float *weight;        /* weight of the node; dimension Obj_Weight_Dim */
   float *summed_weight; /* sum of the weights in the subtree rooted at
                            this node */
   float *my_sum_weight; /* sum of weights of nodes assigned to this proc */
   int num_vertex;       /* the number of vertices in the corresponding
                            element */
   ZOLTAN_ID_PTR vertices; /* the vertices of the corresponding element;
                            local to this processor */
   ZOLTAN_ID_PTR in_vertex; /* starting vertex for determining the path
                            through the children */
   ZOLTAN_ID_PTR out_vertex; /* ending vertex for determining the path
                            through the children */
   int assigned_to_me;   /* for a leaf, 1 if this element is assigned to
                            this processor, 0 if not.  for nonleaves, 1 if
                            the entire subtree is assigned to this proc,
                            0 if none of the subtree, -1 if part */
   int known_to_me;      /* for coarse grid objects, 1 if it is known to this
                            processor (i.e., returned by get_coarse_obj) and
                            0 if not */
   int partition;        /* partition to which this node is assigned;
                            meaningful only during the partition algorithm */
};

typedef struct Zoltan_Reftree_Struct ZOLTAN_REFTREE;

/* Hash table structures */

struct Zoltan_Reftree_hash_node {
  ZOLTAN_ID_PTR gid;            /* Global id */
  ZOLTAN_REFTREE *reftree_node; /* pointer to a node of the refinement tree */
  struct Zoltan_Reftree_hash_node *next;
};

struct Zoltan_Reftree_inthash_node {
  ZOLTAN_ID_PTR gid;            /* Global id */
  int lid;                      /* integer corresponding to the gid */
  struct Zoltan_Reftree_inthash_node *next;
};

/* data structure pointed to by zz->Data_Structure */

struct Zoltan_Reftree_data_struct {
  ZOLTAN_REFTREE *reftree_root;
  struct Zoltan_Reftree_hash_node **hash_table;
  int hash_table_size;
};

/* Prototypes */

extern int Zoltan_Reftree_Init(ZZ *zz);
extern int Zoltan_Reftree_Build(ZZ *zz);
extern void Zoltan_Reftree_Print(ZZ *zz,ZOLTAN_REFTREE *subroot, int level);

extern int Zoltan_Reftree_Coarse_Grid_Path(int nobj, int *num_vert,
                               ZOLTAN_ID_PTR vertices, ZOLTAN_ID_PTR in_vertex,
                               ZOLTAN_ID_PTR out_vertex, double *coords,
                               int *order, ZOLTAN_ID_PTR gids,
                               ZOLTAN_ID_PTR lids, char *initpath_method,
                               ZZ *zz);

extern ZOLTAN_REFTREE* Zoltan_Reftree_hash_lookup(ZZ *zz, 
                                      struct Zoltan_Reftree_hash_node **hashtab,
                                      ZOLTAN_ID_PTR key, int n);
extern int Zoltan_Reftree_inthash_lookup(ZZ *zz, 
                                   struct Zoltan_Reftree_inthash_node **hashtab,
                                   ZOLTAN_ID_PTR key, int n);
extern void Zoltan_Reftree_Hash_Insert(ZZ *zz, ZOLTAN_REFTREE *reftree_node,
                          struct Zoltan_Reftree_hash_node **hashtab, int size);
extern void Zoltan_Reftree_IntHash_Insert(ZZ *zz, ZOLTAN_ID_PTR gid, int lid,
                        struct Zoltan_Reftree_inthash_node **hashtab, int size);
extern void Zoltan_Reftree_Hash_Remove(ZZ *zz, ZOLTAN_REFTREE *reftree_node,
                          struct Zoltan_Reftree_hash_node **hashtab, int size);
extern void Zoltan_Reftree_Clear_Hash_Table(
                          struct Zoltan_Reftree_hash_node **hashtab, int size);
extern void Zoltan_Reftree_Clear_IntHash_Table(
                       struct Zoltan_Reftree_inthash_node **hashtab, int size);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* __REFTREE_CONST_H */
