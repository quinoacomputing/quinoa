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
#include "zoltan_types.h"
#include "zoltan_id.h"
#include "zoltan_util.h"
#include "zoltan_mem.h"

/*****************************************************************************/
/*
 *  This file contains routines for manipulating 
 *  the global and local IDs used by Zoltan.
 *
 *  Some manipulations are performed via macros.  In particular, macros
 *  specifying whether global or local IDs are to be manipulated are
 *  provided.  Use of these macros is recommended over use of these
 *  basic functions.
 *  See zoltan_id.h for definitions of these macros.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routines for allocating and initializing IDs.
 */

ZOLTAN_ID_PTR ZOLTAN_Malloc_ID(int n, char *file, int line)
{
/*
 * Allocates an array of size n of ZOLTAN_ID_TYPEs and initializes them.
 */

ZOLTAN_ID_PTR tmp;
char *yo = "ZOLTAN_Malloc_ID";

  /* 
   * Don't use ZOLTAN_MALLOC macro here; prefer to pass file and line 
   * where ZOLTAN_Malloc_ID was called.
   */
  tmp = (ZOLTAN_ID_PTR) Zoltan_Malloc(n * sizeof(ZOLTAN_ID_TYPE), file, line);

  if (tmp != NULL) {
    ZOLTAN_INIT_ID(n,tmp);
  }
  else if (n > 0) {
    char msg[256];
    sprintf(msg, "NULL pointer returned; malloc called from %s, line %d.",
            file, line);
    ZOLTAN_PRINT_ERROR(-1, yo, msg);
  }

  return tmp;
}

/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routines for printing IDs.
 */

void ZOLTAN_PRINT_ID(int n, ZOLTAN_ID_PTR a)
{
/* Prints a single ID. */
int i;

  printf("(");
  for (i = 0; i < n; i++)
    printf( ZOLTAN_ID_SPEC " ",a[i]);
  printf(") ");
}
  
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Routines to compare Global IDs.  
 *  Functions are provided to test whether two IDs are equal (EQ),
 *  less than (LT), and greater than (GT).
 *  The negation operator can be used to test whether two IDs are 
 *  not equal (!ZOLTAN_EQ_ID(n,a,b)), less than or equal (!ZOLTAN_GT_GID(n,a,b))
 *  or greater than or equal (!ZOLTAN_LT_GID(n,a,b)).
 */
/*****************************************************************************/

int ZOLTAN_EQ_ID(int n, ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b)
{
/* 
 * Returns 1 if a == b; 0 otherwise.
 * a == b if for all i, a[i] == b[i].
 */
int i;
  for (i = 0; i < n; i++)
    if (a[i] != b[i])
      return(0);
  return(1);
}

/*****************************************************************************/

#ifdef ZOLTAN_NEEDED
/* Commented out since never used. */

int ZOLTAN_LT_ID(int n, ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b)
{
/* 
 * Returns 1 if a < b; 0 otherwise.
 * a < b if for some i, a[i] < b[i] and a[j] == b[j] for all j < i.
 */
int i;
  for (i = 0; i < n; i++)
    if (a[i] == b[i])
      continue;
    else if (a[i] > b[i])
      return(0);
    else /* a[i] < b[i] */
      return(1);

  return(0); /* because a == b */
}

/*****************************************************************************/

int ZOLTAN_GT_ID(int n, ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b)
{
/* 
 * Returns 1 if a < b; 0 otherwise.
 * a > b if for some i, a[i] > b[i] and a[j] == b[j] for all j < i.
 */
int i;
  for (i = 0; i < n; i++)
    if (a[i] == b[i])
      continue;
    else if (a[i] < b[i])
      return(0);
    else /* a[i] > b[i] */
      return(1);

  return(0); /* because a == b */
}
#endif /* ZOLTAN_NEEDED */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
