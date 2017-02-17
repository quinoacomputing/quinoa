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
extern "C" {
#endif

/* Communication routines for parallel hypergraph partitioning. */
/* UNDER CONSTRUCTION */

#include "phg.h"
#include "zz_const.h"

#if 0
/*
 * note: we need to change this so that we keep the plan around for
 * multiple uses.  what is the best way?
 */

/*
 * gather_by_list
 * input:   a list of processors, a message to be sent to these processors,
 *          and a receive buffer to get data back.
 *
 *  output: concatenated messages from other processors in rbuff, with
 *          appropriate size info in rbuff_size.
 */
int gather_by_list(int procs_length, int   *procs,
		   int sbuff_size,   char  *sbuff,
		   int *rbuff_size,  char **rbuff, MPI_Comm *comm)
{
    int i, receive_size, mtag, err, myProc;
    char            *send;
    ZOLTAN_COMM_OBJ *plan;
    static char     *yo = "gather_by_list";

    MPI_Comm_rank(*comm, &myProc);

    if (!(send  = (char*) ZOLTAN_MALLOC(procs_length * sizeof(int)))) {
	ZOLTAN_PRINT_ERROR(myProc, yo, "Insufficient memory");
	return ZOLTAN_MEMERR;
    }

    mtag = 0;
    for (i = 0; i < procs_length; ++i)
	((int*)send)[i] = sbuff_size;

    /* create and resize communication plan */
    err = Zoltan_Comm_Create(&plan, procs_length, procs, *comm,
	    mtag++, &receive_size);
    if (err < 0)
	ZOLTAN_PRINT_ERROR(myProc, yo, "Zoltan_Comm_Create failed.");
    err = Zoltan_Comm_Resize(plan, (int*)send, mtag++, rbuff_size);
    if (err < 0)
	ZOLTAN_PRINT_ERROR(myProc, yo, "Zoltan_Comm_Resize failed.");

    ZOLTAN_FREE(&send);

    /* allocate send and receive buffer */
    if ((!(*rbuff = (char*) ZOLTAN_MALLOC(*rbuff_size)) && (*rbuff_size != 0))
      ||!(  send  = (char*) ZOLTAN_MALLOC(sbuff_size * procs_length))) {
	printf("failed allocating %d + %d * %d bytes. . .", *rbuff_size, sbuff_size, procs_length);
	fflush(NULL);
	ZOLTAN_PRINT_ERROR(myProc, yo, "Insufficient Memory");
	return ZOLTAN_MEMERR;
    }

    /* copy message for each processor */
    for (i = 0; i < procs_length; ++i)
	memcpy(send + i * sbuff_size, sbuff, sbuff_size);

    /* finally, we can do the communication */
    err = Zoltan_Comm_Do(plan, mtag++, send, 1, *rbuff);
    if (err < 0)
	ZOLTAN_PRINT_ERROR(myProc, yo, "Zoltan_Comm_Do failed.");

    /* clean up */
    Zoltan_Comm_Destroy(&plan);
    ZOLTAN_FREE(&send);
    return err;
}

/*
 * gather_row
 * input:   our processor location, a message, and space to hold a return
 *          message.
 * output:  received messages from each processor in our row in rbuff.
 */
int gather_row(int nProc_x, int nProc_y, int myProc_x, int myProc_y,
	       int sbuff_size, char *sbuff, int *rbuff_size, char **rbuff,
	       MPI_Comm *comm)
{
    int  i, my_proc, err;
    int  *procs;
    static char *yo = "gather_row";

    /* create list of processors to send to (everyone in our row but us)     */
    /* this assumes a linear left->right numbering of processors in our grid */
    myProc_y *= nProc_x;
    my_proc = myProc_x + myProc_y;

    if (!(procs = (int*)  ZOLTAN_MALLOC((nProc_x) * sizeof(int)))) {
	ZOLTAN_PRINT_ERROR(my_proc, yo, "Insufficient memory");
	return ZOLTAN_MEMERR;
    }

    /*
      for (i = 0; i < myProc_x; ++i)
	  procs[i] = myProc_y + i;
      for (i = myProc_x; i < nProc_x - 1; ++i)
	  procs[i] = myProc_y + i + 1;
    */

    /* do self communication */
    for (i = 0; i < nProc_x; ++i)
      procs[i] = myProc_y + i;

    err = gather_by_list(nProc_x, procs,
			 sbuff_size, sbuff, rbuff_size, rbuff, comm);

    ZOLTAN_FREE(&procs);
    return err;
}

/*
 * gather_col
 * input:   our processor location, a message, and space to hold a return
 *          message.
 * output:  received messages from each processor in our row in rbuff.
 */
int gather_col(int nProc_x, int nProc_y, int myProc_x, int myProc_y,
	       int sbuff_size, char *sbuff, int *rbuff_size, char **rbuff,
	       MPI_Comm *comm)
{
    int  i, my_proc, err;
    int  *procs;
    static char *yo = "gather_col";

    my_proc = myProc_y * nProc_x + myProc_x;

    if (!(procs = (int*) ZOLTAN_MALLOC(nProc_y * sizeof(int)))) {
	ZOLTAN_PRINT_ERROR(my_proc, yo, "Insufficient memory");
	return ZOLTAN_MEMERR;
    }

    /* create list of processors to send to (everyone in our row but us)     */
    /* this assumes a linear left->right numbering of processors in our grid */
    /*
      for (i = 0; i < myProc_y; ++i)
	  procs[i] = nProc_x * i + myProc_x;
      for (i = myProc_y; i < nProc_y - 1; ++i)
	  procs[i] = nProc_x * (i + 1) + myProc_x;
    */

    /* do self communication */
    for (i = 0; i < nProc_y; ++i)
	procs[i] = nProc_x * i + myProc_x;

    err = gather_by_list(nProc_y, procs,
			 sbuff_size, sbuff, rbuff_size, rbuff, comm);

    ZOLTAN_FREE(&procs);
    return err;
}

/*
 * gather_col_root
 * input:   our processor location, a message, and space to hold a return
 *          message.
 * output:  on non-root processors (ones not at the base of a column), nothing.
 *          on root processors, received messages from each processor in the
 *          column in rbuff.
 */
int gather_col_root(int nProc_x, int nProc_y, int myProc_x, int myProc_y,
   int sbuff_size, char *sbuff, int *rbuff_size, char **rbuff, MPI_Comm *comm)
{
    int  *procs;
    static char *yo = "gather_col_root";

    if (!(procs = (int*) ZOLTAN_MALLOC(sizeof(int)))) {
	ZOLTAN_PRINT_ERROR(myProc_y * nProc_x + myProc_x, yo,
		"Insufficient Memory");
	return ZOLTAN_MEMERR;
    }
    *procs = myProc_x;
    return gather_by_list(1, procs, sbuff_size, sbuff, rbuff_size, rbuff, comm);
}

/*
 * gather_row_root
 * input:   our processor location, a message, and space to hold a return
 *          message.
 * output:  on non-root processors (ones not at the start of a row), nothing.
 *          on root processors, received messages from each processor in the
 *          row in rbuff.
 */
int gather_row_root(int nProc_x, int nProc_y, int myProc_x, int myProc_y,
   int sbuff_size, char *sbuff, int *rbuff_size, char **rbuff, MPI_Comm *comm)
{
    int *procs;
    static char *yo = "gather_row_root";

    if(!(procs = (int*) ZOLTAN_MALLOC(sizeof(int)))) {
	ZOLTAN_PRINT_ERROR(myProc_y * nProc_x + myProc_x, yo,
		"Insufficient Memory");
	return ZOLTAN_MEMERR;
    }
    *procs = myProc_y;
    return gather_by_list(1, procs, sbuff_size, sbuff, rbuff_size, rbuff, comm);
}

/*
 * gather_slice
 * input:   our processor location, a message, space to hold a return message,
 *          and a flag indicating whether to communicate by row or by column
 *          (0 indicates column, non-zero indicates row)
 * output:  received messages concatenated in rbuff.
 */
int Zoltan_PHG_gather_slice(int nProc_x, int nProc_y, int myProc_x,
	int myProc_y,int sbuff_size, char *sbuff, int *rbuff_size,
	char ** rbuff, MPI_Comm *comm, int do_row)
{
    if (do_row)
	return gather_row(nProc_x, nProc_y, myProc_x, myProc_y, sbuff_size,
			  sbuff, rbuff_size, rbuff, comm);
    else
	return gather_col(nProc_x, nProc_y, myProc_x, myProc_y, sbuff_size,
			  sbuff, rbuff_size, rbuff, comm);
}

/*
 * gather_slice_root
 * input:   our processor location, a message, space to hold a return message,
 *          and a flag indicating whether to communicate by row or by column.
 * output:  in non-root processors, nothing.
 *          in root processors, messages from all processors in the appropriate
 *          row or column, concatenated in rbuff.
 */
int Zoltan_PHG_gather_slice_root(int nProc_x, int nProc_y, int myProc_x,
	int myProc_y, int sbuff_size, char *sbuff, int *rbuff_size,
	char **rbuff, MPI_Comm *comm, int do_row)
{
    if (do_row)
	return gather_row_root(nProc_x, nProc_y, myProc_x, myProc_y,
		sbuff_size, sbuff, rbuff_size, rbuff, comm);
    else
	return gather_col_root(nProc_x, nProc_y, myProc_x, myProc_y,
		sbuff_size, sbuff, rbuff_size, rbuff, comm);
}

#endif


void
Zoltan_PHGComm_Destroy(PHGComm* comm)
{
  if (comm == NULL)
    return;
  if(comm->Communicator != 0 && comm->Communicator != MPI_COMM_NULL)
    MPI_Comm_free (&comm->Communicator);
  if(comm->row_comm != 0 && comm->row_comm != MPI_COMM_NULL)
    MPI_Comm_free (&comm->row_comm);
  if(comm->col_comm != 0 && comm->col_comm != MPI_COMM_NULL)
    MPI_Comm_free (&comm->col_comm);
}

void
Zoltan_PHGComm_Init(PHGComm* comm)
{
  comm->Communicator = MPI_COMM_NULL;
  comm->row_comm = MPI_COMM_NULL;
  comm->col_comm = MPI_COMM_NULL;
}


#ifdef __cplusplus
} /* close extern "C" */
#endif
