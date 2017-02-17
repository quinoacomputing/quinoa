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
//
// ************************************************************************
// 
// C++ wrappers for Zoltan communication library.  
//
// Two styles of initialization: 
//
//   C++ style:  Zoltan_Comm comm(nvals, assign, comm, tag, pnvals_recv);
//
//   C style: Zoltan_Comm comm;
//            comm.Create(nvals, assign, comm, tag, pnvals_recv);
// 
// ************************************************************************

#ifndef ZOLTAN_COMM_CPP_H_
#define ZOLTAN_COMM_CPP_H_

#include "zoltan_comm.h"

class Zoltan_Comm {

public:
  
  Zoltan_Comm(const int &nvals, int *assign, const MPI_Comm &comm, 
              const int &tag, int *pnvals_recv)
    {
    // Assumption: MPI has been initialized prior to this call.
    Zoltan_Comm_Create(&this->Plan, nvals, assign, comm, tag, pnvals_recv);
    }

  Zoltan_Comm()
    {
    this->Plan = NULL;
    // Caller will have to call Create to finish initialization of object
    }

  int Create(const int &nvals, int *assign, const MPI_Comm &comm, 
             const int &tag, int *pnvals_recv)
    {
    if (this->Plan)
      {
      Zoltan_Comm_Destroy(&this->Plan);
      this->Plan = NULL;
      }

    int rc = Zoltan_Comm_Create(&this->Plan, nvals, assign, comm, tag, pnvals_recv);

    return rc;
    }

  Zoltan_Comm (const Zoltan_Comm &plan) // Copy constructor
   {
   this->Plan = Zoltan_Comm_Copy(plan.Plan);
   }

  ~Zoltan_Comm()
    {
    Zoltan_Comm_Destroy(&this->Plan);
    }

  Zoltan_Comm & operator= (const Zoltan_Comm &plan) // Copy operator
  {
    Zoltan_Comm_Copy_To(&this->Plan, plan.Plan);

    return *this;
  }
      
  int Resize(int *sizes, const int &tag, int *sum_recv_sizes)
    {
    return Zoltan_Comm_Resize( this->Plan, sizes, tag, sum_recv_sizes);
    }
  
  int Do(const int &tag, char *send_data, const int &nbytes, char *recv_data)
    {
    return Zoltan_Comm_Do(this->Plan, tag, send_data, nbytes, recv_data);
    }

  int Do_Post( const int &tag, char *send_data, const int &nbytes, char *recv_data)
    {
    return Zoltan_Comm_Do_Post(this->Plan, tag, send_data, nbytes, recv_data);
    }

  int Do_Wait(const int &tag, char *send_data, const int &nbytes, char *recv_data)
    {
    return Zoltan_Comm_Do_Wait(this->Plan, tag, send_data, nbytes, recv_data);
    }
  
  int Do_Reverse(const int &tag, char *send_data, const int &nbytes, int *sizes, char *recv_data)
    {
    return Zoltan_Comm_Do_Reverse(this->Plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }

  int Do_Reverse_Post(const int &tag, char *send_data, const int &nbytes, int *sizes, 
    char *recv_data)
    {
    return Zoltan_Comm_Do_Reverse_Post(this->Plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }

  int Do_Reverse_Wait(const int &tag, char *send_data, const int &nbytes, int *sizes, 
    char *recv_data)
    {
    return Zoltan_Comm_Do_Reverse_Wait(this->Plan, tag, send_data, nbytes, sizes, 
        recv_data);
    }
  
  int Info( int *nsends, int *send_procs,
    int *send_lengths, int *send_nvals, int *send_max_size, int *send_list,
    int *nrecvs, int *recv_procs, int *recv_lengths, int *recv_nvals,
    int *recv_total_size, int *recv_list, int *self_msg) const
    {
      return Zoltan_Comm_Info( this->Plan, nsends, send_procs, send_lengths, 
        send_nvals, send_max_size, send_list, nrecvs, recv_procs, recv_lengths, 
        recv_nvals, recv_total_size, recv_list, self_msg);
    }
  
  int Invert_Plan()
    {
    return Zoltan_Comm_Invert_Plan(&this->Plan);
    }

  // Static methods

  static int Invert_Map( int *lengths_to, int *procs_to, 
    const int &nsends, const int &self_msg,
    int * &plengths_from, int * &pprocs_from, int &pnrecvs, 
    const int &my_proc, const int &nprocs, const int &out_of_mem, 
    const int &tag, const MPI_Comm &comm) 
    {
    return Zoltan_Comm_Invert_Map( lengths_to, procs_to, nsends, self_msg,
        &plengths_from, &pprocs_from, &pnrecvs, my_proc, nprocs, out_of_mem, 
        tag, comm);
    }
      
  static int Exchange_Sizes( int *sizes_to, int *procs_to, 
    const int &nsends, const int &self_msg,
    int *sizes_from, int *procs_from, const int &nrecvs, int *total_recv_size,
    const int &my_proc, const int &tag, const MPI_Comm &comm) 
    {
    return Zoltan_Comm_Exchange_Sizes(sizes_to, procs_to, nsends, self_msg,
        sizes_from, procs_from, nrecvs, total_recv_size, my_proc, tag, 
        comm);
    }

private:

  ZOLTAN_COMM_OBJ *Plan;
};

#endif
