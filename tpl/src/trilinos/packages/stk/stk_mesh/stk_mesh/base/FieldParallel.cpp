// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_util/stk_config.h>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll, CommBuffer
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/util/PairIter.hpp>   // for PairIter
#include <stk_util/util/SortAndUnique.hpp>   // for PairIter

#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc

#include <utility>                      // for pair
#include <sstream>                      // for basic_ostream::operator<<, etc

namespace stk {
namespace mesh {


//----------------------------------------------------------------------
void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const BulkData & mesh = ghosts.mesh();
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();
  const unsigned ghost_id = ghosts.ordinal();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.internal_comm_list().begin() , iend = mesh.internal_comm_list().end(); i != iend ; ++i ) {
    Entity e = i->entity;
    const MeshIndex meshIdx = mesh.mesh_index(e);
    const unsigned bucketId = meshIdx.bucket->bucket_id();

    const bool owned = i->owner == parallel_rank ;

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;

      if(is_matching_rank(f, *meshIdx.bucket)) {
        e_size += field_bytes_per_entity( f , bucketId );
      }
    }

    if (e_size == 0) {
      continue;
    }

    const EntityCommInfoVector& infovec = i->entity_comm->comm_map;
    PairIterEntityComm ec(infovec.begin(), infovec.end());
    if ( owned ) {
      for ( ; ! ec.empty() ; ++ec ) {
        if (ec->ghost_id == ghost_id) {
          send_size[ ec->proc ] += e_size ;
        }
      }
    }
    else {
      for ( ; ! ec.empty() ; ++ec ) {
        if (ec->ghost_id == ghost_id) {
          recv_size[ i->owner ] += e_size ;
          break;//jump out since we know we're only recving 1 msg from the 1-and-only owner
        }
      }
    }
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const snd_size = send_size.data() ;
    const unsigned * const rcv_size = recv_size.data() ;
    sparse.allocate_buffers( mesh.parallel(), snd_size, rcv_size);
  }

  // Send packing:

  for (int phase = 0; phase < 2; ++phase) {

    for ( EntityCommListInfoVector::const_iterator i =  mesh.internal_comm_list().begin(), iend = mesh.internal_comm_list().end() ; i != iend ; ++i ) {
      if ( (i->owner == parallel_rank && phase == 0) || (i->owner != parallel_rank && phase == 1) ) {
        Entity e = i->entity;
        const MeshIndex meshIdx = mesh.mesh_index(e);
        const unsigned bucketId = meshIdx.bucket->bucket_id();

        for ( fi = fb ; fi != fe ; ++fi ) {
          const FieldBase & f = **fi ;

          if(!is_matching_rank(f, e)) continue;

          const unsigned size = field_bytes_per_entity( f , e );

          if ( size ) {
            unsigned char * ptr =
              reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , bucketId, meshIdx.bucket_ordinal, size ));

            const EntityCommInfoVector& infovec = i->entity_comm->comm_map;
            PairIterEntityComm ec(infovec.begin(), infovec.end());
            if (phase == 0) { // send
              for ( ; !ec.empty() ; ++ec ) {
                if (ec->ghost_id == ghost_id) {
                  CommBuffer & b = sparse.send_buffer( ec->proc );
                  b.pack<unsigned char>( ptr , size );
                }
              }
            }
            else { //recv
              for ( ; !ec.empty(); ++ec ) {
                if (ec->ghost_id == ghost_id) {
                  CommBuffer & b = sparse.recv_buffer( i->owner );
                  b.unpack<unsigned char>( ptr , size );
                  break;
                }
              }
            }
          }
        }
      }
    }
    if (phase == 0) { sparse.communicate(); }
  }
}

//----------------------------------------------------------------------

/** Sum (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */

namespace {

enum Operation
{
  SUM,
  MIN,
  MAX
};

template <typename T, Operation OP>
struct DoOp;

template <typename T>
struct DoOp<T, SUM>
{
  T operator()(T lhs, T rhs) const
  { return lhs + rhs; }
};

template <typename T>
struct DoOp<T, MIN>
{
  T operator()(T lhs, T rhs) const
  { return lhs < rhs ? lhs : rhs; }
};

template <typename T>
struct DoOp<T, MAX>
{
  T operator()(T lhs, T rhs) const
  { return lhs > rhs ? lhs : rhs; }
};

template<typename T>
struct CommMsgs {
  std::vector<int> comm_procs;
  std::vector<bool> msg_recvd;
  std::vector<std::vector<T> > send_data;
  std::vector<std::vector<T> > recv_data;
};

template<typename T, typename MsgPacker, typename MsgUnpacker>
void parallel_data_exchange_sym_pack_unpack(MPI_Comm mpi_communicator,
                                            const std::vector<int>& comm_procs,
                                            MsgPacker& pack_msg,
                                            MsgUnpacker& unpack_msg,
                                            bool deterministic)
{
  //
  //  Determine the number of processors involved in this communication
  //
#if defined( STK_HAS_MPI)
  const int msg_tag = 10242;
  int class_size = sizeof(T);

  int num_comm_procs = comm_procs.size();
  std::vector<std::vector<T> > send_data(num_comm_procs);
  std::vector<std::vector<T> > recv_data(num_comm_procs);
  std::vector<MPI_Request> send_requests(num_comm_procs);
  std::vector<MPI_Request> recv_requests(num_comm_procs);
  std::vector<MPI_Status> statuses(num_comm_procs);

  for(int i=0; i<num_comm_procs; ++i) {
    int iproc = comm_procs[i];
    pack_msg(iproc, send_data[i]);
    recv_data[i].resize(send_data[i].size());

    char* recv_buffer = (char*)recv_data[i].data();
    int buf_size = recv_data[i].size()*class_size;
    MPI_Irecv(recv_buffer, buf_size, MPI_CHAR, iproc, msg_tag, mpi_communicator, &recv_requests[i]);

    char* send_buffer = (char*)send_data[i].data();
    MPI_Isend(send_buffer, buf_size, MPI_CHAR, iproc, msg_tag, mpi_communicator, &send_requests[i]);
  }

  MPI_Status status;
  for(int i = 0; i < num_comm_procs; ++i) {
      int idx = i;
      if (deterministic) {
          MPI_Wait(&recv_requests[i], &status);
      }
      else {
          MPI_Waitany(num_comm_procs, recv_requests.data(), &idx, &status);
      }
      unpack_msg(comm_procs[idx], recv_data[idx]);
  }

  MPI_Waitall(num_comm_procs, send_requests.data(), statuses.data());
#endif
}

template <typename T, Operation OP>
void parallel_op_impl(const BulkData& mesh, std::vector<FieldBase*> fields, bool deterministic = false)
{
  if (fields.empty()) {
    return;
  }

  std::vector<int> comm_procs = mesh.all_sharing_procs(fields[0]->entity_rank());
  stk::mesh::EntityRank first_field_rank = fields[0]->entity_rank();
  for(size_t i=1; i<fields.size(); ++i) {
      const FieldBase* f = fields[i];
      if (f->entity_rank() != first_field_rank) {
          const std::vector<int>& sharing_procs = mesh.all_sharing_procs(f->entity_rank());
          for(int p : sharing_procs) {
              stk::util::insert_keep_sorted_and_unique(p, comm_procs);
          }
      }
  }

  auto msgPacker = [&fields, &mesh](int proc, std::vector<T>& send_data)
  {
    send_data.clear();
    size_t reserve_len = 0;
    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j];
        ThrowRequireMsg(f.type_is<T>(),
                      "Please don't mix fields with different primitive types in the same parallel assemble operation");

        const BucketIndices& bktIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank())[proc];
        for(size_t i=0; i<bktIndices.bucket_info.size(); ++i) {
            unsigned bucket = bktIndices.bucket_info[i].bucket_id;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                reserve_len += bktIndices.bucket_info[i].num_entities_this_bucket*num_Ts_per_field;
            }
        }
    }
    send_data.reserve(reserve_len);

    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j];
        const BucketIndices& bktIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank())[proc];

        const std::vector<unsigned>& ords = bktIndices.ords;
        unsigned offset = 0;

        for(size_t i=0; i<bktIndices.bucket_info.size(); ++i) {
            unsigned bucket = bktIndices.bucket_info[i].bucket_id;
            unsigned num_entities = bktIndices.bucket_info[i].num_entities_this_bucket;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                const T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket));

                for (unsigned iord=0; iord<num_entities; ++iord) {
                    const unsigned idx = ords[offset++]*num_Ts_per_field;

                    for (int d = 0; d < num_Ts_per_field; ++d) {
                      send_data.push_back(data[idx+d]);
                    }
                }
            }
            else {
                offset += num_entities;
            }
        }
    }
  };

  auto msgUnpacker = [&fields, &mesh](int iproc, std::vector<T>& recv_data)
  {
    DoOp<T, OP> do_op;

    unsigned offset = 0;
    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j] ;
        const BucketIndices& bktIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank())[iproc];
        const std::vector<unsigned>& ords = bktIndices.ords;

        size_t ords_offset = 0;
        for(size_t i=0; i<bktIndices.bucket_info.size(); ++i) {
            unsigned bucket = bktIndices.bucket_info[i].bucket_id;
            unsigned num_entities = bktIndices.bucket_info[i].num_entities_this_bucket;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket));

                for (unsigned iord=0; iord<num_entities; ++iord) {
                    const unsigned idx = ords[ords_offset++]*num_Ts_per_field;

                    for (int d = 0; d < num_Ts_per_field; ++d) {
                      data[idx+d] = do_op(data[idx+d], recv_data[offset + d]);
                    }
                    offset += num_Ts_per_field;
                }
            }
            else {
                ords_offset += num_entities;
            }
        }
    }
  };

  MPI_Comm comm = mesh.parallel();
  parallel_data_exchange_sym_pack_unpack<T>(comm, comm_procs, msgPacker, msgUnpacker, deterministic);
}

template <Operation OP>
inline
void parallel_op(const BulkData& mesh, const std::vector<FieldBase*>& fields, bool deterministic)
{
  if (mesh.parallel_size() == 1 || fields.empty()) return;

  if (fields[0]->type_is<double>()) {
    parallel_op_impl<double, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<float>()) {
    parallel_op_impl<float, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<int>()) {
    parallel_op_impl<int, OP>(mesh, fields, deterministic);
  }
  else {
    ThrowRequireMsg(false, "Error, parallel_max only operates on fields of type double, float or int.");
  }
}

}

void parallel_sum(const BulkData& mesh, const std::vector<FieldBase*>& fields, bool deterministic)
{
  parallel_op<SUM>(mesh, fields, deterministic);
}

//----------------------------------------------------------------------

/** Communicate and take the maximum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (maximum) field values on each sharing proc.
 */
void parallel_max(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  parallel_op<MAX>(mesh, fields, false);
}

/** Communicate and take the minimum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (minimum) field values on each sharing proc.
 */
void parallel_min(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  parallel_op<MIN>(mesh, fields, false);
}

} // namespace mesh
} // namespace stk
