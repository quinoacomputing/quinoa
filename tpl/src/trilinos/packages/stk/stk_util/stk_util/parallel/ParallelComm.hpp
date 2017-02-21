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

#ifndef stk_util_parallel_ParallelComm_hpp
#define stk_util_parallel_ParallelComm_hpp

#include <cstddef>                      // for size_t, ptrdiff_t
#include <vector>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/environment/ReportHandler.hpp> // for ThrowAssertMsg

namespace stk { template <unsigned int N> struct CommBufferAlign; }

//------------------------------------------------------------------------

namespace stk {

/** Perform collective all-to-all communication with individually
 *  varying message sizes.  The collective operation uses an
 *  all-to-all if the maximum number of sends or receives from
 *  any one processor is greater than the given bounds.
 *
 *  This is a work-horse communication for mesh data structures
 *  with parallel domain decomposition.
 */
class CommAll ;
class CommSparse;

/** Pack and unpack buffers for the sparse all-to-all communication.
 */
class CommBuffer ;

/** Given the send sizes determine the receive sizes.
 *  Send and receive size arrays are dimensioned to
 *  the size of the parallel machine.
 *  Return global parallel logical OR of the input local flag.
 *  This parallel reduction is aggregated into the required
 *  communication for determining the sparse sizes.
 *  Output the receive sizes and maximum number of send or
 *  receive messages for a single processor.
 *  A dense all-to-all communication is used if:
 *     num_msg_bound < num_msg_maximum
 *  otherwise a set of point-to-point messages are used.
 */
bool comm_sizes( ParallelMachine ,
                 const unsigned   num_msg_bound ,
                       unsigned & num_msg_maximum ,
                 const unsigned * const send_size ,
                       unsigned * const recv_size ,
                 bool local_flag = false );

/** If the communication is known to be dense.
 */
bool comm_dense_sizes( ParallelMachine ,
                       const unsigned * const send_size ,
                             unsigned * const recv_size ,
                       bool local_flag = false );

//------------------------------------------------------------------------

class CommBuffer {
public:

  /** Pack a value to be sent:  buf.pack<type>( value ) */
  template<typename T> CommBuffer &pack( const T & value );

private:
  /** Do not try to pack a pointer for global communication */
  template<typename T> CommBuffer &pack( const T* value ) {
    ThrowAssertMsg(false,"CommBuffer::pack(const T* value) not allowed. Don't pack a pointer for communication!");
    return *this;
  }

public:

  /** Pack an array of values to be sent:  buf.pack<type>( ptr , num ) */
  template<typename T> CommBuffer &pack( const T * value , size_t number );

  /** Unpack a received value:  buf.unpack<type>( value ) */
  template<typename T> CommBuffer &unpack( T & value );

  /** Unpack an array of received values:  buf.unpack<type>( ptr , num ) */
  template<typename T> CommBuffer &unpack( T * value , size_t number );

  /** Peek at a received value (don't advance buffer): buf.peek<type>(value) */
  template<typename T> CommBuffer &peek( T & value );

  /** Peek at an array of received values: buf.peek<type>( ptr , num ) */
  template<typename T> CommBuffer &peek( T * value , size_t number );

  /** Skip buffer ahead by a number of values. */
  template<typename T> CommBuffer &skip( size_t number );

  /** Reset the buffer to the beginning so that size() == 0 */
  void reset();

  /** Reset the buffer pointers to NULL */
  void reset_to_null();

  /** Size, in bytes, of the buffer.
   *  If the buffer is not yet allocated this is zero.
   */
  size_t capacity() const ;

  // TODO - terribly misinforming when used on recv buffer, returns 0!
  /** Size, in bytes, of the buffer that has been processed.
   *  If the buffer is not yet allocated then this is the
   *  number of bytes that has been attempted to pack.
   */
  size_t size() const ;
  void set_size(size_t newsize_bytes);

  /** Size, in bytes, of the buffer remaining to be processed.
   *  Equal to 'capacity() - size()'.  A negative result
   *  indicates either the buffer is not allocated or an
   *  overflow has occurred.  An overflow will have thrown
   *  an exception.
   */
  ptrdiff_t remaining() const ;

  /** Pointer to base of buffer. */
  void * buffer() const ;

  ~CommBuffer();
  CommBuffer();

private:
  friend class CommAll ;
  friend class CommSparse ;
  friend class CommGather ;
  friend class CommBroadcast ;

  static CommBuffer * allocate( const unsigned, const unsigned * const );
  static void deallocate( const unsigned , CommBuffer * );

  void pack_overflow() const ;
  void unpack_overflow() const ;

  typedef unsigned char * ucharp ;

  ucharp m_beg ;
  ucharp m_ptr ;
  ucharp m_end ;
};

//------------------------------------------------------------------------

class CommAll {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  int             parallel_size() const { return m_size ; }
  int             parallel_rank() const { return m_rank ; }

  /** Obtain the message buffer for a given processor */
  CommBuffer & send_buffer( int p ) const
  {
#ifndef NDEBUG
    if ( m_size <= p ) { rank_error("send_buffer",p); }
#endif
    return m_send[p] ;
  }

  /** Obtain the message buffer for a given processor */
  CommBuffer & recv_buffer( int p ) const
  {
#ifndef NDEBUG
    if ( m_size <= p ) { rank_error("recv_buffer",p); }
#endif
    return m_recv[p] ;
  }

  //----------------------------------------
  /** Construct for undefined communication.
   *  No buffers are allocated.
   */
  CommAll(bool propagate_local_error_flags=true);

  /** Allocate send and receive buffers based upon input sizes.
   *  If recv_size == NULL then the receive size
   *  is determined by communicating the send sizes.
   *  Symmetry is given by passing the same data for both
   *  send and receive sizes.
   *  Return global parallel OR of local flags.
   */
  bool allocate_buffers( ParallelMachine ,
                         const unsigned num_msg_bounds ,
                         const unsigned * const send_size ,
                         const unsigned * const recv_size ,
                         const bool local_flag = false );

  bool allocate_buffers( ParallelMachine comm ,
                         const unsigned * const send_sizes,
                         const unsigned * const recv_sizes );

  /**
   * Allocate symmetric buffers, no communication required. buf_sizes should
   * have lenth = parallel_size(comm).
   */
  bool allocate_symmetric_buffers( ParallelMachine comm ,
                                   const unsigned * const buf_sizes );

  //----------------------------------------
  /** Construct for a to-be-sized communication.
   *  Allocate surrogate send buffers to enable
   *  no-op packing for the purpose of send sizing.
   *  Surrogate send scenario:
   *  1) Surrogate send buffers are "packed" for sizing where
   *     packing sizes are recorded but no data is copied.
   *  2) 'allocate_buffers(symmetric,flag)' is called to allocate
   *     buffers.  The symmetric flag guarantees that the send
   *     sizes matches the receive sizes.
   *  3) Send buffers are identically packed; however, this
   *     packing copies data into the send buffers.
   */
  explicit CommAll( ParallelMachine, bool propagate_local_error_flags=true );

  /** Allocate asymmetric communication based upon
   *  sizing from the surrogate send buffer packing.
   *  If symmetric then the receive sizes are guaranteed
   *  to be identical to the send sizes.
   *  Return global parallel OR of local flags.
   */
  bool allocate_buffers( const unsigned num_msg_bounds ,
                         const bool symmetric = false ,
                         const bool local_flag = false );

  //----------------------------------------
  /** Communicate send buffers to receive buffers.  */
  void communicate();

  //----------------------------------------
  /** Swap send and receive buffers leading to reversed communication. */
  void swap_send_recv();

  /** Reset, but do not reallocate, message buffers for reprocessing.
   *  Sets 'size() == 0' and 'remaining() == capacity()'.
   */
  void reset_buffers();

  ~CommAll();

  static bool sm_verbose;
  static bool sm_verbose_proc0_only;

private:

  CommAll( const CommAll & );
  CommAll & operator = ( const CommAll & );

  void rank_error( const char * , int ) const ;

  bool allocate_buffers( const unsigned * const send_size ,
                         const unsigned * const recv_size ,
                         bool local_flag );

  ParallelMachine m_comm ;
  bool            m_propagate_local_error_flags;
  int             m_size ;
  int             m_rank ;
  unsigned        m_bound ;
  unsigned        m_max ;
  CommBuffer    * m_send ;
  CommBuffer    * m_recv ;
};

//------------------------------------------------------------------------

class CommBroadcast {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  int             parallel_size() const { return m_size ; }
  int             parallel_rank() const { return m_rank ; }

  /** Obtain the message buffer for the root_rank processor */
  CommBuffer & send_buffer();

  /** Obtain the message buffer for the local processor */
  CommBuffer & recv_buffer();

  //----------------------------------------

  CommBroadcast( ParallelMachine , int root_rank );

  void communicate();

  bool allocate_buffer( const bool local_flag = false );

  ~CommBroadcast();

private:

  CommBroadcast();
  CommBroadcast( const CommBroadcast & );
  CommBroadcast & operator = ( const CommBroadcast & );

  ParallelMachine m_comm ;
  int             m_size ;
  int             m_rank ;
  int             m_root_rank ;
  CommBuffer      m_buffer ;
};

//----------------------------------------------------------------------

class CommGather {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  int             parallel_size() const { return m_size ; }
  int             parallel_rank() const { return m_rank ; }

  ~CommGather();

  CommGather( ParallelMachine , int root_rank , unsigned send_size );

  CommBuffer & send_buffer() { return m_send ; }

  void communicate();

  CommBuffer & recv_buffer( int );

  void reset();

private:

  CommGather();
  CommGather( const CommGather & );
  CommGather & operator = ( const CommGather & );

  ParallelMachine m_comm ;
  int             m_size ;
  int             m_rank ;
  int             m_root_rank ;
  CommBuffer      m_send ;
  CommBuffer    * m_recv ;
  int           * m_recv_count ;
  int           * m_recv_displ ;
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Inlined template implementations for the CommBuffer

namespace stk {


template<>
struct CommBufferAlign<1> {
  static size_t align( size_t ) { return 0 ; }
};

template<unsigned N>
struct CommBufferAlign {
  static size_t align( size_t i ) { i %= N ; return i ? ( N - i ) : 0 ; }
};

template<typename T>
inline
CommBuffer &CommBuffer::pack( const T & value )
{
  enum { Size = sizeof(T) };
  size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  if ( m_beg ) {
    if ( m_end < m_ptr + nalign + Size ) { pack_overflow(); }
    while ( nalign ) { --nalign ; *m_ptr = 0 ; ++m_ptr ; }
    T * tmp = reinterpret_cast<T*>(m_ptr);
    *tmp = value ;
    m_ptr = reinterpret_cast<ucharp>( ++tmp );
  }
  else {
    m_ptr += nalign + Size ;
  }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::pack( const T * value , size_t number )
{
  enum { Size = sizeof(T) };
  size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  if ( m_beg ) {
    if ( m_end < m_ptr + nalign + number * Size ) { pack_overflow(); }
    while ( nalign ) { --nalign ; *m_ptr = 0 ; ++m_ptr ; }
    T * tmp = reinterpret_cast<T*>(m_ptr);
    while ( number ) { --number ; *tmp = *value ; ++tmp ; ++value ; }
    m_ptr = reinterpret_cast<ucharp>( tmp );
  }
  else {
    m_ptr += nalign + number * Size ;
  }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::skip( size_t number )
{
  enum { Size = sizeof(T) };
  m_ptr += CommBufferAlign<Size>::align( m_ptr - m_beg ) + Size * number ;
  if ( m_beg && m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::unpack( T & value )
{
  enum { Size = sizeof(T) };
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  value = *tmp ;
  m_ptr = reinterpret_cast<ucharp>( ++tmp );
  if ( m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::unpack( T * value , size_t number )
{
  enum { Size = sizeof(T) };
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  while ( number ) { --number ; *value = *tmp ; ++tmp ; ++value ; }
  m_ptr = reinterpret_cast<ucharp>( tmp );
  if ( m_end < m_ptr ) { unpack_overflow(); }
  return *this;
}
template<typename item>
inline
item unpack(stk::CommBuffer& buf)
{
    item object;
    buf.unpack<item>(object);
    return object;
}

template<typename T>
inline
CommBuffer &CommBuffer::peek( T & value )
{
  enum { Size = sizeof(T) };
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  value = *tmp ;
  if ( m_end < reinterpret_cast<ucharp>(++tmp) ) { unpack_overflow(); }
  return *this;
}

template<typename T>
inline
CommBuffer &CommBuffer::peek( T * value , size_t number )
{
  enum { Size = sizeof(T) };
  const size_t nalign = CommBufferAlign<Size>::align( m_ptr - m_beg );
  T * tmp = reinterpret_cast<T*>( m_ptr + nalign );
  while ( number ) { --number ; *value = *tmp ; ++tmp ; ++value ; }
  if ( m_end < reinterpret_cast<ucharp>(tmp) ) { unpack_overflow(); }
  return *this;
}

inline
void CommBuffer::reset()
{ m_ptr = m_beg ; }

inline
void CommBuffer::reset_to_null()
{ m_beg = nullptr; m_ptr = nullptr; m_end = nullptr; }

inline
size_t CommBuffer::capacity() const
{ return m_end - m_beg ; }

inline
size_t CommBuffer::size() const
{ return m_ptr - m_beg ; }

inline
void CommBuffer::set_size(size_t newsize_bytes)
{ m_beg = nullptr;  m_ptr = nullptr; m_ptr += newsize_bytes ; m_end = nullptr; }

inline
ptrdiff_t CommBuffer::remaining() const
{ return m_end - m_ptr ; }

inline
void * CommBuffer::buffer() const
{ return static_cast<void*>( m_beg ); }

std::vector<int> ComputeReceiveList(std::vector<int>& sendSizeArray, MPI_Comm &mpi_communicator);

//
//  Parallel_Data_Exchange: General object exchange template with unknown comm plan
//
template<typename T>
void parallel_data_exchange_t(std::vector< std::vector<T> > &send_lists,
                              std::vector< std::vector<T> > &recv_lists,
                              MPI_Comm &mpi_communicator ) {
  //
  //  Determine the number of processors involved in this communication
  //
  const int msg_tag = 10242;
  int num_procs;
  MPI_Comm_size(mpi_communicator, &num_procs);
  int my_proc;
  MPI_Comm_rank(mpi_communicator, &my_proc);
  ThrowRequire((unsigned int) num_procs == send_lists.size() && (unsigned int) num_procs == recv_lists.size());
  int class_size = sizeof(T);
  //
  //  Determine number of items each other processor will send to the current processor
  //
  std::vector<int> global_number_to_send(num_procs);
  for(int iproc=0; iproc<num_procs; ++iproc) {
    global_number_to_send[iproc] = send_lists[iproc].size();
  }
  std::vector<int> numToRecvFrom = ComputeReceiveList(global_number_to_send, mpi_communicator);
  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    recv_lists[iproc].resize(numToRecvFrom[iproc]);
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)recv_lists[iproc].data();
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR, iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)send_lists[iproc].data();
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR,
               iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
}

//
//  Generalized comm plans
//
//  This plan assumes the send and recv lists have identical sizes so no extra sizing communications are needed
//
template<typename T>
void parallel_data_exchange_sym_t(std::vector< std::vector<T> > &send_lists,
                                  std::vector< std::vector<T> > &recv_lists,
                                  MPI_Comm &mpi_communicator )
{
  //
  //  Determine the number of processors involved in this communication
  //
#if defined( STK_HAS_MPI)
  const int msg_tag = 10242;
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    recv_lists[iproc].resize(send_lists[iproc].size());
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)recv_lists[iproc].data();
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)send_lists[iproc].data();
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR,
               iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
#endif
}

template<typename T>
inline void parallel_data_exchange_nonsym_known_sizes_t(std::vector< std::vector<T> > &send_lists,
                                                std::vector< std::vector<T> > &recv_lists,
                                                MPI_Comm mpi_communicator )
{
#if defined( STK_HAS_MPI)
  const int msg_tag = 10243; //arbitrary tag value, anything less than 32768 is legal
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)recv_lists[iproc].data();
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR, iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }

  MPI_Barrier(mpi_communicator);

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)send_lists[iproc].data();
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR, iproc, msg_tag, mpi_communicator);
    }
  }

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
#endif
}

//
//  This plan assumes the send and recv lists are matched, but that the actual ammount of data to send is unknown.
//  A processor knows which other processors it will be receiving data from, but does not know who much data.
//  Thus the comm plan is known from the inputs, but an additional message sizing call must be done.
//
template<typename T>
void parallel_data_exchange_sym_unknown_size_t(std::vector< std::vector<T> > &send_lists,
                                               std::vector< std::vector<T> > &recv_lists,
                                               MPI_Comm &mpi_communicator )
{
#if defined( STK_HAS_MPI)
  const int msg_tag = 10242;
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the message sizes
  //
  std::vector<int> send_msg_sizes(num_procs);
  std::vector<int> recv_msg_sizes(num_procs);
  std::vector<MPI_Request> recv_handles(num_procs);

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    send_msg_sizes[iproc] = send_lists[iproc].size();
  }    
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size()>0) {
      MPI_Irecv(&recv_msg_sizes[iproc], 1, MPI_INT, iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size()>0) {
      MPI_Send(&send_msg_sizes[iproc], 1, MPI_INT, iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
      recv_lists[iproc].resize(recv_msg_sizes[iproc]);
    }
  }
  //
  //  Send the actual messages as raw byte streams.
  //
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)recv_lists[iproc].data();
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)send_lists[iproc].data();
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR,
               iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
#endif
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

