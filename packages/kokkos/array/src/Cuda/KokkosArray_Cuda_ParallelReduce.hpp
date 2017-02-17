/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CUDA_PARALLELREDUCE_HPP
#define KOKKOS_CUDA_PARALLELREDUCE_HPP

#if defined( __CUDACC__ )

#include <KokkosArray_ParallelReduce.hpp>

#include <vector>
#include <iostream>
#include <stdexcept>

#include <Cuda/KokkosArray_Cuda_Parallel.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------
// Maximize shared memory and minimize L1 cache:
//   cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferShared );
// For 2.0 capability: 48 KB shared and 16 KB L1
//----------------------------------------------------------------------------
// Must have consistent '__shared__' statement across all device kernels.
// Since there may be more than one kernel in a file then have to make this
// a simple array of words.
//----------------------------------------------------------------------------

template< class DriverType >
__device__ __noinline__
void cuda_reduce_shared( const Cuda::size_type used_warp_count )
{
  typedef          Cuda::size_type  size_type ;
  typedef typename DriverType::value_type value_type ;

  typedef volatile value_type * vvp ;
  typedef volatile const value_type * cvvp ;

  enum { ValueWordStride = DriverType::ValueWordStride };
  enum { WarpStride      = DriverType::WarpStride };
  enum { HalfWarpSize    = Impl::CudaTraits::WarpSize >> 1 };

  extern __shared__ size_type shared_data[];

  // threadIdx.x == index within warp [ 0 .. WarpSize - 1 ]
  // threadIdx.y == which warp        [ 0 .. used_warp_count - 1 ]

  // Phase A: Reduce within my warp:
  //          Warp's reads occur before joins and writes
  //          so there is no race condition.
  //          Declare shared data to be volatile to
  //          prevent compiler from introducing a race condition.
  //
  if ( threadIdx.y < used_warp_count && threadIdx.x < HalfWarpSize ) {
    enum { n1  = ValueWordStride * 1 };
    enum { n2  = ValueWordStride * 2 };
    enum { n4  = ValueWordStride * 4 };
    enum { n8  = ValueWordStride * 8 };
    enum { n16 = ValueWordStride * 16 };

    size_type * const data = shared_data + DriverType::shared_data_offset();

    DriverType::join( *((vvp) data), *((cvvp)( data + n16 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n8 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n4 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n2 )) );
    DriverType::join( *((vvp) data), *((cvvp)( data +  n1 )) );
  }

  // Phase B: Use a single warp to reduce results from each warp.
  //          This requires: used_warp_count <= WarpSize
  //

  __syncthreads();

  if ( 0 == threadIdx.y && threadIdx.x + 1 < used_warp_count ) {
    enum { n1  = WarpStride * ValueWordStride * 1 };
    enum { n2  = WarpStride * ValueWordStride * 2 };
    enum { n4  = WarpStride * ValueWordStride * 4 };
    enum { n8  = WarpStride * ValueWordStride * 8 };
    enum { n16 = WarpStride * ValueWordStride * 16 };

    size_type * const data = shared_data + DriverType::shared_data_offset( 0 , threadIdx.x );

    if ( threadIdx.x + 2 < used_warp_count ) {
      if ( threadIdx.x + 4 < used_warp_count ) {
        if ( threadIdx.x + 8 < used_warp_count ) {
          if ( threadIdx.x + 16 < used_warp_count ) {
            DriverType::join( *((vvp) data) , *((cvvp)( data + n16 )) );
          }
          DriverType::join( *((vvp) data) , *((cvvp)( data + n8 )) );
        }
        DriverType::join( *((vvp) data) , *((cvvp)( data + n4 )) );
      }
      DriverType::join( *((vvp) data) , *((cvvp)( data + n2 )) );
    }
    DriverType::join( *((vvp) data) , *((cvvp)( data + n1 )) );
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class ReduceTraits , class FinalizeType >
class ParallelReduce< FunctorType , ReduceTraits , FinalizeType , Cuda > {
public:
  typedef ParallelReduce< FunctorType , ReduceTraits , FinalizeType , Cuda > this_type ;
  typedef          Cuda               device_type ;
  typedef          Cuda::size_type    size_type ;
  typedef typename ReduceTraits::value_type value_type ;

  enum { WarpSize   = Impl::CudaTraits::WarpSize };
  enum { WarpStride = WarpSize + 1 };
  enum { WarpIndexMask  = Impl::CudaTraits::WarpIndexMask };
  enum { WarpIndexShift = Impl::CudaTraits::WarpIndexShift };

  enum { SharedMemoryBanks = Impl::CudaTraits::SharedMemoryBanks_13 };

  enum { ValueWordCount = ( sizeof(value_type) + sizeof(size_type) - 1 )
                          / sizeof(size_type) };

  /** \brief  If the reduction value occupies an
   *          exact multiple of shared memory banks
   *          then it must be padded to avoid bank conflicts.
   */
  enum { ValueWordStride = ValueWordCount +
          ( ValueWordCount % SharedMemoryBanks ? 0 : 2 ) };

  enum { WordsPerWarp       = ValueWordStride * WarpSize };
  enum { WordsPerWarpStride = ValueWordStride * WarpStride };

  //----------------------------------------------------------------------

  const FunctorType  m_work_functor ;
  const FinalizeType m_work_finalize ;
  const size_type    m_work_count ;
  const size_type    m_work_stride ;

  const size_type    m_work_block_offset ;
  const size_type    m_global_block_count ;
  const size_type    m_stream_block_offset ;
  const size_type    m_stream_block_recycle ;

  // Scratch space for multi-block reduction
  // m_scratch_warp  == number of warps required
  // m_scratch_upper == upper bound of reduction scratch space used.

  size_type * const m_scratch_space ;
  size_type * const m_scratch_flag ;
  const size_type   m_scratch_warp ;
  const size_type   m_scratch_upper ;
  
  //----------------------------------------------------------------------

  static inline
  __device__
  void join( volatile       value_type & update ,
             volatile const value_type & input )
    { ReduceTraits::join( update , input ); }

  static inline
  __device__
  void init( value_type & update )
    { ReduceTraits::init( update ); }

  static inline
  __device__
  size_type shared_data_offset()
  { return ValueWordStride * ( threadIdx.x + WarpStride * threadIdx.y ); }

  static inline
  __device__
  size_type shared_data_offset( size_type x , size_type y )
  { return ValueWordStride * ( x + WarpStride * y ); }

  static inline
  __device__
  size_type shared_flag_offset()
  { return ValueWordStride * ( WarpStride * blockDim.y - 1 ); }

  inline
  __device__
  size_type * scratch_data_for_block() const
  {
    const size_type stream_block_id = m_stream_block_offset + blockIdx.x ;

    return m_scratch_space +
      shared_data_offset(
        stream_block_id &  WarpIndexMask  /* for threadIdx.x */ ,
        stream_block_id >> WarpIndexShift /* for threadIdx.y */ );
  }


public:

  ParallelReduce( const FunctorType  & functor ,
                  const FinalizeType & finalize ,
                  const size_type      work_count ,
                  const size_type      work_stride ,
                  const size_type      work_block_offset ,
                  const size_type      global_block_count ,
                  const size_type      stream_block_count ,
                  const size_type      stream_block_offset ,
                  const size_type      stream_block_recycle )
    : m_work_functor(  functor )
    , m_work_finalize( finalize )
    , m_work_count(    work_count )
    , m_work_stride(   work_stride )
    , m_work_block_offset( work_block_offset )
    , m_global_block_count(   global_block_count )
    , m_stream_block_offset(  stream_block_offset )
    , m_stream_block_recycle( stream_block_recycle )
    , m_scratch_space( cuda_internal_reduce_multiblock_scratch_space() )
    , m_scratch_flag(  cuda_internal_reduce_multiblock_scratch_flag() )
    , m_scratch_warp( ( stream_block_count >> WarpIndexShift ) +
                      ( stream_block_count &  WarpIndexMask ? 1 : 0 ) )
    , m_scratch_upper( ValueWordStride * ( stream_block_count + m_scratch_warp - 1 ) )
  {}

public:

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // Warp #0 of this block read this block's previous reduction value
  // from global memory.
  // This will be the initial reduction value for Warp #0 Thread #0.
  inline
  __device__
  void reduce_global_recycle() const
  {
    enum { WarpSize = Impl::CudaTraits::WarpSize };

    extern __shared__ size_type shared_data[];

    if ( 0 == threadIdx.y ) {
      // Coalesced global memory read

      size_type * const scratch = scratch_data_for_block();

      for ( size_type i = threadIdx.x ; i < ValueWordCount ; i += WarpSize ) {
        shared_data[i] = scratch[i] ;
      }
    }
  }

  //--------------------------------------------------------------------------
  // The last block to finish reduces the results from all blocks.
  inline
  __device__
  void reduce_global_complete() const
  {
    enum { WarpSize = Impl::CudaTraits::WarpSize };

    extern __shared__ size_type shared_data[];

    // Each warp does a coalesced read of its own data.

    if ( threadIdx.y < m_scratch_warp ) {

      // Coalesced global memory read for this warp's data.

      size_type i = shared_data_offset( 0 , threadIdx.y );
      size_type j = i + WordsPerWarp ;

      if ( m_scratch_upper < j ) {
        j = m_scratch_upper ;

        // Only partial data will be read by this warp
        // so initialize the values before reading.
        init( *((value_type *)( shared_data + shared_data_offset() )) );
      }

      for ( i += threadIdx.x ; i < j ; i += WarpSize ) {
        shared_data[i] = m_scratch_space[i] ;
      }
    }

    // Reduce these contributions
    cuda_reduce_shared< this_type >( m_scratch_warp );
  }

  //--------------------------------------------------------------------------
  // Warp #0 of this block write this block's reduction value to global memory

  __device__
  size_type reduce_global_contribute() const
  {
    enum { WarpSize = Impl::CudaTraits::WarpSize };

    extern __shared__ size_type shared_data[];

    const size_type flag_offset = shared_flag_offset();

    if ( 0 == threadIdx.y ) {

      // Coalesced global memory write

      size_type * const scratch = scratch_data_for_block();

      for ( size_type i = threadIdx.x ; i < ValueWordCount ; i += WarpSize ) {
        scratch[i] = shared_data[i] ;
      }

      __threadfence(); // Wait for write to complete

      //----------------------------------
      // Warp #0 Thread #0 : Check if this is the last block to finish:

      if ( 0 == threadIdx.x ) {

        // atomicInc returns value prior to increment.

        const size_type last_block_flag =
          m_global_block_count ==
          1 + atomicInc(m_scratch_flag, m_global_block_count+1);

        // Inform the entire block of the last_block status:
        shared_data[ flag_offset ] = last_block_flag ;

        if ( last_block_flag ) {
          // Reset the flag for the next reduction operation
          *m_scratch_flag = 0 ;
        }
      }
    }

    __syncthreads(); // All threads of block wait for last_block flag to be set.

    return shared_data[ flag_offset ];
  }
  //--------------------------------------------------------------------------

  inline
  __device__
  void execute_on_device() const
  {
    typedef Cuda::size_type size_type ;

    extern __shared__ size_type shared_data[];

    value_type & value = *( (value_type *)( shared_data + shared_data_offset() ) );
    init( value );

    // If this block is not the first block in the stream
    // then must read the previous value into
    // the to-be-updated thread #0 value.

    if ( m_stream_block_recycle ) {
      reduce_global_recycle();
    }

    // Phase 1: Reduce to per-thread contributions
    {
      const size_type work_count  = m_work_count ;
      const size_type work_stride = m_work_stride ;

      size_type iwork =
        threadIdx.x + blockDim.x * (
        threadIdx.y + blockDim.y * (
        blockIdx.x + m_work_block_offset ) );

      for ( ; iwork < work_count ; iwork += work_stride ) {
        m_work_functor( iwork , value );
      }
    }

    // Phase 2: Reduce this block's thread's contributions
    //          to a single reduction value.
    cuda_reduce_shared< this_type >( blockDim.y );

    // Phase 3: Reduce contributions from multiple blocks

    int last_block = 1 == m_global_block_count ;

    if ( ! last_block ) {

      // This block contribute its partial result
      // and determine if it is the last block.

      last_block = reduce_global_contribute();

      // If the last block then reduce partial result
      // contributions from all blocks to a single result.
      if ( last_block ) {

        reduce_global_complete();
      }
    }

    // Phase 4: Last block Warp #0 Thread #0 performs serial finalization
    if ( last_block && 0 == threadIdx.x && 0 == threadIdx.y ) {
      m_work_finalize( value );
    }
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  static
  void execute( const size_t         work_count ,
                const FunctorType  & functor ,
                const FinalizeType & finalize )
  {
    typedef MemoryManager< Cuda > memory_manager ;

    typedef ParallelReduce< FunctorType , ReduceTraits , FinalizeType , Cuda > self_type ;

    const size_type maximum_shared_words = cuda_internal_maximum_shared_words();

    dim3 block( Impl::CudaTraits::WarpSize , 
                cuda_internal_maximum_warp_count() , 1 );

    while ( maximum_shared_words < block.y * WordsPerWarpStride ) {
      block.y >>= 1 ;
    }

    dim3 grid( 1 , 1 , 1 );

    if ( 0 == work_count ) { // Avoid infinite loop.
      block.y = 1 ;
    }
    else if ( work_count <= WarpSize * block.y ) { // Need at most one block
      while ( work_count <= WarpSize * ( block.y >> 1 ) ) { block.y >>= 1 ; }
    }
    else {
      const size_t threads_per_block = WarpSize * block.y ;

      grid.x = ( work_count + threads_per_block - 1 ) / threads_per_block ;

      // At most one block per thread so that the final reduction
      // operation can process one reduction value per thread.
      if ( grid.x > threads_per_block ) { grid.x = threads_per_block ; }
    }

    const size_type work_stride = block.x * block.y * grid.x ;
    const size_type shmem_size =
      sizeof(size_type) * ( ValueWordStride * (WarpStride * block.y - 1) + 1 );

    memory_manager::disable_memory_view_tracking();

    self_type driver( functor , finalize ,
                      work_count , work_stride , 0 ,
                      grid.x , grid.x , 0 , 0 );

    memory_manager::enable_memory_view_tracking();

    CudaParallelLaunch< self_type >::execute( driver , grid , block , shmem_size );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class FunctorType , class ReduceTraits >
class ParallelReduce< FunctorType , ReduceTraits , void , Cuda > 
{
public:
  typedef typename ReduceTraits::value_type     value_type ;
  typedef Value< value_type , Cuda >  view_type ;

  static
  void execute( const size_t        work_count ,
                const FunctorType & work_functor ,
                      value_type  & result )
  {
    view_type tmp =
      create_value< view_type >(
        std::string("parallel_reduce_temporary_result") );

    ParallelReduce< FunctorType , ReduceTraits , view_type , Cuda >
      ::execute( work_count , work_functor , tmp );

    deep_copy( result , tmp );
  }
};


} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if 1

namespace KokkosArray {
namespace Impl {

template < class FunctorType , class ReduceTraits , class FinalizeType >
class CudaMultiFunctorParallelReduceMember ;

template < class ReduceTraits , class FinalizeType >
class CudaMultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType> {
public:

  typedef Cuda::size_type size_type ;

protected:
  CudaMultiFunctorParallelReduceMember( size_type stream_count )
    : m_stream_count( stream_count )
    {}
public:

  const size_type m_stream_count ; ///< Number of streams needed

  virtual ~CudaMultiFunctorParallelReduceMember() {}

  virtual
  void execute( const FinalizeType & finalize ,
                const size_type      block_count ,
                const size_type      block_offset ,
                const size_type      warp_count ,
                const size_type      shmem_size ,
                cudaStream_t       & cuda_stream ,
                const size_type      global_block_count ,
                const size_type      stream_block_count ,
                const size_type      stream_block_offset ,
                const size_type      stream_block_recycle ) const = 0 ;
};

template < class FunctorType , class ReduceTraits , class FinalizeType >
class CudaMultiFunctorParallelReduceMember :
  public CudaMultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType>
{
public:
  typedef CudaMultiFunctorParallelReduceMember<void,ReduceTraits,FinalizeType> base_type ;
  
  typedef Cuda            device_type ;
  typedef Cuda::size_type size_type ;

  FunctorType m_functor ;
  size_type   m_work_count ;   ///< Parallel work items

  CudaMultiFunctorParallelReduceMember(
    const FunctorType & functor ,
    const size_type     work_count ,
    const size_type     stream_count )
  : base_type( stream_count )
  , m_functor( functor )
  , m_work_count( work_count )
  {}

  virtual
  void execute( const FinalizeType & finalize ,
                const size_type      block_count ,
                const size_type      block_offset ,
                const size_type      warp_count ,
                const size_type      shmem_size ,
                cudaStream_t       & cuda_stream ,
                const size_type      global_block_count ,
                const size_type      stream_block_count ,
                const size_type      stream_block_offset ,
                const size_type      stream_block_recycle ) const
  {
    typedef MemoryManager< Cuda > memory_manager ;

    typedef ParallelReduce< FunctorType , ReduceTraits , FinalizeType , Cuda >  driver_type ;

    const dim3 grid( block_count , 1 , 1 );
    const dim3 block( Impl::CudaTraits::WarpSize , warp_count , 1 );
    const size_type work_stride =
      block.x * block.y * grid.x * base_type::m_stream_count ;

    memory_manager::disable_memory_view_tracking();

    driver_type driver( m_functor , finalize ,
                        m_work_count , work_stride , block_offset ,
                        global_block_count ,
                        std::min( global_block_count , stream_block_count ),
                        stream_block_offset , stream_block_recycle );

    memory_manager::enable_memory_view_tracking();

    // Currently must be local memory launch for
    // multiple kernel launches into multiple streams.
    // For constant memory utilization each copy of a functors into
    // constant memory would require a unique and managed location.

    cuda_parallel_launch_local_memory< driver_type ><<< grid , block , shmem_size , cuda_stream >>>( driver );
  }
};

} // namespace Impl

template < class ReduceTraits , class FinalizeType >
class MultiFunctorParallelReduce< ReduceTraits , FinalizeType , Cuda > {
public:
  typedef          Cuda               device_type ;
  typedef          Cuda::size_type    size_type ;
  typedef typename ReduceTraits::value_type value_type ;

  enum { WarpSize   = Impl::CudaTraits::WarpSize };
  enum { WarpStride = WarpSize + 1 };
  enum { WarpIndexMask  = Impl::CudaTraits::WarpIndexMask };
  enum { WarpIndexShift = Impl::CudaTraits::WarpIndexShift };

  enum { SharedMemoryBanks = Impl::CudaTraits::SharedMemoryBanks_13 };

  enum { ValueWordCount = ( sizeof(value_type) + sizeof(size_type) - 1 )
                          / sizeof(size_type) };

  /** \brief  If the reduction value occupies an
   *          exact multiple of shared memory banks
   *          then it must be padded to avoid bank conflicts.
   */
  enum { ValueWordStride = ValueWordCount +
          ( ValueWordCount % SharedMemoryBanks ? 0 : 2 ) };

  enum { WordsPerWarp       = ValueWordStride * WarpSize };
  enum { WordsPerWarpStride = ValueWordStride * WarpStride };

private:

  typedef Impl::CudaMultiFunctorParallelReduceMember< void , ReduceTraits , FinalizeType > MemberType ;
  typedef std::vector< MemberType * > MemberVector ;

  MemberVector m_member_functors ;
  FinalizeType m_finalize ;
  size_type    m_shmem_size ;
  size_type    m_stream_count ; ///< Total number of streams
  size_type    m_warps_per_block ;
  size_type    m_threads_per_block ;
  size_type    m_blocks_per_stream ;

public:

  MultiFunctorParallelReduce( const FinalizeType & finalize )
    : m_member_functors()
    , m_finalize( finalize )
    , m_shmem_size( 0 )
    , m_warps_per_block( Impl::cuda_internal_maximum_warp_count() )
    , m_stream_count( Impl::cuda_internal_stream_count() )
    , m_threads_per_block( 0 )
    , m_blocks_per_stream( 0 )
  {
    typedef MultiFunctorParallelReduce< ReduceTraits , FinalizeType , Cuda > self_type ;

    // Consistent block sizes and shared memory so that any block
    // can perform the final reduction.

    const size_type maximum_shared_words = Impl::cuda_internal_maximum_shared_words();

    while ( maximum_shared_words < m_warps_per_block * WordsPerWarpStride ) {
      m_warps_per_block >>= 1 ;
    }

    m_shmem_size = sizeof(size_type) * ( ValueWordStride * (WarpStride * m_warps_per_block - 1) + 1 );

    // Each stream has a range of block offsets

    m_threads_per_block = WarpSize * m_warps_per_block ;
    m_blocks_per_stream = m_threads_per_block / m_stream_count ;
  }

  ~MultiFunctorParallelReduce()
  {
    while ( ! m_member_functors.empty() ) {
      delete m_member_functors.back();
      m_member_functors.pop_back();
    }
  }

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & f )
  {
    const size_type blocks_requested = ( work_count + m_threads_per_block - 1 ) / m_threads_per_block ;

    size_type streams_requested = ( blocks_requested + m_blocks_per_stream - 1 ) / m_blocks_per_stream ;

    if ( m_stream_count < streams_requested ) {
      streams_requested = m_stream_count ;
    }

    // Functor will be executed on stream_count_requested streams

    MemberType * member =
      new Impl::CudaMultiFunctorParallelReduceMember< FunctorType , ReduceTraits , FinalizeType >
        ( f , work_count , streams_requested );

    m_member_functors.push_back( member );
  }

  void execute()
  {
    // Dispatch to streams, second and subsequent dispatches set the 
    // recycling flag so previous reduction results will be read.

    // Each dispatch has a fixed block count, dispatch the functor
    // on up to all streams until the requested block count is met.
    // When the stream_count * blocks_per_stream < blocks_requested
    // then the functor will iterate within the dispatch.

    typename MemberVector::iterator m ;

    const size_type stream_block_count = m_blocks_per_stream * m_stream_count ;

    size_type global_block_count = 0 ;

    for ( m = m_member_functors.begin() ; m != m_member_functors.end() ; ++m ) {
      global_block_count += (*m)->m_stream_count ;
    }

    global_block_count *= m_blocks_per_stream ;

    size_type k = 0 ;
    for ( m = m_member_functors.begin() ; m != m_member_functors.end() ; ++m ) {

      MemberType & member = **m ;

      for ( size_type j = 0 ; j < member.m_stream_count ; ++j , ++k ) {

        // This functor may span multiple streams,
        // provide a block offset for the multiple streams
        const size_type work_block_offset    = m_blocks_per_stream * j ;

        const size_type stream_offset = k % m_stream_count ;
        const size_type stream_block_offset =
          m_blocks_per_stream * stream_offset ;

        const size_type stream_block_recycle = m_stream_count <= k ;

        cudaStream_t & s = Impl::cuda_internal_stream( stream_offset );

        member.execute( m_finalize ,
                        m_blocks_per_stream , work_block_offset ,
                        m_warps_per_block , m_shmem_size , s ,
                        global_block_count ,
                        stream_block_count ,
                        stream_block_offset , stream_block_recycle );
      }
    }
  }
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#endif /* defined( __CUDACC__ ) */

#endif /* KOKKOS_CUDA_PARALLELREDUCE_HPP */

