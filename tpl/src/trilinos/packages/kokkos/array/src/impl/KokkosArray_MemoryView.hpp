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

#ifndef KOKKOS_MEMORYVIEW_HPP
#define KOKKOS_MEMORYVIEW_HPP

#include <impl/KokkosArray_forward.hpp>

namespace KokkosArray {
namespace Impl {

template< typename ValueType , class DeviceType >
class MemoryView {
private:

  MemoryView( const MemoryView & );
  MemoryView & operator = ( const MemoryView & );

public:

  /** \brief  Only defined on the device */
  ValueType * ptr_on_device() const ;

  MemoryView();
  ~MemoryView();

  /** \brief  If a non-null view */
  operator bool() const ;

  /** \brief  View to the same memory */
  bool operator == ( const MemoryView & ) const ;

  /** \brief  Not view to the same memory */
  bool operator != ( const MemoryView & ) const ;
};

template< class DeviceType > class MemoryManager ;

//----------------------------------------------------------------------------

template< class ValueType , class Device >
struct Factory< MemoryView< ValueType , Device > , Impl::MirrorUseView >
{
  typedef MemoryView< ValueType , Device > output_type ;

  static inline
  const output_type & create( const output_type & input , const size_t )
  { return input ; }

  template< class DeviceInput >
  static inline
  output_type create( const MemoryView< ValueType , DeviceInput > & input ,
                      const size_t count )
  {
    output_type output ;
    output.allocate( count , std::string() );
    return output ;
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* KOKKOS_MEMORYVIEW_HPP */

