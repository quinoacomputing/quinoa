/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   August 2007
 */

#ifndef util_FixedPoolAlloc_hpp
#define util_FixedPoolAlloc_hpp

#include <cstddef>

namespace phdmesh {

void ** fixed_pool_buffer_init( const std::size_t nbyte_total ,
                                const std::size_t nbyte ,
                                void ** );

void throw_fixed_pool_buffer_bad_size( const std::size_t nbyte_total ,
                                       const std::size_t nbyte_first ,
                                       const std::size_t nbyte );

void throw_fixed_pool_buffer_exhausted( const std::size_t nbyte_total ,
                                        const std::size_t nbyte );

void throw_fixed_pool_buffer_bad_deallocate( const std::size_t ,
                                             void * const );

template<unsigned NBYTE>
class FixedPoolBuffer {
private:

  enum { NWORD = NBYTE / sizeof(void*) };

  std::size_t size ;
  void     ** link ;
  void      * mem[ NWORD ];

  FixedPoolBuffer( const FixedPoolBuffer<NBYTE> & );
  FixedPoolBuffer<NBYTE> & operator = ( const FixedPoolBuffer<NBYTE> & );

public:

  ~FixedPoolBuffer() {}

  FixedPoolBuffer() : size(0) { link = mem ; }

  bool full() const { return ! link ; }

  void * allocate( std::size_t nbyte )
    {
      void * p = NULL ;
      if ( nbyte ) {
        if ( ! size ) { // Initialize on the first call
          size = nbyte ;
          fixed_pool_buffer_init(NBYTE,size,link);
        }
        else if ( nbyte != size ) {
          throw_fixed_pool_buffer_bad_size(NBYTE,size,nbyte);
        }
        if ( ! link ) {
          throw_fixed_pool_buffer_exhausted(NBYTE,nbyte);
        }
        p = link ;
        link = reinterpret_cast<void**>( *link );
      }
      return p ;
    }

  void deallocate( void * p , std::size_t nbyte )
    {
      if ( nbyte != size ) {
        throw_fixed_pool_buffer_bad_size(NBYTE,size,nbyte);
      }
      void ** const vp = reinterpret_cast<void**>(p);
      if ( vp < mem || ( mem + NWORD - size ) < vp ) {
        throw_fixed_pool_buffer_bad_deallocate(NBYTE,p);
      }
      *vp = link ;
      link = reinterpret_cast<void**>(vp);
      return ;
    }
};


template<unsigned NBYTE,typename T=void*>
class FixedPoolAllocator {
private:

  FixedPoolBuffer<NBYTE> * buffer ;

  FixedPoolAllocator<NBYTE,T> &
    operator = ( const FixedPoolAllocator<NBYTE,T> & );

  template<unsigned M,typename U> friend class FixedPoolAllocator ;
 
public:

  typedef FixedPoolBuffer<NBYTE> buffer_type ;

  typedef           T value_type;
  typedef std::size_t size_type;
  typedef          T* pointer;
  typedef    const T* const_pointer;
  typedef          T& reference;
  typedef    const T& const_reference;
  typedef std::ptrdiff_t difference_type;

  pointer       address( reference       value ) const { return & value ; }
  const_pointer address( const_reference value ) const { return & value ; }

  ~FixedPoolAllocator() {}

  FixedPoolAllocator() : buffer(NULL) {}

  FixedPoolAllocator( const FixedPoolAllocator<NBYTE,T> & a )
    : buffer(a.buffer) {}

  explicit FixedPoolAllocator( buffer_type & b ) : buffer( & b ) {}

  template<typename U>
  struct rebind { typedef FixedPoolAllocator<NBYTE,U> other ; };

  template<typename U>
  FixedPoolAllocator( const FixedPoolAllocator<NBYTE,U> & a )
    : buffer(a.buffer) {}

  pointer allocate( size_type n , const void * = NULL )
    { return reinterpret_cast<pointer>( buffer->allocate( n * sizeof(T) ) ); }

  void deallocate( pointer p , size_type n )
    { buffer->deallocate(p, n * sizeof(T) ); }

  void construct( pointer p , const T & val ) { new(p) T(val); }
  void destroy(   pointer p ) { p->~T(); }

  template<typename U>
  void construct( U * p , const U & val ) { new(p) U(val); }

  template<typename U>
  void destroy( U * p ) { p->~U(); }

  size_type max_size() const { return NBYTE / sizeof(T) ; }
};

template<unsigned NBYTE,typename T>
bool operator == ( const FixedPoolAllocator<NBYTE,T> & lhs ,
                   const FixedPoolAllocator<NBYTE,T> & rhs )
{ return lhs.buffer == rhs.buffer ; }

template<unsigned NBYTE,typename T>
bool operator != ( const FixedPoolAllocator<NBYTE,T> & lhs ,
                   const FixedPoolAllocator<NBYTE,T> & rhs )
{ return lhs.buffer != rhs.buffer ; }

}

#endif

