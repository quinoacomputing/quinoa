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
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <sstream>
#include <util/FixedPoolAlloc.hpp>

namespace phdmesh {

void throw_fixed_pool_buffer_bad_size( const std::size_t nbyte_total ,
                                       const std::size_t nbyte_first ,
                                       const std::size_t nbyte )
{
  std::ostringstream msg ;
  msg << "phdmesh::FixedPoolBuffer<" ;
  msg << nbyte_total ;
  msg << ">::{de}allocate( nbyte = " ;
  msg << nbyte ;
  msg << " ) != ( nbyte_first = " ;
  msg << nbyte_first ;
  msg << " ) ) FAILED TO MAINTAIN FIXED ALLOCATION SIZE" ;
  throw std::invalid_argument( msg.str() );
}

void throw_fixed_pool_buffer_exhausted( const std::size_t nbyte_total ,
                                        const std::size_t nbyte )
{
  std::ostringstream msg ;
  msg << "phdmesh::FixedPoolBuffer<" ;
  msg << nbyte_total ;
  msg << ">::allocate(" ;
  msg << nbyte ;
  msg << ") FAILED: MEMORY EXHAUSTED" ;
  throw std::runtime_error( msg.str() );
}

void throw_fixed_pool_buffer_bad_deallocate( const std::size_t nbyte_total ,
                                             void * const p )
{
  std::ostringstream msg ;
  msg << "phdmesh::FixedPoolBuffer<" ;
  msg << nbyte_total ;
  msg << ">::deallocate(" ;
  msg << p ;
  msg << ") FAILED: POINTER NOT IN RANGE" ;
  throw std::runtime_error( msg.str() );
}

void **fixed_pool_buffer_init( const std::size_t nbyte_total ,
                               const std::size_t nbyte ,
                               void ** p )
{
  const std::size_t rinc = nbyte % sizeof(void*);
  const std::size_t ninc = nbyte / sizeof(void*) +
                           ( rinc ? ( sizeof(void*) - rinc ) : 0 );
  const std::size_t rtot = nbyte_total % sizeof(void*);
  const std::size_t ntot = nbyte_total / sizeof(void*) +                            ( rtot ? ( sizeof(void*) - rtot ) : 0 );
  const std::size_t nelem = ntot / ninc ;

  if ( ! nelem ) {
    std::ostringstream msg ;
    msg << "phdmesh::FixedPoolBuffer<" ;
    msg << nbyte_total ;
    msg << ">::allocate(" ;
    msg << nbyte ;
    msg << ") FAILED: CANNOT HOLD EVEN ONE MEMBER" ;
    throw std::invalid_argument( msg.str() );
  }

  void ** const ep = p + ninc * ( nelem - 1 );
  while ( ep != p ) { void ** const n = p + ninc ; *p = n ; p = n ; }
  *p = NULL ;
  return p ;
}


}

