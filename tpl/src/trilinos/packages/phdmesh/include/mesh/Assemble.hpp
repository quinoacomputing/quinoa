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

#ifndef phdmesh_Assemble_hpp
#define phdmesh_Assemble_hpp

#include <cstddef>
#include <vector>

#include <util/Basics.hpp>

#include <mesh/Types.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/Comm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

typedef void (*assemble_function_ptr)( void * dst ,
                                       const void * src ,
                                       const unsigned id );

namespace {

template< template<unsigned,unsigned> class Op , unsigned N , typename T>
void assemble_function( void * v_dst , const void * v_src , const unsigned id )
{
  Op<N,0>( reinterpret_cast<      T*>(v_dst) ,
           reinterpret_cast<const T*>(v_src) + id * N );
}

}

class Assemble {
private:

  assemble_function_ptr op ;
  const FieldBase * dst ;
  const FieldBase * src ;
  unsigned          src_rank ;

public:

  EntityType src_type() const { return EntityType( src_rank ); }
  const FieldBase & dst_field() const { return *dst ; }
  const FieldBase & src_field() const { return *src ; }

  void operator()( void * v_dst, const void * v_src, const unsigned id ) const
    { (*op)( v_dst , v_src , id ); }

  ~Assemble() {}

  Assemble() : op( NULL ), dst( NULL ), src( NULL ) {}

  Assemble( const Assemble & rhs )
    : op( rhs.op ), dst( rhs.dst ), src( rhs.src ) {}

  Assemble & operator = ( const Assemble & rhs )
    { op = rhs.op ; dst = rhs.dst ; src = rhs.src ; return *this ; }

  Assemble( assemble_function_ptr arg_op ,
            const FieldBase & arg_dst ,
            const FieldBase & arg_src )
    : op( arg_op ), dst( & arg_dst ), src( & arg_src ) {}

  template< template<unsigned,unsigned> class Op , unsigned N , typename T >
  Assemble( const Op<N,0> & ,
            const FieldBase & arg_dst ,
            const FieldBase & arg_src )
    : op( & assemble_function< Op<N,0> , T > ),
      dst( & arg_dst ), src( & arg_src ) {}
};

/** Assemble partial sum contributions from 'src_field' to 'dst_field'.
 *  The dimensions of these fields must be compatible such that
 *    Op<NT>( dst[N].data( dst_field ) , src.data( src_field )[N]
 *    src.data( src_field )(NT,N) -> dst[N].data( dst_field )[NT]
 *  Prior to performing assembly copy the source field values
 *  from the owning processor to the aura processor.
 *  Then assemble contributions ordered by the source entities' identifiers.
 *  This guarantees a consistent result no matter what the
 *  domain decomposition.
 */

/** Order and parallel consistent assembly of fields.
 *  Updates the contributed field values from aura entities.
 *  Usage:
 *    std::vector<Assemble> a(2) ;
 *    a[0] = Assemble( Sum<3>() , dst_field_x , src_field_x );
 *    a[1] = Assemble( Max<1>() , dst_field_y , src_field_y );
 *    assemble( a );
 */
void assemble( const std::vector<Assemble> & );

//----------------------------------------------------------------------

}

#endif

