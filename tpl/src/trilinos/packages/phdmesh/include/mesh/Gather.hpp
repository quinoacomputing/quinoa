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

#ifndef phdmesh_Gather_hpp
#define phdmesh_Gather_hpp

#include <mesh/FieldData.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

template< unsigned NType , enum EntityType EType ,
          unsigned NRel , class field_type >
bool gather_field_data( const field_type & field ,
                        const Entity     & entity ,
                        typename FieldTraits< field_type >::data_type * dst )
{
  typedef typename FieldTraits< field_type >::data_type T ;

  PairIterRelation rel = entity.relations( EType );

  bool result = NRel == (unsigned) rel.size();

  if ( result ) {
    T * const dst_end = dst + NType * NRel ;
    for ( const T * src ;
          ( dst < dst_end ) &&
          ( src = field_data( field , * rel->entity() ) ) ;
          ++rel , dst += NType ) {
      Copy<NType>( dst , src );
    }
    result = dst == dst_end ;
  }
  return result ;
}

//----------------------------------------------------------------------

}

#endif

