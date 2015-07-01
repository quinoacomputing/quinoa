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
#include <iostream>
#include <sstream>
#include <algorithm>

#include <mesh/Entity.hpp>
#include <mesh/BulkData.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

std::ostream &
print_entity_key( std::ostream & os , EntityType type , entity_id_type id )
{
  const char * const name = entity_type_name( type );
  return os << name << "[" << id << "]" ;
}

std::ostream &
print_entity_key( std::ostream & os , entity_key_type key )
{
  const EntityType type   = entity_type(key);
  const entity_id_type id = entity_id(key);
  return print_entity_key( os , type , id );
}

std::ostream &
print_entity( std::ostream & os , const std::string & lead , const Entity & e )
{
  print_entity_key( os , e.key() );
  os << " Owner(P" << e.owner_rank() << ") Relationships {" ;

  for ( PairIterRelation con = e.relations() ; con ; ++con ) {
    os << std::endl << lead << "  " << *con ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

namespace {
const entity_key_type & zero_key()
{ static entity_key_type z = 0 ; return z ; }
}

Entity::Entity()
  : SetvMember< entity_key_type >( zero_key() ),
    m_relation(), m_kernel(), m_kernel_ord(0), m_owner_rank(0)
{}

Entity::~Entity()
{
  if ( m_kernel ) {
    std::ostringstream msg ;
    msg << "phdmesh::Entity::~Entity() ERROR: " ;
    print_entity_key( msg , key() );
    msg << " is still a member of a phdmesh::Kernel" ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

