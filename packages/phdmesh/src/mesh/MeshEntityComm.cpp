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

#include <iterator>
#include <stdexcept>
#include <sstream>

#include <mesh/BulkData.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/EntityComm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

EntityComm::~EntityComm() {}

const char * EntityComm::name() const
{
  static const char my_name[] = "phdmesh::EntityComm" ;
  return my_name ;
}

//----------------------------------------------------------------------

void EntityComm::send_entity(
  CommBuffer & buffer ,
  const BulkData & receive_mesh ,
  const std::vector<EntityProc>::const_iterator ibeg ,
  const std::vector<EntityProc>::const_iterator iend ) const
{
  pack_entity( buffer , receive_mesh , ibeg , iend );
}

void EntityComm::receive_entity(
  CommBuffer              & buffer ,
  BulkData                    & receive_mesh ,
  const unsigned            send_source ,
  std::vector<EntityProc> & receive_info ) const
{
  entity_key_type       key ;
  unsigned              owner_rank ;
  std::vector<Part*>    add ;
  std::vector<Relation> relations ;
  std::vector<unsigned> send_destinations ;

  unpack_entity( buffer , receive_mesh ,
                 key , owner_rank ,
                 add , relations , send_destinations );

  EntityProc ep ;

  ep.first = & receive_mesh.declare_entity( key , add , owner_rank );

  receive_mesh.declare_relation( *ep.first , relations );

  ep.second = send_source ;

  receive_info.push_back( ep );
}

//----------------------------------------------------------------------

void EntityComm::pack_entity(
  CommBuffer & buf ,
  const BulkData & recv_mesh ,
  const std::vector<EntityProc>::const_iterator ibeg ,
  const std::vector<EntityProc>::const_iterator iend ) const
{
  const MetaData      & recv_mesh_meta_data = recv_mesh.mesh_meta_data();
  const Entity      & entity      = * ibeg->first ;
  const entity_key_type key       = entity.key();
  const unsigned      owner_rank  = entity.owner_rank();
  const Kernel      & kernel      = entity.kernel();
  const MetaData      & send_mesh_meta_data = kernel.mesh().mesh_meta_data();
  const bool          same_mesh_meta_data = & send_mesh_meta_data == & recv_mesh_meta_data ;
  const unsigned      dest_size   = std::distance( ibeg , iend );
        PairIterRelation  rel         = entity.relations();

  const std::pair<const unsigned * , const unsigned * >
    tmp = entity.kernel().superset_part_ordinals();

  std::vector<unsigned> part_ordinals( tmp.first , tmp.second );

  if ( ! same_mesh_meta_data ) { // Map the parts by name
    std::vector<unsigned>::iterator i = part_ordinals.begin();

    while ( i != part_ordinals.end() ) {
      Part &       send_p = send_mesh_meta_data.get_part( *i );
      Part * const recv_p = recv_mesh_meta_data.get_part( send_p.name() );
      if ( recv_p != NULL ) {
        *i = recv_p->mesh_meta_data_ordinal(); ++i ;
      }
      else {
        i = part_ordinals.erase(i);
      }
    }
  }

  buf.pack<entity_key_type>( key );
  buf.pack<unsigned>( owner_rank );

  // Parts:
  {
    const unsigned n = part_ordinals.size();
    buf.pack<unsigned>( n );
    buf.pack<unsigned>( & part_ordinals[0] , n );
  }

  // Relationships:
  {
    const unsigned rel_size = rel.size();
    buf.pack<unsigned>( rel_size );
    uint64_type rel_data[2] ;
    for ( ; rel ; ++rel ) {
      rel_data[0] = rel->entity()->key();
      rel_data[1] = rel->attribute();
      buf.pack<uint64_type>( rel_data , 2 );
    }
  }

  // Processors:
  buf.pack<unsigned>( dest_size );
  for ( std::vector<EntityProc>::const_iterator i = ibeg ; i != iend ; ++i ) {
    buf.pack<unsigned>( i->second );
  }
}

void EntityComm::unpack_entity(
  CommBuffer            & recv_buf ,
  const BulkData            & recv_mesh ,
  entity_key_type       & key ,
  unsigned              & owner_rank ,
  std::vector<Part*>    & parts ,
  std::vector<Relation> & relations ,
  std::vector<unsigned> & send_destinations ) const
{
  const MetaData & recv_mesh_meta_data = recv_mesh.mesh_meta_data();

  parts.clear();
  relations.clear();
  send_destinations.clear();

  recv_buf.unpack<entity_key_type>( key );

  recv_buf.unpack<unsigned>( owner_rank );

  // Parts:
  {
    unsigned n ; recv_buf.unpack<unsigned>( n );
    for ( ; n ; --n ) {
      unsigned ordinal ; recv_buf.unpack<unsigned>( ordinal );
      Part * const p = & recv_mesh_meta_data.get_part( ordinal );
      parts.push_back( p );
    }
  }

  // Relationships:
  {
    unsigned rel_size ; recv_buf.unpack<unsigned>( rel_size );

    relations.reserve( rel_size );

    uint64_type rel_data[2] ;

    for ( unsigned j = 0 ; j < rel_size ; ++j ) {
      recv_buf.unpack<uint64_type>( rel_data , 2 );

      Entity * rel_entity = recv_mesh.get_entity( rel_data[0] );

      if ( rel_entity != NULL ) {
        Relation rel( (relation_attr_type) rel_data[1] , *rel_entity );
        relations.push_back( rel );
      }
    }
  }

  // Processors:
  {
    unsigned dest_size ; recv_buf.unpack<unsigned>( dest_size );

    const unsigned u_zero = 0 ;
    send_destinations.assign( dest_size , u_zero );

    unsigned * const tmp = & send_destinations[0] ;
    recv_buf.unpack<unsigned>( tmp , dest_size );
  }
}

//----------------------------------------------------------------------

void EntityComm::pack_field_values(
  CommBuffer & buf , Entity & entity ) const
{
  const Kernel & kernel = entity.kernel();
  const BulkData   & mesh   = kernel.mesh();
  const MetaData & mesh_meta_data = mesh.mesh_meta_data();

  const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields();

  const std::vector< FieldBase * >::const_iterator i_end = fields.end();
  const std::vector< FieldBase * >::const_iterator i_beg = fields.begin();

  std::vector< FieldBase * >::const_iterator i ;

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    buf.pack<unsigned>( size );
  }

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    if ( size ) {
      unsigned char * const ptr = (unsigned char *) field_data( f , entity );
      buf.pack<unsigned char>( ptr , size );
    }
  }
}

void EntityComm::unpack_field_values(
  CommBuffer & buf , Entity & entity ) const
{
  static const char method[] = "phdmesh::EntityComm::unpack_field_values" ;

  const Kernel & kernel = entity.kernel();
  const BulkData   & mesh   = kernel.mesh();
  const MetaData & mesh_meta_data = mesh.mesh_meta_data();

  const std::vector< FieldBase * > & fields = mesh_meta_data.get_fields();

  const std::vector< FieldBase * >::const_iterator i_end = fields.end();
  const std::vector< FieldBase * >::const_iterator i_beg = fields.begin();

  std::vector< FieldBase * >::const_iterator i ;

  std::ostringstream msg ;
  bool ok = true ;

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    unsigned recv_data_size ; buf.unpack<unsigned>( recv_data_size );
    if ( size != recv_data_size ) {
      if ( ok ) {
        msg << "P" << mesh.parallel_rank();
        msg << ": " << method ;
        msg << "( " ;
        print_entity_key( msg , entity.key() );
        msg << " ) FAILED, incompatible size for field {" ;
      }
      msg << " " << (*i)->name();
      ok = false ;
    }
  }

  if ( ! ok ) {
    msg << " }" ;
    throw std::runtime_error( msg.str() );
  }

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    if ( size ) {
      unsigned char * ptr = (unsigned char *) field_data( f , entity );
      buf.unpack<unsigned char>( ptr , size );
    }
  }
}

//----------------------------------------------------------------------

}

