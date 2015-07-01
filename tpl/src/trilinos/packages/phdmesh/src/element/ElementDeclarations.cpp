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

#include <mesh/MetaData.hpp>
#include <mesh/Part.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/Entity.hpp>
#include <element/CellTopology.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

const CellTopology * get_cell_topology( const Part & p )
{ return p.attribute<CellTopology>(); }

void set_cell_topology( Part & p , const CellTopology * singleton )
{
  static const char method[] = "phdmesh::set_cell_topology" ;

  MetaData & m = p.mesh_meta_data();

  const CellTopology * t = NULL ;

  if ( singleton == NULL ||
       singleton != ( t = m.declare_attribute_no_delete(p,singleton) ) ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( " ;
    msg << p.name();
    msg << " , " ;
    if ( singleton ) { msg << singleton->name ; }
    else             { msg << "NULL" ; }
    msg << " ) ERROR" ;
    if ( t ) { msg << "Existing topology = " << t->name ; }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

const CellTopology * get_cell_topology( const Kernel & kernel )
{
  const CellTopology * top = NULL ;
  PartSet parts ;
  kernel.supersets( parts );

  PartSet::iterator i = parts.begin() ;

  for ( ; NULL == top && i != parts.end() ; ++i ) {
    top = get_cell_topology( **i );
  }

  bool ok = true ;

  for ( ; ok && i != parts.end() ; ++i ) {
    const CellTopology * const tmp = get_cell_topology( **i );
    ok = tmp == NULL || tmp == top ;
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "phdmesh::get_cell_topology( Kernel[" ;
    for ( i = parts.begin() ; i != parts.end() ; ++i ) {
      const CellTopology * const tmp = get_cell_topology( **i );
      msg << " " << (*i)->name();
      if ( top ) { msg << "->" << tmp->name ; }
      msg << " ] ) FAILED WITH MULTIPLE LOCAL TOPOLOGIES" ;
      throw std::runtime_error( msg.str() );
    }
  }

  return top ;
}

const CellTopology * get_cell_topology( const Entity & entity )
{ return get_cell_topology( entity.kernel() ); }

//----------------------------------------------------------------------

Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const unsigned elem_id ,
                          Entity * node[] )
{
  static const char method[] = "phdmesh::declare_element" ;

  const CellTopology * const top = get_cell_topology( part );

  if ( top == NULL ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( mesh , " ;
    msg << part.name();
    msg << " , " ;
    msg << elem_id ;
    msg << " , node[] ) ERROR, Part does not have a local topology" ;
    throw std::runtime_error( msg.str() );
  }

  const EntityType type = Element ;

  PartSet add ;

  { Part * const tmp = & part ; add.push_back( tmp ); }

  const entity_key_type key = entity_key( type , elem_id );

  Entity & elem = mesh.declare_entity( key , add );

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    mesh.declare_relation( elem , * node[i] , i );
  }
  return elem ;
}

//----------------------------------------------------------------------

Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const unsigned elem_id ,
                          const unsigned node_id[] )
{
  static const char method[] = "phdmesh::declare_element" ;

  const CellTopology * const top = get_cell_topology( part );

  if ( top == NULL ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( mesh , " ;
    msg << part.name();
    msg << " , " ;
    msg << elem_id ;
    msg << " , node_id[] ) ERROR, Part does not have a local topology" ;
    throw std::runtime_error( msg.str() );
  }

  const EntityType type = Element ;

  PartSet add ;

  { Part * const tmp = & part ; add.push_back( tmp ); }

  const entity_key_type key = entity_key( type , elem_id );

  Entity & elem = mesh.declare_entity( key , add );

  for ( unsigned i = 0 ; i < top->node_count ; ++i ) {
    const entity_key_type node_key = entity_key( Node , node_id[i] );

    Entity & node = mesh.declare_entity( node_key , add );

    mesh.declare_relation( elem , node , i );
  }
  return elem ;
}

//----------------------------------------------------------------------

Entity & declare_element_side(
  BulkData & mesh ,
  const unsigned global_side_id ,
  Entity & elem ,
  const unsigned local_side_id )
{
  static const char method[] = "phdmesh::declare_element_side" ;

  const CellTopology * const elem_top = get_cell_topology( elem );

  const CellTopology * const side_top =
    ( elem_top && local_side_id < elem_top->side_count )
    ? elem_top->side[ local_side_id ].topology : NULL ;

  if ( NULL == side_top ) {
     std::ostringstream msg ;
     msg << method << "( mesh , "
         << global_side_id
         << " , " ;
     print_entity_key( msg , elem.key() );
     msg << " , "
         << local_side_id
         << " ) FAILED" ;
     if ( NULL == elem_top ) {
       msg << " Cannot discern element topology" ;
     }
     else {
       msg << " Cell side id exceeds " ;
       msg << elem_top->name ;
       msg << ".side_count = " ;
       msg << elem_top->side_count ;
     }
     throw std::runtime_error( msg.str() );
   }

  const unsigned * const side_node_map = elem_top->side[ local_side_id ].node ;

  const EntityType      side_type = (EntityType) side_top->dimension ;
  const entity_key_type side_key  = entity_key( side_type, global_side_id );

  PartSet parts ;

  elem.kernel().supersets( parts );

  Entity & side = mesh.declare_entity( side_key , parts , elem.owner_rank() );

  PairIterRelation rel = elem.relations( Node );

  for ( unsigned i = 0 ; i < side_top->node_count ; ++i ) {
    Entity & node = * rel[ side_node_map[i] ].entity();
    mesh.declare_relation( side , node , i );
  }

  mesh.declare_relation( elem , side , local_side_id );

  return side ;
}

}

