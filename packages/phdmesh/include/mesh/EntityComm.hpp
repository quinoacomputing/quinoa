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

#ifndef phdmesh_EntityComm_hpp
#define phdmesh_EntityComm_hpp

//----------------------------------------------------------------------

#include <vector>

#include <util/ParallelComm.hpp>
#include <mesh/Types.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

//----------------------------------------------------------------------
/** Send local mesh entities to 'recv_mesh' according to 'send'.
 *  Received mesh entities are identified in 'recv'.
 */

class EntityComm {
public:
  EntityComm() {}
  virtual ~EntityComm();
  virtual const char * name() const ;

  /** Send entity to receive_mesh and processors in the range [ibeg,iend).
   *  Default behavior is to call pack_entity.
   */
  virtual void send_entity(
    CommBuffer & buffer ,
    const BulkData & receive_mesh ,
    const std::vector<EntityProc>::const_iterator ibeg ,
    const std::vector<EntityProc>::const_iterator iend ) const ;

  /** Receive entity from send_source and update receive_info accordingly.
   *  Default behavior is to invoke unpack_entity, declare the entity
   *  via declare_entity, and then introduce the (entity,send_source)
   *  pair into receieve_info.
   */
  virtual void receive_entity(
    CommBuffer & buffer ,
    BulkData & receive_mesh ,
    const unsigned send_source ,
    std::vector<EntityProc> & receive_info ) const ;

  /** Pack entity information for declaration in the receive mesh.
   *  Parts are mapped from the send mesh to the receive mesh.
   *  Include the span of destination processors.
   */
  void pack_entity( CommBuffer & , const BulkData & ,
                    const std::vector<EntityProc>::const_iterator ,
                    const std::vector<EntityProc>::const_iterator ) const ;

  /** Unpack entity information filled by pack_entity. */
  void unpack_entity( CommBuffer & , const BulkData & ,
                      entity_key_type & entity_key ,
                      unsigned        & owner_rank ,
                      std::vector<Part*> & parts ,
                      std::vector<Relation> & relations ,
                      std::vector<unsigned> & send_destinations ) const ;

  /** Pack an entity's field values into a buffer */
  void pack_field_values( CommBuffer & , Entity & ) const ;

  /** Unpack an entity's field values from a buffer */
  void unpack_field_values( CommBuffer & , Entity & ) const ;

private:
  EntityComm( const EntityComm & );
  EntityComm & operator = ( const EntityComm & );
};


} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

