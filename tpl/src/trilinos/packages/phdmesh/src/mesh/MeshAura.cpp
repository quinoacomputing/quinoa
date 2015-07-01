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

#include <util/ParallelComm.hpp>

#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/Comm.hpp>
#include <mesh/EntityComm.hpp>


namespace phdmesh {

//----------------------------------------------------------------------
// Regenerate the parallel aura of mesh entities that are
// relationed to the shared node mesh entities.

namespace {

bool not_member( const std::vector<EntityProc> & aura_domain ,
                 const EntityProc & ep )
{
  const std::vector<EntityProc>::const_iterator
    i = lower_bound( aura_domain , ep );

  return i == aura_domain.end() || ep != *i ;
}

//----------------------------------------------------------------------

class RegenAuraComm : public EntityComm {
private:
  RegenAuraComm( const RegenAuraComm & );
  RegenAuraComm & operator = ( const RegenAuraComm & );
public:
  RegenAuraComm() {}
  ~RegenAuraComm() {}

  const char * name() const ;

  void receive_entity(
    CommBuffer & buffer ,
    BulkData & receive_mesh ,
    const unsigned send_source ,
    std::vector<EntityProc> & receive_info ) const ;
};

const char * RegenAuraComm::name() const
{ static const char n[] = "phdmesh::RegenAuraComm" ; return n ; }

void RegenAuraComm::receive_entity(
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

  const MetaData & S = receive_mesh.mesh_meta_data();
  Part * const owns_part = & S.locally_owned_part();
  Part * const uses_part = & S.locally_used_part();

  // This is an aura, not a member of the owns part or uses part

  std::vector<Part*>::iterator ip = add.end();

  while ( ip != add.begin() ) {
    --ip ;
    if ( *ip == owns_part || *ip == uses_part ) { ip = add.erase( ip ); }
  }

  EntityProc ep ;

  ep.first = receive_mesh.get_entity( key );

  if ( NULL != ep.first ||
       send_source != owner_rank ||
       send_source == receive_mesh.parallel_rank() ) {
    std::string msg( "phdmesh::RegenAuraComm::receive_entity FAILED" );
    throw std::logic_error( msg );
  }

  ep.first = & receive_mesh.declare_entity( key , add , owner_rank );

  receive_mesh.declare_relation( *ep.first , relations );

  ep.second = owner_rank ;

  receive_info.push_back( ep );
}

//----------------------------------------------------------------------

void scrub( BulkData & M , std::vector<EntityProc> & new_domain ,
                       std::vector<EntityProc> & new_range )
{
  static const char method[] = "phdmesh::comm_mesh_regenerate_aura::scrub" ;

  Part & uses_part = M.mesh_meta_data().locally_used_part();
  const unsigned p_rank = M.parallel_rank();
  const unsigned p_size = M.parallel_size();

  std::vector<EntityProc>::iterator i ;

  for ( i = new_range.end() ; i != new_range.begin() ; ) {
    --i ;
    Entity & aura_entity = *i->first ;

    const unsigned aura_type = aura_entity.entity_type();

    if ( aura_entity.kernel().has_superset( uses_part ) ) {
      i->first = NULL ;
    }
    else {

      bool destroy_it = true ;

      for ( PairIterRelation aura_con = aura_entity.relations();
            destroy_it && aura_con ; ++aura_con ) {
        Entity & e = * aura_con->entity();
        const unsigned e_type = aura_con->entity_type();

        if ( e_type < aura_type ) { // Aura Uses this entity
          if ( e.sharing() ) {
            destroy_it = false ; // Uses a shared entity
          }
        }
        else if ( aura_type < e_type ) { // Aura UsedBy this entity
          if ( ! e.kernel().has_superset( uses_part ) ) {
            destroy_it = false ; // Used by a retained aura entity
          }
          else {
            std::ostringstream msg ;
            msg << "P" << p_rank ;
            msg << " " << method ;
            msg << " FAILED Aura " ;
            print_entity_key( msg , aura_entity.key() );
            msg << " UsedBy non-aura " ;
            print_entity_key( msg , e.key() );
            throw std::runtime_error( msg.str() );
          }
        }
      }

      if ( destroy_it ) {
        M.destroy_entity( i->first );
        i->first = NULL ;
      }
    }
  }

  // Inform the owners of the destroyed aura entities

  CommAll all( M.parallel() );

  bool change = false ;

  for ( i = new_range.begin() ; i != new_range.end() ; ++i ) {
    const unsigned char flag = NULL == i->first ;
    if ( flag ) { change = true ; }
    all.send_buffer( i->second ).skip<unsigned char>(1);
  }

  change = all.allocate_buffers( p_size / 4 , false , change );

  if ( change ) {

    for ( i = new_range.begin() ; i != new_range.end() ; ++i ) {
      const unsigned char flag = NULL == i->first ;
      all.send_buffer( i->second ).pack<unsigned char>(flag);
    }

    all.communicate();

    for ( i = new_domain.begin() ; i != new_domain.end() ; ++i ) {
      unsigned char flag ;
      all.recv_buffer( i->second ).unpack<unsigned char>(flag);
      if ( flag ) { i->first = NULL ; }
    }

    for ( i = new_domain.end() ; i != new_domain.begin() ; ) {
      --i ;
      if ( NULL == i->first ) { i = new_domain.erase(i); }
    }

    for ( i = new_range.end() ; i != new_range.begin() ; ) {
      --i ;
      if ( NULL == i->first ) { i = new_range.erase(i); }
    }
  }
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void comm_mesh_regenerate_aura( BulkData & M )
{
  static const char method[] = "phdmesh::comm_mesh_regenerate_aura" ;

  const RegenAuraComm mgr ;
  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  const std::vector<EntityProc> & shares = M.shared_entities();

  std::vector<EntityProc> old_aura_domain( M.ghost_source() );
  std::vector<EntityProc> old_aura_range(  M.ghost_destination() );

  scrub( M , old_aura_domain , old_aura_range );

  // Run the shared entities and push all 'UsedBy' entities
  // into the send vector.

  std::vector<EntityProc> new_aura_domain ;
  std::vector<EntityProc> new_aura_range ;

  for ( std::vector<EntityProc>::const_iterator
        i = shares.begin() ; i != shares.end() ; ) {

    const std::vector<EntityProc>::const_iterator ib = i ;
    for ( ; i != shares.end() && ib->first == i->first ; ++i );
    const std::vector<EntityProc>::const_iterator ie = i ;

    Entity * const shared_entity = ib->first ;
    const EntityType shared_type = shared_entity->entity_type();

    for ( PairIterRelation shared_entity_con = shared_entity->relations() ;
          shared_entity_con ; ++shared_entity_con ) {

      Entity * const e = shared_entity_con->entity();
      const EntityType e_type = shared_entity_con->entity_type();

      if ( shared_type < e_type && p_rank == e->owner_rank() ) {

        // "shared_entity->UsedBy->e" and 'e' is locally owned
        // Send to each shares processor

        for ( std::vector<EntityProc>::const_iterator
              j = ib ; j != ie ; ++j ) {
          EntityProc entry ;

          entry.first = e ;
          entry.second = j->second ;

          if ( not_member( shares , entry ) &&
               not_member( old_aura_domain , entry ) ) {
            new_aura_domain.push_back( entry );
          }

          // and all of the shared_entity->UsedBy->entity->Uses->{entities}

          for ( PairIterRelation jc = e->relations(); jc ; ++jc ) {

            if ( jc->entity_type() < e_type && shared_entity != jc->entity() ) {

              // e->Uses->(jc->entity())

              entry.first = jc->entity();
              if ( not_member( shares , entry ) ) {
                if ( p_rank == entry.first->owner_rank() ) {
                  if ( not_member( old_aura_domain , entry ) ) {
                    new_aura_domain.push_back( entry );
                  }
                }
                else if ( entry.second != entry.first->owner_rank() ) {
                  // Inform the owner of this need
                  new_aura_range.push_back( entry );
                }
              }
            }
          }
        }
      }
    }
  }

  {
    CommAll all( M.parallel() );

    std::vector<EntityProc>::iterator j ;

    for ( j = new_aura_range.begin() ; j != new_aura_range.end() ; ++j ) {
      all.send_buffer( j->first->owner_rank() ).skip<entity_key_type>(2);
    }

    all.allocate_buffers( p_size / 4 , false /* not symmetric */ );

    for ( j = new_aura_range.begin() ; j != new_aura_range.end() ; ++j ) {
      entity_key_type data[2] ;
      data[0] = j->first->key();
      data[1] = j->second ;
      all.send_buffer( j->first->owner_rank() ).pack<entity_key_type>(data,2);
    }

    new_aura_range.clear();

    all.communicate();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = all.recv_buffer(p);
      while ( buf.remaining() ) {
        entity_key_type data[2] ;
        buf.unpack<entity_key_type>( data , 2 );
        EntityProc tmp ;
        tmp.first  = M.get_entity( data[0] , method );
        tmp.second = data[1] ;
        if ( not_member( old_aura_domain , tmp ) ) {
          new_aura_domain.push_back( tmp );
        }
      }
    }
  }

  sort_unique( new_aura_domain );

  communicate_entities( mgr , M , M , new_aura_domain , new_aura_range , false );

  new_aura_domain.insert( new_aura_domain.end() , old_aura_domain.begin() ,
                                                  old_aura_domain.end() );

  new_aura_range.insert( new_aura_range.end() , old_aura_range.begin() ,
                                                old_aura_range.end() );

  sort_unique( new_aura_domain );
  sort_unique( new_aura_range );

  M.set_ghosting( new_aura_domain , new_aura_range );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

