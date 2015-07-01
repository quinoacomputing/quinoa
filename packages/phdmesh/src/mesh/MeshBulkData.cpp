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
#include <util/ParallelReduce.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/Comm.hpp>
#include <mesh/FieldData.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

BulkData::BulkData( const MetaData & mesh_meta_data ,
            ParallelMachine parallel ,
             unsigned kernel_capacity )
  : m_mesh_meta_data( mesh_meta_data ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_kernel_capacity( kernel_capacity ),
    m_sync_change_count( 0 )
{
  static const char method[] = "phdmesh::BulkData::Mesh" ;

  m_mesh_meta_data.assert_committed( method );

  verify_parallel_consistency( mesh_meta_data , parallel );
}

BulkData::~BulkData()
{
  m_shares_all.clear();
  m_aura_domain.clear();
  m_aura_range.clear();

  // Remove entities from the kernels.
  // Destroy entities, which were allocated by the set itself.
  // Destroy kernels, which were *not* allocated by the set.

  for ( unsigned i = EntityTypeEnd ; 0 < i ; ) {
    --i ;

    EntitySet & eset = m_entities[i] ;
    KernelSet & kset = m_kernels[i] ;

    for ( EntitySet::iterator ie = eset.begin() ; ie != eset.end() ; ++ie ) {
      ie->m_kernel = KernelSet::iterator();
    }
    eset.clear();

    while ( kset.size() ) {
      KernelSet::iterator ik = kset.end();
      --ik ;
      destroy_kernel( ik );
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void BulkData::update_state()
{
  for ( unsigned i = 0 ; i < EntityTypeEnd ; ++i ) {

    KernelSet & kset = m_kernels[ i ] ;

    for ( KernelSet::iterator ik = kset.begin() ; ik != kset.end() ; ++ik ) {
      ik->update_state();
    }
  }
}

//----------------------------------------------------------------------

size_t BulkData::synchronize_changes()
{
  int change = 0 ;

  // Synchronize changes to entities or relations. 
  // For example, if changes are accumulated and processed
  // via batch-oriented call-backs then these call-backs are
  // invoked here.




  // The very last operation performed is to sort the kernel entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and kernels
  // is independent of the order in which a set of changes were
  // performed.

  if ( internal_sort_kernel_entities() ) { change = 1 ; }

  // If any change has occured on any process
  // then increment the change counter.

  all_reduce( m_parallel_machine , Max<1>( & change ) );

  if ( change ) { ++m_sync_change_count ; }

  return m_sync_change_count ;
}

//----------------------------------------------------------------------

const EntitySet & BulkData::entities( unsigned entity_type ) const
{
  const unsigned i = entity_type ;

  if ( EntityTypeEnd <= i ) {
    // Error
    std::string msg( "phdmesh::BulkData::entities FAILED with invalid type" );
    throw std::invalid_argument(msg);
  }

  return m_entities[ entity_type ];
}

const KernelSet & BulkData::kernels( unsigned entity_type ) const
{
  const unsigned i = entity_type ;

  if ( EntityTypeEnd <= i ) {
    // Error
    std::string msg( "phdmesh::BulkData::kernels FAILED with invalid type" );
    throw std::invalid_argument(msg);
  }

  return m_kernels[ entity_type ];
}

//----------------------------------------------------------------------

Entity * BulkData::get_entity( entity_key_type key ,
                                   const char * required_by ) const
{
  const EntitySet & es = entities( phdmesh::entity_type( key ) );
  EntitySet::iterator i = es.find( key );
  if ( required_by && i == es.end() ) {
    static const char method[] = "phdmesh::BulkData::get_entity" ;
    std::ostringstream msg ;
    msg << method << "( " ;
    print_entity_key( msg , key );
    msg << " , " << required_by << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }
  return i != es.end() ? & * i : (Entity*) NULL ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

Entity & BulkData::declare_entity( entity_key_type key ,
                                       const std::vector<Part*> & parts ,
                                       int new_owner )
{
  const char method[] = "phdmesh::BulkData::declare_entity" ;

  const bool bad_key = 0 == entity_id( key );
  const bool bad_own = ((int) m_parallel_size ) <= new_owner ;

  if ( bad_key || bad_own ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( " ;
    print_entity_key( msg , key );
    if ( bad_key ) { msg << " : BAD VALUE" ; }
    msg << " , {" ;
    for ( std::vector<Part*>::const_iterator
          i = parts.begin() ; i != parts.end() ; ++i ) {
      msg << " " << (*i)->name();
    }
    msg << " } , " ;
    msg << new_owner ;
    if ( bad_own ) { msg << " : BAD VALUE" ; }
    msg << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }

  EntitySet & eset = m_entities[ entity_type( key ) ] ;

  const std::pair< EntitySet::iterator , bool > result = eset.insert( key );

  // Determine whether to change the parts,
  // which is an expensive operation.

  const KernelSet::const_iterator k = result.first->m_kernel ;

  std::vector<Part*> add ;
  std::vector<Part*> rem ;

  Part * const owns = & m_mesh_meta_data.locally_owned_part();

  if ( result.second ) { // New entity, must add all parts
    add.reserve( parts.size() + 1 );
    add.assign( parts.begin() , parts.end() );

    // If owner is not given then it is the local processor

    if ( new_owner < 0 ) { new_owner = m_parallel_rank ; }

    if ( new_owner == (int) m_parallel_rank ) { add.push_back( owns ); }
  }
  else { // Existing entity

    if ( ! k->member_all( parts ) ) { // Must add parts
      add.reserve( parts.size() + 1 );
      add.assign( parts.begin() , parts.end() );
    }

    const unsigned current_owner = result.first->m_owner_rank ;

    if ( new_owner < 0 ) { new_owner = current_owner ; }

    if ( new_owner != (int) current_owner ) {
      if ( new_owner == (int) m_parallel_rank ) { // Change to local owner
        add.push_back( owns );
      }
      else { // Change to remote owner
        rem.push_back( owns );
      }
    }
  }

  result.first->m_owner_rank = new_owner ;

  if ( ! add.empty() || ! rem.empty() ) {
    change_entity_parts( * result.first , add , rem );
  }

  return * result.first ;
}

//----------------------------------------------------------------------

void BulkData::change_entity_parts(
  Entity & e ,
  const std::vector<Part*> & add_parts ,
  const std::vector<Part*> & remove_parts )
{
  const char method[] = "phdmesh::BulkData::change_entity_parts" ;

  // Change required if:
  // (1) Entity not a member of a kernel
  // (2) Kernel is not a member of all add_parts
  // (3) Kernel is a member of any remove_parts

  const KernelSet::iterator k = e.m_kernel ;  

  const bool require_change = ( ! k ) ||
                              ( ! k->member_all( add_parts ) ) ||
                              (   k->member_any( remove_parts ) );

  if ( require_change ) {

    // Include supersets of add parts and subsets of removed parts.

    PartSet a_parts( add_parts );
    PartSet r_parts( remove_parts );

    order( a_parts );
    order( r_parts );

    for ( PartSet::const_iterator
          ir = remove_parts.begin(); ir != remove_parts.end() ; ++ir ) {
      Part & rp = **ir ;

      for ( std::vector<Part*>::const_iterator
              jp =  rp.subsets().begin() ;
              jp != rp.subsets().end() ; ++jp ) {
        insert( r_parts , **jp );
      }
    }

    for ( PartSet::const_iterator
          ia = add_parts.begin(); ia != add_parts.end() ; ++ia ) {
      Part & ap = **ia ;

      for ( std::vector<Part*>::const_iterator
              jp =  ap.supersets().begin() ;
              jp != ap.supersets().end() ; ++jp ) {
        insert( a_parts , **jp );
      }
    }

    // Do not allow an intersection-part to explicitly appear

    for ( std::vector<Part*>::const_iterator
          i = a_parts.begin() ; i != a_parts.end() ; ++i ) {
      if ( ! (*i)->intersection_of().empty() ) {
        std::ostringstream msg ;
        msg << method ;
        msg << "( " ;
        print_entity_key( msg , e.key() );
        msg << " , { " ;
        msg << (*i)->name();
        msg << " } , ) FAILED due to explicit use of this intersection-part." ;
        throw std::runtime_error( msg.str() );
      }
    }

    for ( std::vector<Part*>::const_iterator
          i = r_parts.begin() ; i != r_parts.end() ; ++i ) {
      if ( ! (*i)->intersection_of().empty() ) {
        std::ostringstream msg ;
        msg << method ;
        msg << "( " ;
        print_entity_key( msg , e.key() );
        msg << " , , { " ;
        msg << (*i)->name();
        msg << " } ) FAILED due to explicit use of this intersection-part." ;
        throw std::runtime_error( msg.str() );
      }
    }

    // Do not allow a PartRelation::m_target part to appear in the input

    const std::vector<PartRelation> & part_rel = m_mesh_meta_data.get_part_relations();

    for ( std::vector<PartRelation>::const_iterator
          i = part_rel.begin() ; i != part_rel.end() ; ++i ) {
      bool error_add = contain( a_parts , * i->m_target );
      bool error_rem = contain( r_parts , * i->m_target );
      if ( error_add || error_rem ) {
        std::ostringstream msg ;
        msg << method ;
        msg << "( " ;
        print_entity_key( msg , e.key() );
        msg << " , " ;
        if ( error_add ) { msg << "{ " << i->m_target->name() << " }" ; }
        msg << " , " ;
        if ( error_rem ) { msg << "{ " << i->m_target->name() << " }" ; }
        msg << " ) FAILED due to explicit use of this relation-part." ;
        throw std::runtime_error( msg.str() );
      }
    }

    internal_change_entity_parts( e , a_parts , r_parts );
  }
}

//----------------------------------------------------------------------

namespace {

void merge_in( std::vector<unsigned> & vec , const PartSet & parts )
{
  std::vector<unsigned>::iterator i = vec.begin();
  PartSet::const_iterator ip = parts.begin() ;

  for ( ; i != vec.end() && ip != parts.end() ; ++i ) {

    const unsigned ord = (*ip)->mesh_meta_data_ordinal();

    if ( ord <= *i ) {
      if ( ord < *i ) { i = vec.insert( i , ord ); }
      // Now have: ord == *i
      ++ip ;
    }
  }

  for ( ; ip != parts.end() ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    vec.push_back( ord );
  }
}

}

// The 'add_parts' and 'rem_parts' are completely disjoint.

void BulkData::internal_change_entity_parts(
  Entity & e ,
  const PartSet & add_parts ,
  const PartSet & rem_parts )
{
  // Change if:
  // (1) no kernel
  // (2) any add_parts are not in the kernel
  // (3) any rem_parts are in the kernel

  const KernelSet::iterator ik_old = e.m_kernel ;
  const unsigned            i_old  = e.m_kernel_ord ;

  const bool require_change =
    ( ! ik_old ) ||
    ( ! ik_old->member_all( add_parts ) ) ||
    (   ik_old->member_any( rem_parts ) ) ;

  if ( require_change ) {

    Part * const univ_part = & m_mesh_meta_data.universal_part();

    //----------------------------------
#if 0
    // Verify ordered and complete:
    for ( PartSet::const_iterator
          i = add_parts.begin() ; i != add_parts.end() ; ++i ) {

      bool ok = add_parts.begin() == i ||
                (i[-1])->mesh_meta_data_ordinal() < (*i)->mesh_meta_data_ordinal() ;

      const PartSet super = (*i)->supersets();

      for ( PartSet::const_iterator
            ip = super.begin() ; ok && ip != super.end() ; ++ip ) {
        if ( univ_part != *ip ) { ok = contain( add_parts , **ip ); }
      }

      if ( ! ok ) {
        std::ostringstream msg ;
        msg << "internal_change_entity_parts:" ;
        for ( PartSet::const_iterator j =
              add_parts.begin() ; j != add_parts.end() ; ++j ) {
          if ( j == i ) {
            msg << " FAILED( " ;
            msg << " " << (*j)->name() << "[" << (*j)->mesh_meta_data_ordinal() << "]" ;
            msg << " )" ;
          }
          else {
            msg << " " << (*j)->name() << "[" << (*j)->mesh_meta_data_ordinal() << "]" ;
          }
        }
        throw std::logic_error( msg.str() );
      }
    }

    //----------------------------------

    for ( PartSet::const_iterator
          i = rem_parts.begin() ; i != rem_parts.end() ; ++i ) {
      const PartSet & sub = (*i)->subsets();
      if ( ! sub.empty() && ! contain( rem_parts , sub ) ) {
        std::ostringstream msg ;
        msg << "internal_change_entity_parts:" ;
        for ( PartSet::const_iterator
              j = rem_parts.begin() ; j != rem_parts.end() ; ++j ) {
          if ( j == i ) {
            msg << " FAILED( " << (*j)->name() << " )" ;
          }
          else {
            msg << " " << (*j)->name();
          }
        }
        throw std::logic_error( msg.str() );
      }
    }

#endif
    //----------------------------------

    PartSet parts_removed ;

    std::vector<unsigned> parts_total ;

    if ( ! ik_old ) {
      const unsigned univ_ord = univ_part->mesh_meta_data_ordinal();

      parts_total.reserve( add_parts.size() + 1 );

      parts_total.push_back( univ_ord );
    }
    else {
      // Keep any of the existing kernel's parts
      // that are not a remove part.
      // This will include the 'intersection' parts.
      //
      // These parts are properly ordered and unique.

      const std::pair<const unsigned *, const unsigned*>
        kernel_parts = ik_old->superset_part_ordinals();

      const unsigned number_parts = kernel_parts.second - kernel_parts.first ;

      parts_total.reserve( number_parts + add_parts.size() );

      if ( rem_parts.empty() ) {
        parts_total.assign( kernel_parts.first , kernel_parts.second );
      }
      else {
        for ( const unsigned * ik = kernel_parts.first ;
                               ik < kernel_parts.second ; ++ik ) {

          Part * const p = & m_mesh_meta_data.get_part( *ik );

          if ( ! contain( rem_parts , *p ) ) {
            parts_total.push_back( *ik );
          }
          else {
            parts_removed.push_back( p );
          }
        }
      }
    }

    merge_in( parts_total , add_parts );

    //--------------------------------
    // Move the entity to the new kernel.

    KernelSet::iterator ik =
      declare_kernel( e.entity_type(), parts_total.size(), & parts_total[0] );

    // If changing kernels then copy its field values from old to new kernel

    if ( ik_old ) {
      Kernel::copy_fields( *ik , ik->m_size , *ik_old , i_old );
    }
    else {
      Kernel::zero_fields( *ik , ik->m_size );
    }

    // Set the new kernel
    e.m_kernel     = ik ;
    e.m_kernel_ord = ik->m_size ;
    ik->m_entities[ ik->m_size ] = & e ;
    ++( ik->m_size );

    // If changing kernels then remove the entity from the kernel,
    if ( ik_old ) { remove_entity( ik_old , i_old ); }

    // Propagate part changes through the entity's relations.

    internal_propagate_part_changes( e , parts_removed );
  }
}

//----------------------------------------------------------------------

void BulkData::change_entity_identifier( Entity & e , entity_id_type id )
{
  static const char method[] =
    "phdmesh::BulkData::change_entity_identifier" ;

  const bool valid_id = id != 0 ;
  bool ok = valid_id ;

  if ( ok ) {
    const EntityType      type = e.entity_type();
    const entity_key_type key  = entity_key( type , id );
    Entity * const ptr = & e ;

    EntitySet & eset = m_entities[ type ] ;

    const std::pair< EntitySet::iterator , bool >
      result = eset.insert( key , ptr );

    ok = result.second || ( ptr == & * result.first );
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "P" << parallel_size() ;
    msg << ": " << method ;
    msg << "( " ;
    print_entity_key( msg , e.key() );
    msg << " , " ;
    msg << id ;
    msg << " ) FAILED" ;
    if ( valid_id ) { msg << " identifier is in use" ; }
    else            { msg << " identifier is invalid" ; }
    throw std::invalid_argument( msg.str() );
  }
}

void BulkData::change_entity_owner( Entity & e , unsigned owner_rank )
{
  static const char method[] = "phdmesh::BulkData::change_entity_owner" ;

  if ( parallel_size() <= owner_rank ) {
    std::ostringstream msg ;
    msg << "P" << parallel_size() ;
    msg << ": " << method ;
    msg << "( " ;
    print_entity_key( msg , e.key() );
    msg << " , " ;
    msg << owner_rank ;
    msg << " ) FAILED due to invalid rank >= " ;
    msg << parallel_size() ;
    throw std::invalid_argument( msg.str() );
  }
  e.m_owner_rank = owner_rank ;
}

//----------------------------------------------------------------------

void BulkData::destroy_entity( Entity * e )
{
  while ( ! e->m_relation.empty() ) {
    destroy_relation( * e , * e->m_relation.back().entity() );
  }

  remove_entity( e->m_kernel , e->m_kernel_ord );

  e->m_kernel     = KernelSet::iterator();
  e->m_kernel_ord = 0 ;

  const unsigned  entity_type = e->entity_type();

  EntitySet & es = m_entities[ entity_type ];

  es.erase( e );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void verify_set_shares( const BulkData & M )
{
  static const char method[] = "phdmesh::BulkData::set_shares" ;

  std::string msg ;

  const std::vector<EntityProc> & shares = M.shared_entities();

  // Parallel verification

  bool ok = comm_verify( M.parallel() , shares , msg );

  if ( ok ) {

    const unsigned p_rank = M.parallel_rank();
    Part & owns_part = M.mesh_meta_data().locally_owned_part();
    Part & uses_part = M.mesh_meta_data().locally_used_part();

    std::ostringstream os ;

    os << "P" << p_rank << ": " << method << " FAILED " << std::endl ;

    const std::vector<EntityProc>::const_iterator es = shares.end();
          std::vector<EntityProc>::const_iterator is = shares.begin();

    for ( ; is != es ; ++is ) {

      // Verify not attempting to share with self

      if ( p_rank == is->second ) {
        print_entity_key( os , is->first->key() );
        os << " paired with P" << p_rank ;
        os << std::endl ;
        ok = false ;
      }

      // Verify each member has the uses part

      if ( ! is->first->kernel().has_superset( uses_part ) ) {
        print_entity_key( os , is->first->key() );
        os << " to be shared with P" ;
        os << is->second ;
        os << " does not have part " ;
        os << uses_part.name();
        os << std::endl ;
        ok = false ;
      }
    }

    // Verify all uses but not owned entities are shared
    // and their owner is one of the shared.

    for ( unsigned t = 0 ; t < EntityTypeEnd ; ++t ) {

      const KernelSet & kernels = M.kernels( t );

      const KernelSet::iterator ek = kernels.end();
            KernelSet::iterator ik = kernels.begin();

      for ( ; ek != ik ; ++ik ) {

        if ( ik->has_superset( uses_part ) ) {

          const Kernel::iterator e = ik->end();
                Kernel::iterator i = ik->begin();

          if ( ik->has_superset( owns_part ) ) {

            for ( ; e != i ; ++i ) {

              if ( (*i)->owner_rank() != p_rank ) {
                print_entity_key( os , (*i)->key() ); 
                os << " member of " ;
                os << owns_part.name();
                os << " but non-local owner = P" ;
                os << (*i)->owner_rank() ;
                os << std::endl ;
                ok = false ;
              }
            }
          }
          else {

            for ( ; e != i ; ++i ) {

              if ( (*i)->owner_rank() == p_rank ) {
                print_entity_key( os , (*i)->key() ); 
                os << " not member of " ;
                os << owns_part.name();
                os << " but local owner = P" ;
                os << p_rank ;
                os << std::endl ;
                ok = false ;
              }

              PairIterEntityProc ss = (*i)->sharing();

              if ( ss.empty() ) {
                print_entity_key( os , (*i)->key() );
                os << " not member of " ;
                os << owns_part.name();
                os << ", is member of " ;
                os << uses_part.name();
                os << ", but is not shared." ;
                os << std::endl ;
                msg.append( os.str() );
                ok = false ;
              }

              for ( ; ss && (*i)->owner_rank() != ss->second ; ++ss );

              if ( ss.empty() ) {
                print_entity_key( os , (*i)->key() );
                os << " not member of " ;
                os << owns_part.name();
                os << ", is member of " ;
                os << uses_part.name();
                os << ", but does not share with owner P" ;
                os << (*i)->owner_rank();
                os << std::endl ;
                msg.append( os.str() );
                ok = false ;
              }
            }
          } // end not member of owns_spart
        }   // end not member of uses_part
      }     // end kernel loop
    }       // end entity type loop

    if ( ! ok ) { msg = os.str(); }
  }

  {
    // Global reduce on the error flag, are any of the flag false ?
    unsigned flag = ok ;
    all_reduce( M.parallel() , Min<1>( & flag ) );
    ok = flag ;
  }

  if ( ! ok ) {
    throw std::runtime_error( msg );
  }
}

}

void BulkData::set_shares( const std::vector<EntityProc> & s )
{
  m_shares_all = s ;

  // Clear the entities' sharing in preparation for setting it.

  {
    const PairIterEntityProc ss_empty ;

    for ( unsigned t = 0 ; t < EntityTypeEnd ; ++t ) {
      const EntitySet::iterator e = m_entities[t].end();
            EntitySet::iterator i = m_entities[t].begin();
      for ( ; e != i ; ++i ) { i->m_sharing = ss_empty ; }
    }
  }

  // Set the entities' sharing.

  const std::vector<EntityProc>::iterator es = m_shares_all.end();
        std::vector<EntityProc>::iterator is = m_shares_all.begin();

  while ( is != es ) {
    const std::vector<EntityProc>::iterator js = is ;
    for ( ; is != es && js->first == is->first ; ++is );
    js->first->m_sharing = PairIterEntityProc( js , is );
  }

  // Verify 

  verify_set_shares( *this );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void verify_set_ghosting( const BulkData & M )
{
  static const char method[] = "phdmesh::BulkData::set_ghosting" ;

  std::string msg ;

  const std::vector<EntityProc> & domain = M.ghost_source();
  const std::vector<EntityProc> & range  = M.ghost_destination();

  // Parallel verification

  bool ok = comm_verify( M.parallel(), domain , range , msg );

  // Local verification of range ordering

  if ( ok ) { ok = verify( range , msg ); }

  if ( ok ) {

    Part & uses_part = M.mesh_meta_data().locally_used_part();
    Part & owns_part = M.mesh_meta_data().locally_owned_part();

    const unsigned p_rank = M.parallel_rank();

    std::ostringstream os ;

    os << "P" << M.parallel_rank() << ": " << method << " FAILED " ;

    // Verify each entity in the domain has the owns part
    {
      Entity * e = NULL ;

      const std::vector<EntityProc>::const_iterator es = domain.end();
            std::vector<EntityProc>::const_iterator is = domain.begin();

      for ( ; is != es ; ++is ) {
        if ( e != is->first ) {
          e = is->first ;
        }

        if ( ! e->kernel().has_superset( owns_part ) ||
             p_rank != e->owner_rank() ) {
          print_entity_key( os , is->first->key() );
          os << " does not have part " ;
          os << owns_part.name();
          os << std::endl ;
          msg = os.str();
          ok = false ;
        }
      }
    }

    // Verify
    // (1) each entity in the range does not have the uses part
    // (2) an entity only appears in the range once

    {
      unsigned entity_count[ EntityTypeEnd ];
      unsigned kernel_count[ EntityTypeEnd ];

      for ( unsigned i = 0 ; i < EntityTypeEnd ; ++i ) {
        entity_count[i] = 0 ;
        kernel_count[i] = 0 ;
      }

      Entity * e = NULL ;

      const std::vector<EntityProc>::const_iterator es = range.end();
            std::vector<EntityProc>::const_iterator is = range.begin();

      for ( ; is != es ; ++is ) {
        const bool bad_domain = e == is->first ;
        e = is->first ;

        ++( entity_count[ e->entity_type() ] );

        const bool bad_part = e->kernel().has_superset( uses_part );
        const bool bad_rank = is->second != e->owner_rank();

        if ( bad_part || bad_rank || bad_domain ) {
          print_entity_key( os , e->key() );
          if ( bad_part ) {
            os << " member of " ;
            os << uses_part.name();
          }
          if ( bad_domain ) {
            os << " Received from multiple procesors " ;
          }
          if ( bad_rank ) {
            os << " Received from P" ;
            os << is->second ;
            os << " instead of owner P" ;
            os << e->owner_rank();
          }
          msg = os.str();
          ok = false ;
        }
      }

      // Verify aura range count equals not member of uses count.

      for ( unsigned i = 0 ; ok && i < EntityTypeEnd ; ++i ) {
        const KernelSet & kset = M.kernels( i );
        const KernelSet::const_iterator ek = kset.end();
              KernelSet::const_iterator ik = kset.begin();
        for ( ; ik != ek ; ++ik ) {
          if ( ! ik->has_superset( uses_part ) ) {
            kernel_count[i] += ik->size();
          }
        }

        if ( entity_count[i] != kernel_count[i] ) {
          os << " aura entity count = " << entity_count[i] ;
          os << " != " << kernel_count[i] << " = not 'uses' entity count" ;
          msg = os.str();
          ok = false ;
        }
      }
    }
  }

  {
    // Global reduce on the error flag, are any of the flag false ?
    unsigned flag = ok ;
    all_reduce( M.parallel() , Min<1>( & flag ) );
    ok = flag ;
  }

  if ( ! ok ) {
    throw std::runtime_error( msg );
  }

}

}

void BulkData::set_ghosting(
  const std::vector<EntityProc> & d ,
  const std::vector<EntityProc> & r )
{
  m_aura_domain = d ;
  m_aura_range  = r ;

  verify_set_ghosting( *this );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void get_kernels( const BulkData & mesh ,
                  EntityType type ,
                  Part & part ,
                  std::vector<const Kernel*> & kernels )
{
  const KernelSet & ks = mesh.kernels( type );
  kernels.clear();

  const KernelSet::const_iterator ie = ks.end();
        KernelSet::const_iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    const Kernel * const k = & * ik ;
    if ( k->has_superset( part ) ) {
      kernels.push_back( k );
    }
  }
}

void get_kernels_intersect(
  const BulkData & mesh ,
  EntityType type ,
  const std::vector<Part*> & parts ,
  std::vector<const Kernel*> & kernels )
{
  const KernelSet & ks = mesh.kernels( type );
  kernels.clear();

  const KernelSet::const_iterator ie = ks.end();
        KernelSet::const_iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    const Kernel * const k = & * ik ;
    if ( k->has_superset( parts ) ) {
      kernels.push_back( k );
    }
  }
}

void get_kernels_union(
  const BulkData & mesh ,
  EntityType type ,
  const std::vector<Part*> & parts ,
  std::vector<const Kernel*> & kernels )
{
  const KernelSet & ks = mesh.kernels( type );
  kernels.clear();

  const KernelSet::const_iterator ie = ks.end();
        KernelSet::const_iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    const Kernel * const k = & * ik ;
    bool result = false ;
    const std::vector<Part*>::const_iterator ep = parts.end();
          std::vector<Part*>::const_iterator ip = parts.begin();
    for ( ; ! result && ip != ep ; ++ip ) {
      result = k->has_superset( **ip );
    }
    if ( result ) {
      kernels.push_back( k );
    }
  }
}

//----------------------------------------------------------------------

void count_entities(
  BulkData & mesh , Part & part , unsigned * const count )
{
  static const char method[] = "phdmesh::count_entities" ;

  mesh.mesh_meta_data().assert_same_mesh_meta_data( method , part.mesh_meta_data() );

  for ( unsigned i = 0 ; i < EntityTypeEnd ; ++i ) {
    count[i] = 0 ;

    const KernelSet & ks = mesh.kernels( i );

    KernelSet::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( ik->has_superset( part ) ) {
        count[i] += ik->size();
      }
    }
  }
}

void count_entities(
  BulkData & mesh , const PartSet & parts , unsigned * const count )
{
  static const char method[] = "phdmesh::count_entities" ;

  for ( unsigned i = 0 ; i < EntityTypeEnd ; ++i ) {
    count[i] = 0 ;
  }

  if ( ! parts.empty() ) {
    PartSet tmp( parts );

    order( tmp );

    {
      const PartSet::iterator j = tmp.end();
            PartSet::iterator i = tmp.begin();
      for ( ; i != j ; ++i ) {
        mesh.mesh_meta_data().assert_same_mesh_meta_data( method , (*i)->mesh_meta_data() );
      }
    }

    for ( unsigned i = 0 ; i < EntityTypeEnd ; ++i ) {
      const KernelSet & ks = mesh.kernels( i );

      KernelSet::const_iterator ik ;

      for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
        if ( ik->has_superset( tmp ) ) {
          count[i] += ik->size();
        }
      }
    }
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

