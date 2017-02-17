/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

//----------------------------------------------------------------------

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Trace.hpp>

//----------------------------------------------------------------------

namespace stk_classic {
namespace mesh {

bool comm_mesh_verify_parallel_consistency(
  BulkData & M , std::ostream & error_log );

//----------------------------------------------------------------------

unsigned BulkData::determine_new_owner( Entity & entity ) const
{
  // We will decide the new owner by looking at all the processes sharing
  // this entity. The new owner will be the sharing process with lowest rank.

  // The local process is a candidate only if the entity is not destroyed.
  unsigned new_owner =
    EntityLogDeleted == entity.log_query() ? ~0u : m_parallel_rank ;

  for ( PairIterEntityComm
        share = m_entity_comm_map.sharing(entity.key()); ! share.empty() ; ++share ) {
    if ( share->proc < m_parallel_size &&
         ( new_owner < share->proc || m_parallel_size <= new_owner ) ) {
      new_owner = share->proc ;
    }
  }

  return new_owner ;
}

//----------------------------------------------------------------------

namespace {

// A method for quickly finding an entity
Entity* find_entity(const EntityVector& entities, const EntityKey& key,
                    bool expect_success = false)
{
  EntityVector::const_iterator itr =
    std::lower_bound(entities.begin(),
                     entities.end(),
                     key,
                     EntityLess());
  if (itr == entities.end() || (*itr)->key() != key) {
    ThrowRequireMsg(!expect_success,
                    "Expected to be able to find entity of type: " <<
                    key.type() << " and rank: " << key.rank());
    return NULL;
  }
  return *itr;
}

struct EntityProcState {
  EntityProc entity_proc;
  EntityModificationLog state;

  bool operator<(const EntityProcState& rhs) const
  {
    EntityLess el;
    return el(entity_proc, rhs.entity_proc);
  }
};

bool pack_entity_modification( const BulkData & mesh ,
                               const bool pack_shared ,
                               CommAll & comm )
{
  bool flag = false ;

  const std::vector<Entity*> & entity_comm = mesh.entity_comm();

  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity & entity = **i ;

    if ( entity.log_query() == EntityLogModified ||
         entity.log_query() == EntityLogDeleted ) {

      for ( PairIterEntityComm ec = mesh.entity_comm(entity.key()); ! ec.empty() ; ++ec ) {
        const bool shared = 0 == ec->ghost_id ;
        if ( pack_shared == shared ) {
          comm.send_buffer( ec->proc )
              .pack<EntityKey>( entity.key() )
              .pack<EntityModificationLog>( entity.log_query() );

          flag = true ;
        }
      }
    }
  }

  return flag ;
}

void communicate_entity_modification( const BulkData & mesh ,
                                      const bool shared ,
                                      std::vector<EntityProcState > & data )
{
  CommAll comm( mesh.parallel() );

  // Sizing send buffers:
  const bool local_mod = pack_entity_modification( mesh , shared , comm );

  // Allocation of send and receive buffers:
  const bool global_mod =
    comm.allocate_buffers( comm.parallel_size() / 4 , false , local_mod );

  if ( global_mod ) {
    const std::vector<Entity*> & entity_comm = mesh.entity_comm();

    // Packing send buffers:
    pack_entity_modification( mesh , shared , comm );

    comm.communicate();

    for ( unsigned p = 0 ; p < comm.parallel_size() ; ++p ) {
      CommBuffer & buf = comm.recv_buffer( p );
      EntityKey key ;
      EntityProcState tmp ;

      while ( buf.remaining() ) {

        buf.unpack<EntityKey>( key )
           .unpack<EntityModificationLog>( tmp.state );

        // search through entity_comm, should only receive info on entities
        // that are communicated.
        tmp.entity_proc.first  = find_entity(entity_comm, key, true);
        tmp.entity_proc.second = p ;

        data.push_back( tmp );
      }
    }
  }

  std::sort( data.begin() , data.end() );
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

// Postconditions:
//  * DistributedIndex is updated based on entity creation/deletions in the
//    last modification cycle.
//  * Comm lists for shared entities are up-to-date.
//  * shared_new contains all entities that were modified/created on a
//    different process
void BulkData::internal_update_distributed_index(
  std::vector<Entity*> & shared_new )
{
  Trace_("stk_classic::mesh::BulkData::internal_update_distributed_index");

  std::vector< parallel::DistributedIndex::KeyType >
    local_created_or_modified , // only store locally owned/shared entities
    del_entities_keys ;

  // Iterate over all entities known to this process, putting
  // locally deleted entities in del_entities_keys, and putting
  // modified shared/owned entities in local_created_or_modified.
  for ( impl::EntityRepository::iterator
        i = m_entity_repo.begin() ; i != m_entity_repo.end() ; ++i ) {

    Entity & entity = * i->second ;

    if ( EntityLogDeleted == entity.log_query() ) {
      // Has been destroyed
      del_entities_keys.push_back( entity.key().raw_key() );
    }
    else if ( entity.log_query() != EntityLogNoChange &&
              in_owned_closure( entity , m_parallel_rank ) ) {
      // Has been changed and is in owned closure, may be shared
      local_created_or_modified.push_back( entity.key().raw_key() );
    }
  }

  // Update distributed index. Note that the DistributedIndex only
  // tracks ownership and sharing information.
  m_entities_index.update_keys( local_created_or_modified , del_entities_keys );

  if (parallel_size() > 1) {
    // Retrieve data regarding which processes use the local_created_or_modified
    // including this process.
    std::vector< parallel::DistributedIndex::KeyProc >
      global_created_or_modified ;
    m_entities_index.query_to_usage( local_created_or_modified ,
                                     global_created_or_modified );

    //------------------------------
    // Take the usage data and update the sharing comm lists
    {
      Entity * entity = NULL ;

      // Iterate over all global modifications to this entity, this vector is
      // sorted, so we're guaranteed that all modifications to a particular
      // entities will be adjacent in this vector.
      for ( std::vector< parallel::DistributedIndex::KeyProc >::iterator
              i =  global_created_or_modified.begin() ;
            i != global_created_or_modified.end() ; ++i ) {

        EntityKey key( & i->first );
        unsigned modifying_proc = i->second;

        // key should not be in del_entities_keys
        ThrowAssertMsg( !std::binary_search(del_entities_keys.begin(),
                                            del_entities_keys.end(),
                                            i->first),
                        "Key: " << print_entity_key(mesh_meta_data(), key) <<
                        " was locally deleted, but somehow was included in global_created_or_modified; " <<
                        " this probably means there's problem in DistributedIndex." );

        if ( m_parallel_rank != modifying_proc ) {
          // Another process also created or updated this entity.

          // Only want to look up entities at most once
          if ( entity == NULL || entity->key() != key ) {
            // Have not looked this entity up by key
            entity = get_entity( key );

            shared_new.push_back( entity );
          }

          // Add the other_process to the entity's sharing info.
          m_entity_comm_map.insert(entity->key(), EntityCommInfo( 0, // sharing
                                                                   modifying_proc ) );
        }
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

// Enforce that shared entities must be in the owned closure:

void destroy_dependent_ghosts( BulkData & mesh , Entity * entity )
{
  for ( ; ; ) {
    PairIterRelation rel = entity->relations();

    if ( rel.empty() ) { break ; }

    Entity * e = rel.back().entity();

    if ( e->entity_rank() < entity->entity_rank() ) { break ; }

    ThrowRequireMsg( !in_owned_closure( *e , mesh.parallel_rank()),
        "Entity " << print_entity_key(e) << " should not be in closure." );

    destroy_dependent_ghosts( mesh , e );
  }

  mesh.destroy_entity( entity );
}

// Entities with sharing information that are not in the owned closure
// have been modified such that they are no longer shared.
// These may no longer be needed or may become ghost entities.
// There is not enough information so assume they are to be deleted
// and let these entities be re-ghosted if they are needed.

// Open question: Should an owned and shared entity that does not
// have an upward relation to an owned entity be destroyed so that
// ownership transfers to another process?

void resolve_shared_removed_from_owned_closure( BulkData & mesh )
{
  for ( std::vector<Entity*>::const_reverse_iterator
        i =  mesh.entity_comm().rbegin() ;
        i != mesh.entity_comm().rend() ; ++i) {

    Entity * entity = *i ;

    if ( ! mesh.entity_comm_sharing(entity->key()).empty() &&
         ! in_owned_closure( *entity , mesh.parallel_rank() ) ) {

      destroy_dependent_ghosts( mesh , entity );
    }
  }
}

}

// Resolve modifications for shared entities:
// If not locally destroyed and remotely modified
// then set to locally modified.
// If remotely destroyed then determine the new owner.
//
// Post condition:
//  Shared entities are in-sync with respect to modification state.
//  Shared communication lists are updated to reflect all deletions.
//  Ownership has been re-assigned as necessary for deletion
//  of shared entities.

void BulkData::internal_resolve_shared_modify_delete()
{
  Trace_("stk_classic::mesh::BulkData::internal_resolve_shared_modify_delete");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  resolve_shared_removed_from_owned_closure( *this );

  std::vector< EntityProcState > remote_mod ;

  // Communicate entity modification state for shared entities
  // the resulting vector is sorted by entity and process.
  const bool communicate_shared = true ;
  communicate_entity_modification( *this , communicate_shared , remote_mod );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first.
  for ( std::vector<EntityProcState>::reverse_iterator
        i = remote_mod.rbegin(); i != remote_mod.rend() ; ) {

    Entity * const entity        = i->entity_proc.first ;
    const bool locally_destroyed = EntityLogDeleted == entity->log_query();
    bool remote_owner_destroyed  = false;

    // Iterate over all of this entity's remote changes
    for ( ; i != remote_mod.rend() && i->entity_proc.first == entity ; ++i ) {

      const unsigned remote_proc    = i->entity_proc.second ;
      const bool remotely_destroyed = EntityLogDeleted == i->state ;

      // When a shared entity is remotely modified or destroyed
      // then the local copy is also modified.  This modification
      // status is applied to all related higher ranking entities.

      if ( ! locally_destroyed ) {
        m_entity_repo.log_modified( *entity );
      }

      // A shared entity is being deleted on the remote process.
      // Remove it from the sharing communication list.
      // Ownership changes are processed later, but we'll need
      // to know if the remote owner destroyed the entity in order
      // to correctly resolve ownership (it is not sufficient to just
      // look at the comm list of the entity since there is no
      // guarantee that the comm list is correct or up-to-date).

      if ( remotely_destroyed ) {
        m_entity_comm_map.erase(entity->key(), EntityCommInfo(0,remote_proc) );

        // check if owner is destroying
        if ( entity->owner_rank() == remote_proc ) {
          remote_owner_destroyed = true ;
        }
      }
    }

    // Have now processed all remote changes knowledge for this entity.

    PairIterEntityComm new_sharing = m_entity_comm_map.sharing(entity->key());
    const bool   exists_somewhere = ! ( remote_owner_destroyed &&
                                        locally_destroyed &&
                                        new_sharing.empty() );

    // If the entity has been deleted everywhere, nothing left to do
    if ( exists_somewhere ) {

      const bool old_local_owner = m_parallel_rank == entity->owner_rank();

      // Giving away ownership to another process in the sharing list:
      const bool give_ownership = locally_destroyed && old_local_owner ;

      // If we are giving away ownership or the remote owner destroyed
      // the entity, then we need to establish a new owner
      if ( give_ownership || remote_owner_destroyed ) {

        const unsigned new_owner = determine_new_owner( *entity );

        m_entity_repo.set_entity_owner_rank( *entity, new_owner );
        m_entity_repo.set_entity_sync_count( *entity, m_sync_count );
      }

      if ( ! locally_destroyed ) {

        PartVector add_part , remove_part ;

        if ( new_sharing.empty() ) {
          // Is no longer shared, remove the shared part.
          remove_part.push_back(& m_mesh_meta_data.globally_shared_part());
        }

        const bool new_local_owner = m_parallel_rank == entity->owner_rank();

        const bool local_claimed_ownership =
          ( ! old_local_owner && new_local_owner );

        if ( local_claimed_ownership ) {
          // Changing remotely owned to locally owned
          add_part.push_back( & m_mesh_meta_data.locally_owned_part() );
        }

        if ( ! add_part.empty() || ! remove_part.empty() ) {
          internal_change_entity_parts( *entity , add_part , remove_part );
        }
      } // if ( ! locally_destroyed )
    } // if ( exists_somewhere )
  } // remote mod loop

  // Erase all sharing communication lists for Destroyed entities:
  for ( std::vector<Entity*>::const_reverse_iterator
        i = entity_comm().rbegin() ; i != entity_comm().rend() ; ++i) {
    Entity * entity = *i ;

    if ( EntityLogDeleted == entity->log_query() ) {
      // m_ghosting[0] is the SHARED communication
      m_entity_comm_map.erase(entity->key(), *m_ghosting[0] );
    }
  }
}



//----------------------------------------------------------------------
// Resolve modifications for ghosted entities:
// If a ghosted entity is modified or destroyed on the owning
// process then the ghosted entity must be destroyed.
//
// Post condition:
//  Ghosted entities of modified or deleted entities are destroyed.
//  Ghosted communication lists are cleared to reflect all deletions.

void BulkData::internal_resolve_ghosted_modify_delete()
{
  Trace_("stk_classic::mesh::BulkData::internal_resolve_ghosted_modify_delete");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");
  // Resolve modifications for ghosted entities:

  std::vector<EntityProcState > remote_mod ;

  // Communicate entity modification state for ghost entities
  const bool communicate_shared = false ;
  communicate_entity_modification( *this , communicate_shared , remote_mod );

  const size_t ghosting_count = m_ghosting.size();

  std::vector< int > ghosting_change_flags( ghosting_count , 0 );

  // We iterate backwards over remote_mod to ensure that we hit the
  // higher-ranking entities first. This is important because higher-ranking
  // entities like element must be deleted before the nodes they have are
  // deleted.
  for ( std::vector<EntityProcState>::reverse_iterator
        i = remote_mod.rbegin(); i != remote_mod.rend() ; ++i ) {
    Entity *       entity       = i->entity_proc.first ;
    const unsigned remote_proc  = i->entity_proc.second ;
    const bool     local_owner  = entity->owner_rank() == m_parallel_rank ;
    const bool remotely_destroyed = EntityLogDeleted == i->state ;
    const bool locally_destroyed  = EntityLogDeleted == entity->log_query();

    if ( local_owner ) { // Sending to 'remote_proc' for ghosting

      if ( remotely_destroyed ) {

        // remove from ghost-send list

        for ( size_t j = ghosting_count ; j-- ; ) {
          if ( m_entity_comm_map.erase( entity->key(), EntityCommInfo( j , remote_proc ) ) ) {
            ghosting_change_flags[ j ] = true ;
          }
        }
      }

      // Remotely modified ghosts are ignored

    }
    else { // Receiving from 'remote_proc' for ghosting

      // Owner modified or destroyed, must locally destroy.

      for ( PairIterEntityComm ec = m_entity_comm_map.comm(entity->key()) ; ! ec.empty() ; ++ec ) {
        ghosting_change_flags[ ec->ghost_id ] = true ;
      }

      // This is a receive ghost so the only communication information
      // is the ghosting information, can clear it all out.
      m_entity_comm_map.comm_clear(entity->key());

      if ( ! locally_destroyed ) {

        // If mesh modification causes a ghost entity to become
        // a member of an owned-closure then do not automatically
        // destroy it.  The new sharing status will be resolved
        // in 'internal_resolve_parallel_create'.

        if ( ! in_owned_closure( *entity , m_parallel_rank ) ) {

          const bool destroy_entity_successful = destroy_entity(entity);
          ThrowRequireMsg(destroy_entity_successful,
              "Could not destroy ghost entity " << print_entity_key(entity));
        }
      }
    }
  } // end loop on remote mod

  // Erase all ghosting communication lists for:
  // 1) Destroyed entities.
  // 2) Owned and modified entities.

  for ( std::vector<Entity*>::const_reverse_iterator
        i = entity_comm().rbegin() ; i != entity_comm().rend() ; ++i) {

    Entity & entity = **i ;

    const bool locally_destroyed = EntityLogDeleted == entity.log_query();
    const bool locally_owned_and_modified =
      EntityLogModified == entity.log_query() &&
      m_parallel_rank   == entity.owner_rank() ;

    if ( locally_destroyed || locally_owned_and_modified ) {

      // m_ghosting[0] is the SHARED communication

      for ( size_t j = ghosting_count ; j-- ; ) {
        if ( m_entity_comm_map.erase( entity.key(), *m_ghosting[j] ) ) {
          ghosting_change_flags[ j ] = true ;
        }
      }
    }
  }

  std::vector< int > ghosting_change_flags_global( ghosting_count , 0 );

  all_reduce_sum( m_parallel_machine ,
                  & ghosting_change_flags[0] ,
                  & ghosting_change_flags_global[0] ,
                  ghosting_change_flags.size() );

  for ( unsigned ic = 0 ; ic < ghosting_change_flags_global.size() ; ++ic ) {
    if ( ghosting_change_flags_global[ic] ) {
      m_ghosting[ic]->m_sync_count = m_sync_count ;
    }
  }
}

//----------------------------------------------------------------------

// Postconditions:
//  * All shared entities have parallel-consistent owner
//  * Part membership of shared entities is up-to-date
//  * m_entity_comm is up-to-date
void BulkData::internal_resolve_parallel_create()
{
  Trace_("stk_classic::mesh::BulkData::internal_resolve_parallel_create");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  std::vector<Entity*> shared_modified ;

  // Update the parallel index and
  // output shared and modified entities.
  internal_update_distributed_index( shared_modified );

  // ------------------------------------------------------------
  // Claim ownership on all shared_modified entities that I own
  // and which were not created in this modification cycle. All
  // sharing procs will need to be informed of this claim.
  CommAll comm_all( m_parallel_machine );

  for ( int phase = 0; phase < 2; ++phase ) {
    for ( std::vector<Entity*>::iterator
            i = shared_modified.begin() ; i != shared_modified.end() ; ++i ) {
      Entity & entity = **i ;
      if ( entity.owner_rank() == m_parallel_rank &&
           entity.log_query()  != EntityLogCreated ) {

        for ( PairIterEntityComm
                jc = m_entity_comm_map.sharing(entity.key()) ; ! jc.empty() ; ++jc ) {
          comm_all.send_buffer( jc->proc ) .pack<EntityKey>( entity.key() );
        }
      }
    }

    if (phase == 0) { //allocation phase
      comm_all.allocate_buffers( m_parallel_size / 4 );
    }
    else { // communication phase
      comm_all.communicate();
    }
  }

  for ( unsigned p = 0 ; p < m_parallel_size ; ++p ) {
    CommBuffer & buf = comm_all.recv_buffer( p );
    EntityKey key ;
    while ( buf.remaining() ) {
      buf.unpack<EntityKey>( key );

      Entity & entity = * get_entity( key );

      // Set owner, will correct part membership later
      m_entity_repo.set_entity_owner_rank( entity, p);
    }
  }

  // ------------------------------------------------------------
  // Update shared created entities.
  // - Revise ownership to selected processor
  // - Update sharing.
  // - Work backward so the 'in_owned_closure' function
  //   can evaluate related higher ranking entities.

  std::ostringstream error_msg ;
  int error_flag = 0 ;

  PartVector shared_part , owned_part ;
  shared_part.push_back( & m_mesh_meta_data.globally_shared_part() );
  owned_part.push_back(  & m_mesh_meta_data.locally_owned_part() );

  std::vector<Entity*>::const_reverse_iterator iend = shared_modified.rend();
  for ( std::vector<Entity*>::const_reverse_iterator
        i = shared_modified.rbegin() ; i != iend ; ++i) {

    Entity * entity = *i ;

    if ( entity->owner_rank() == m_parallel_rank &&
         entity->log_query() == EntityLogCreated ) {

      // Created and not claimed by an existing owner

      const unsigned new_owner = determine_new_owner( *entity );

      m_entity_repo.set_entity_owner_rank( *entity, new_owner);
    }

    if ( entity->owner_rank() != m_parallel_rank ) {
      // Do not own it and still have it.
      // Remove the locally owned, add the globally_shared
      m_entity_repo.set_entity_sync_count( *entity, m_sync_count);
      internal_change_entity_parts( *entity , shared_part /*add*/, owned_part /*remove*/);
    }
    else if ( ! m_entity_comm_map.sharing(entity->key()).empty() ) {
      // Own it and has sharing information.
      // Add the globally_shared
      internal_change_entity_parts( *entity , shared_part /*add*/, PartVector() /*remove*/ );
    }
    else {
      // Own it and does not have sharing information.
      // Remove the globally_shared
      internal_change_entity_parts( *entity , PartVector() /*add*/, shared_part /*remove*/);
    }

    // Newly created shared entity had better be in the owned closure
    if ( ! in_owned_closure( *entity , m_parallel_rank ) ) {
      if ( 0 == error_flag ) {
        error_flag = 1 ;
        error_msg
          << "\nP" << m_parallel_rank << ": " << " FAILED\n"
          << "  The following entities were declared on multiple processors,\n"
          << "  cannot be parallel-shared, and were declared with"
          << "  parallel-ghosting information. {\n";
      }
      error_msg << "    " << print_entity_key(entity);
      error_msg << " also declared on" ;
      for ( PairIterEntityComm ec = entity->sharing(); ! ec.empty() ; ++ec ) {
        error_msg << " P" << ec->proc ;
      }
      error_msg << "\n" ;
    }
  }

  // Parallel-consistent error checking of above loop
  if ( error_flag ) { error_msg << "}\n" ; }
  all_reduce( m_parallel_machine , ReduceMax<1>( & error_flag ) );
  ThrowErrorMsgIf( error_flag, error_msg.str() );

  // ------------------------------------------------------------
  // Update m_entity_comm based on shared_modified

  const size_t n_old = m_entity_comm.size();

  m_entity_comm.insert( m_entity_comm.end() ,
                        shared_modified.begin() , shared_modified.end() );

  std::inplace_merge( m_entity_comm.begin() ,
                      m_entity_comm.begin() + n_old ,
                      m_entity_comm.end() ,
                      EntityLess() );

  {
    std::vector<Entity*>::iterator i =
      std::unique( m_entity_comm.begin() , m_entity_comm.end() );

    m_entity_comm.erase( i , m_entity_comm.end() );
  }
}

//----------------------------------------------------------------------

bool BulkData::modification_end()
{
  Trace_("stk_classic::mesh::BulkData::modification_end");

  return internal_modification_end( true );
}

#if 0

namespace {

// Very, very handy for debugging parallel resolution...

void print_comm_list( const BulkData & mesh , bool doit )
{
  if ( doit ) {
    std::ostringstream msg ;

    msg << std::endl ;

    for ( std::vector<Entity*>::const_iterator
          i =  mesh.entity_comm().begin() ;
          i != mesh.entity_comm().end() ; ++i ) {

      Entity & entity = **i ;
      msg << "P" << mesh.parallel_rank() << ": " ;

      print_entity_key( msg , MetaData::get(mesh) , entity.key() );

      msg << " owner(" << entity.owner_rank() << ")" ;

      if ( EntityLogModified == entity.log_query() ) { msg << " mod" ; }
      else if ( EntityLogDeleted == entity.log_query() ) { msg << " del" ; }
      else { msg << "    " ; }

      for ( PairIterEntityComm ec = mesh.entity_comm(entity.key()); ! ec.empty() ; ++ec ) {
        msg << " (" << ec->ghost_id << "," << ec->proc << ")" ;
      }
      msg << std::endl ;
    }

    std::cout << msg.str();
  }
}

}

#endif

bool BulkData::internal_modification_end( bool regenerate_aura )
{
  Trace_("stk_classic::mesh::BulkData::internal_modification_end");

  if ( m_sync_state == SYNCHRONIZED ) { return false ; }

  if (parallel_size() > 1) {
    // Resolve modification or deletion of shared entities
    // which can cause deletion of ghost entities.
    internal_resolve_shared_modify_delete();

    // Resolve modification or deletion of ghost entities
    // by destroying ghost entities that have been touched.
    internal_resolve_ghosted_modify_delete();

    // Resolution of shared and ghost modifications can empty
    // the communication information for entities.
    // If there is no communication information then the
    // entity must be removed from the communication list.
    {
      std::vector<Entity*>::iterator i = m_entity_comm.begin();
      bool changed = false ;
      for ( ; i != m_entity_comm.end() ; ++i ) {
        if ( m_entity_comm_map.comm((*i)->key()).empty() ) { *i = NULL ; changed = true ; }
      }
      if ( changed ) {
        i = std::remove( m_entity_comm.begin() ,
                         m_entity_comm.end() , (Entity *) NULL );
        m_entity_comm.erase( i , m_entity_comm.end() );
      }
    }

    // Resolve creation of entities: discover sharing and set unique ownership.
    internal_resolve_parallel_create();

    // Resolve part membership for shared entities.
    // This occurs after resolving creation so created and shared
    // entities are resolved along with previously existing shared entities.
    internal_resolve_shared_membership();

    // Regenerate the ghosting aura around all shared mesh entities.
    if ( regenerate_aura ) { internal_regenerate_shared_aura(); }

    // ------------------------------
    // Verify parallel consistency of mesh entities.
    // Unique ownership, communication lists, sharing part membership,
    // application part membership consistency.
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = comm_mesh_verify_parallel_consistency( *this , msg );
    ThrowErrorMsgIf( !is_consistent, msg.str() );
  }
  else {
    std::vector<Entity*> shared_modified ;
    internal_update_distributed_index( shared_modified );
  }

  // ------------------------------
  // The very last operation performed is to sort the bucket entities.
  // This does not change the entities, relations, or field data.
  // However, it insures that the ordering of entities and buckets
  // is independent of the order in which a set of changes were
  // performed.
  //
  //optimize_buckets combines multiple buckets in a bucket-family into
  //a single larger bucket, and also does a sort.
  //If optimize_buckets has not been requested, still do the sort.
  if (m_optimize_buckets) m_bucket_repository.optimize_buckets();
  else m_bucket_repository.internal_sort_bucket_entities();

  // ------------------------------

  m_sync_state = SYNCHRONIZED ;

  return true ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

enum { PART_ORD_UNIVERSAL = 0 };
enum { PART_ORD_OWNED     = 1 };
enum { PART_ORD_SHARED    = 2 };

namespace {

void pack_induced_memberships( CommAll & comm ,
                               const std::vector<Entity*> & entity_comm )
{
  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity & entity = **i ;

    if ( in_shared( entity , entity.owner_rank() ) ) {
      // Is shared with owner, send to owner.

      OrdinalVector empty , induced ;

      induced_part_membership( entity , empty , induced );

      CommBuffer & buf = comm.send_buffer( entity.owner_rank() );

      unsigned tmp = induced.size();

      buf.pack<unsigned>( tmp );

      for ( OrdinalVector::iterator
            j = induced.begin() ; j != induced.end() ; ++j ) {
        buf.pack<unsigned>( *j );
      }
    }
  }
}

void generate_send_list( const size_t sync_count ,
                         const unsigned p_rank ,
                         const std::vector<Entity*>    & entity_comm ,
                               std::vector<EntityProc> & send_list )
{
  for ( std::vector<Entity*>::const_iterator
        i = entity_comm.begin() ; i != entity_comm.end() ; ++i ) {

    Entity & entity = **i ;

    if ( entity.owner_rank() == p_rank &&
         entity.synchronized_count() == sync_count ) {

      for ( PairIterEntityComm ec = entity.comm() ; ! ec.empty() ; ++ec ) {
        EntityProc tmp( & entity , ec->proc );
        send_list.push_back( tmp );
      }
    }
  }

  {
    std::sort( send_list.begin() , send_list.end() , EntityLess() );
    std::vector<EntityProc>::iterator i =
      std::unique( send_list.begin() , send_list.end() );
    send_list.erase( i , send_list.end() );
  }
}

void pack_part_memberships( CommAll & comm ,
                            const std::vector<EntityProc> & send_list )
{
  for ( std::vector<EntityProc>::const_iterator
        i = send_list.begin() ; i != send_list.end() ; ++i ) {

    Entity & entity = * i->first ;

    std::pair<const unsigned *, const unsigned *>
      part_ord = entity.bucket().superset_part_ordinals();

    // I am the owner; therefore, the first three members are
    // universal, uses, and owns.  Don't send them.

    // I am the owner.  The first two memberships are
    // universal_part and locally_owned_part.  The third
    // membership may be globally_shared_part ;

    const unsigned count_all  = part_ord.second - part_ord.first ;
    const unsigned count_skip =
      ( 2 < count_all && part_ord.first[2] == PART_ORD_SHARED ) ? 3 : 2 ;

    const unsigned count_send = count_all - count_skip ;

    const unsigned * const start_send = part_ord.first + count_skip ;

    comm.send_buffer( i->second ).pack<EntityKey>( entity.key() )
                                 .pack<unsigned>( count_send )
                                 .pack<unsigned>( start_send , count_send );
  }
}

}

//  Mesh entity membership changes must be synchronized among
//  processes that share mesh entities and propagated to
//  processes that ghost copies of the mesh entities.
//
//  Precondition: correct shared and ghosting lists.
//
//  Part memberships may have been added or removed
//  either explicitly or indirectly via entity relationships
//  being added or removed.

void BulkData::internal_resolve_shared_membership()
{
  Trace_("stk_classic::mesh::BulkData::internal_resolve_shared_membership");

  ThrowRequireMsg(parallel_size() > 1, "Do not call this in serial");

  const MetaData & meta  = m_mesh_meta_data ;
  ParallelMachine p_comm = m_parallel_machine ;
  const unsigned  p_rank = m_parallel_rank ;
  const unsigned  p_size = m_parallel_size ;
  const PartVector & all_parts = meta.get_parts();

  const Part & part_universal = meta.universal_part();
  const Part & part_owned  = meta.locally_owned_part();
  const Part & part_shared = meta.globally_shared_part();

  // Quick verification of part ordinal assumptions

  ThrowRequireMsg(PART_ORD_UNIVERSAL == part_universal.mesh_meta_data_ordinal(),
                  "Universal part ordinal is wrong, expected "
                  << PART_ORD_UNIVERSAL << ", got: "
                  << part_universal.mesh_meta_data_ordinal());

  ThrowRequireMsg(PART_ORD_OWNED == part_owned.mesh_meta_data_ordinal(),
                  "Owned part ordinal is wrong, expected "
                  << PART_ORD_OWNED << ", got: "
                  << part_owned.mesh_meta_data_ordinal());

  ThrowRequireMsg(PART_ORD_SHARED == part_shared.mesh_meta_data_ordinal(),
                  "Shared part ordinal is wrong, expected "
                  << PART_ORD_SHARED << ", got: "
                  << part_shared.mesh_meta_data_ordinal());

  //  Shared entities may have been modified due to relationship changes.
  //  Send just the current induced memberships from the sharing to
  //  the owning processes.
  {
    CommAll comm( p_comm );

    pack_induced_memberships( comm , m_entity_comm );

    comm.allocate_buffers( p_size / 4 );

    pack_induced_memberships( comm , m_entity_comm );

    comm.communicate();

    for ( std::vector<Entity*>::iterator
          i = m_entity_comm.begin() ; i != m_entity_comm.end() ; ++i ) {

      Entity & entity = **i ;

      if ( entity.owner_rank() == p_rank ) {
        // Receiving from all sharing processes

        OrdinalVector empty , induced_parts , current_parts , remove_parts ;

        induced_part_membership( entity , empty , induced_parts );

        for ( PairIterEntityComm
              ec = entity.sharing() ; ! ec.empty() ; ++ec ) {

          CommBuffer & buf = comm.recv_buffer( ec->proc );

          unsigned count = 0 ; buf.unpack<unsigned>( count );
          for ( unsigned j = 0 ; j < count ; ++j ) {
            unsigned part_ord = 0 ; buf.unpack<unsigned>( part_ord );
            insert_ordinal( induced_parts , part_ord );
          }
        }

        // Remove any part that is an induced part but is not
        // in the induced parts list.

        entity.bucket().supersets( current_parts );

        OrdinalVector::const_iterator induced_parts_begin = induced_parts.begin(),
                                      induced_parts_end   = induced_parts.end();

        for ( OrdinalVector::iterator
              p = current_parts.begin() ; p != current_parts.end() ; ++p ) {
          if ( membership_is_induced( *meta.get_parts()[*p] , entity.entity_rank() ) &&
               ! contains_ordinal( induced_parts_begin, induced_parts_end , *p ) ) {
            remove_parts.push_back( *p );
          }
        }

        internal_change_entity_parts( entity, induced_parts, remove_parts );
      }
    }
  }

  //------------------------------
  // The owners have complete knowledge of memberships.
  // Send membership information to sync the shared and ghosted copies.
  // Only need to do this for entities that have actually changed.

  {
    std::vector<EntityProc> send_list ;

    generate_send_list( m_sync_count, p_rank, m_entity_comm, send_list);

    CommAll comm( p_comm );

    pack_part_memberships( comm , send_list );

    comm.allocate_buffers( p_size / 4 );

    pack_part_memberships( comm , send_list );

    comm.communicate();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm.recv_buffer( p );
      while ( buf.remaining() ) {

        PartVector owner_parts , current_parts , remove_parts ;

        EntityKey key ; buf.unpack<EntityKey>( key );
        unsigned count = 0 ; buf.unpack<unsigned>( count );
        for ( unsigned j = 0 ; j < count ; ++j ) {
          unsigned part_ord = 0 ; buf.unpack<unsigned>( part_ord );
          insert( owner_parts , * all_parts[ part_ord ] );
        }

        // Any current part that is not a member of owners_parts
        // must be removed.

        Entity * const entity = find_entity(m_entity_comm, key, true);

        entity->bucket().supersets( current_parts );

        for ( PartVector::iterator
              ip = current_parts.begin() ; ip != current_parts.end() ; ++ip ) {
          Part * const part = *ip ;
          const unsigned part_ord = part->mesh_meta_data_ordinal();
          if ( PART_ORD_UNIVERSAL != part_ord &&
               PART_ORD_OWNED     != part_ord &&
               PART_ORD_SHARED    != part_ord &&
               ! contain( owner_parts , *part ) ) {
            remove_parts.push_back( part );
          }
        }

        internal_change_entity_parts( *entity , owner_parts , remove_parts );
      }
    }
  }
}

} // namespace mesh
} // namespace stk_classic

