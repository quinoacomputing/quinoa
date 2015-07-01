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
#include <util/ParallelIndex.hpp>
#include <util/ParallelReduce.hpp>

#include <mesh/BulkData.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/Comm.hpp>

namespace phdmesh {

namespace {

std::pair< std::vector<EntityProc>::const_iterator ,
           std::vector<EntityProc>::const_iterator >
span( const std::vector<EntityProc> & v , unsigned entity_type )
{
  const unsigned t1 = entity_type ;
  const unsigned t2 = t1 + 1 ;

  std::pair< std::vector<EntityProc>::const_iterator ,
             std::vector<EntityProc>::const_iterator > result ;

  result.first  = lower_bound(v,t1);
  result.second = lower_bound(v,t2);

  return result ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool comm_verify_shared_entity_values(
  const BulkData & M , EntityType t , const FieldBase & f )
{
  const unsigned parallel_size = M.parallel_size();

  const unsigned max_size = f.max_size(t) *
                            NumericEnum<>::size( f.numeric_type_ordinal() );

  const std::pair< std::vector<EntityProc>::const_iterator ,
                   std::vector<EntityProc>::const_iterator >
    shares = span( M.shared_entities() , t );

  std::vector<EntityProc>::const_iterator ic ;

  ParallelMachine comm = M.parallel();

  CommAll sparse ;

  {
    const unsigned zero = 0 ;
    std::vector<unsigned> comm_size( parallel_size , zero );

    for ( ic = shares.first ; ic != shares.second ; ++ic ) {
      comm_size[ ic->second ] += max_size ;
    }

    const unsigned * const p_size = & comm_size[0] ;

    sparse.allocate_buffers( comm, parallel_size / 4 , p_size, p_size );
  }

  std::vector<unsigned char> scratch( max_size );

  for ( ic = shares.first ; ic != shares.second ; ++ic ) {
    Entity & e = * ic->first ;
    Kernel & k = e.kernel();
    const unsigned this_size = field_data_size( f , k );
    CommBuffer & b = sparse.send_buffer( ic->second );

    unsigned char * ptr = (unsigned char *) field_data( f , e );

    b.pack<unsigned char>( ptr , this_size );
    b.skip<unsigned char>( max_size - this_size );
  }

  sparse.communicate();

  unsigned ok = 1 ;

  for ( ic = shares.first ; ic != shares.second ; ++ic ) {
    Entity & e = * ic->first ;
    Kernel & k = e.kernel();
    const unsigned this_size = field_data_size( f , k );
    CommBuffer & b = sparse.recv_buffer( ic->second );

    unsigned char * ptr = (unsigned char *) field_data( f , e );
    unsigned char * scr = & scratch[0] ;

    b.unpack<unsigned char>( scr , this_size );
    b.skip<unsigned char>( max_size - this_size );

    // Compare data and scratch
    for ( unsigned j = 0 ; ok && j < this_size ; ++j ) {
      ok = ptr[j] == scr[j] ;
    }
  }

  all_reduce( comm , Min<1>( & ok ) );

  return ok ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void comm_mesh_discover_sharing( BulkData & M )
{
  static const char method[] = "phdmesh::comm_mesh_discover_sharing" ;

  const MetaData & S = M.mesh_meta_data();
  const unsigned p_rank = M.parallel_rank();

  std::vector< ParallelIndex::key_type > local ;
  std::vector< ParallelIndex::KeyProc > global ;
  std::vector<EntityProc> share ;

  unsigned count = 0 ;
  for ( unsigned k = 0 ; k < EntityTypeEnd ; ++k ) {
    const EntitySet & eset = M.entities( k );
    count += eset.size();
  }

  local.reserve( count );

  for ( unsigned k = 0 ; k < EntityTypeEnd ; ++k ) {
    const EntitySet & eset = M.entities( k );

    const EntitySet::const_iterator e = eset.end();
          EntitySet::const_iterator i ;
    for ( i = eset.begin() ; i != e ; ++i ) { local.push_back( i->key() ); }
  }

  ParallelIndex par_index( M.parallel() , local );

  par_index.query( global );

  for ( std::vector< ParallelIndex::KeyProc >::iterator
        i = global.begin() ; i != global.end() ; ) {

    const entity_key_type key = i->first ;

    EntityProc ep ;
    ep.first = M.get_entity( key , method );

    for ( ; i != global.end() && key == i->first ; ++i ) {
      ep.second = i->second ;
      share.push_back( ep );
    }
  }

  sort_unique( share );

  // Now revise ownership

  Part * const owns_part = & S.locally_owned_part();

  const std::vector<EntityProc>::iterator ipe = share.end();
        std::vector<EntityProc>::iterator ip ;

  for ( ip = share.begin() ; ip != ipe ; ) {
    Entity * const entity = ip->first ;
    const unsigned p_send = ip->second ;

    for ( ; ip != ipe && ip->first == entity ; ++ip );

    std::vector<Part*> add_parts , remove_parts ;

    const unsigned p_owner = p_rank < p_send ? p_rank : p_send ;

    if ( p_owner != p_rank ) {
      remove_parts.push_back( owns_part );
    }
    else {
      add_parts.push_back( owns_part );
    }

    M.change_entity_parts( *entity , add_parts , remove_parts );
    M.change_entity_owner( *entity , p_owner );
  }

  // Set the shares relation
 
  M.set_shares( share );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// If an entity is not owned and not UsedBy an owned entity then remove it
// from the sharing and the processor.  If sharing changes then regenerate
// the aura.

bool comm_mesh_scrub_sharing( BulkData & M )
{
  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  std::vector<EntityProc> shares( M.shared_entities() );

  bool changed = false ;

  std::vector<EntityProc>::iterator i ;

  for ( i = shares.end() ; i != shares.begin() ; ) {
    const std::vector<EntityProc>::iterator ie = i ;

    Entity * const entity = (--i)->first ;
    const EntityType entity_type = entity->entity_type();

    for ( ; i != shares.begin() && entity == i->first ; --i );
    if ( entity != i->first ) { ++i ; }

    bool destroy_it = p_rank != entity->owner_rank();

    for ( PairIterRelation
          con = entity->relations() ; con && destroy_it ; ++con ) {

      destroy_it = ! ( entity_type < con->entity_type() &&
                       p_rank == con->entity()->owner_rank() );
    }

    if ( destroy_it ) {
      std::vector<EntityProc>::iterator ib = i ;
      for ( ; ib != ie ; ++ib ) { ib->first = NULL ; }
      M.destroy_entity( entity );
      changed = true ;
    }
  }

  // Inform sharing processors of the destruction

  CommAll all( M.parallel() );

  for ( i = shares.begin() ; i != shares.end() ; ++i ) {
    all.send_buffer( i->second ).skip<unsigned char>(1);
  }

  const bool symmetric = true ;

  changed = all.allocate_buffers( p_size / 4 , symmetric , changed );

  if ( changed ) {

    for ( i = shares.begin() ; i != shares.end() ; ++i ) {
      const unsigned char flag_remove = NULL == i->first ;
      all.send_buffer( i->second ).pack<unsigned char>( flag_remove );
    }

    all.communicate();

    for ( i = shares.begin() ; i != shares.end() ; ++i ) {
      unsigned char flag_remove ; 
      all.recv_buffer( i->second ).unpack<unsigned char>( flag_remove );
      if ( flag_remove ) { i->first = NULL ; }
    }
  
    for ( i = shares.end() ; i != shares.begin() ; ) {
      if ( NULL == (--i)->first ) { i = shares.erase( i ); }
    }

    M.set_shares( shares );

    comm_mesh_regenerate_aura( M );
  }

  return changed ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void pack_info( CommAll & all , const std::vector<EntityProc> & shares )
{
  const std::vector<EntityProc>::const_iterator i_end = shares.end();
        std::vector<EntityProc>::const_iterator i     = shares.begin();

  while ( i != i_end ) {
    const EntityProc & ep = *i ; ++i ;

    Entity & entity = * ep.first ;
    PairIterRelation relation = entity.relations();
    const entity_key_type key = entity.key();
    const unsigned p_owner  = entity.owner_rank();
    const unsigned numrelation = relation.size();

    const unsigned p_send = ep.second ;

    const std::pair<const unsigned * , const unsigned *> part_ordinals =
      entity.kernel().superset_part_ordinals();

    const unsigned numparts = part_ordinals.second - part_ordinals.first ;

    CommBuffer & b = all.send_buffer( p_send );

    b.pack<entity_key_type>( key );
    b.pack<unsigned>( p_owner );
    b.pack<unsigned>( numparts );
    b.pack<unsigned>( numrelation );

    b.pack<unsigned>( part_ordinals.first , numparts );

    for ( ; relation ; ++relation ) {
      uint64_type con[2] ;
      con[0] = relation->attribute();
      con[1] = relation->entity()->key();
      b.pack<uint64_type>( con , 2 );
    }
  }
}

bool unpack_info_verify(
  CommAll & all ,
  const std::vector<EntityProc> & shares ,
  std::string & error_msg )
{
  static const char method[] =
    "phdmesh::comm_mesh_verify_parallel_consistency" ;

  const unsigned    u_zero = 0 ;
  const uint64_type ul_zero = 0 ;
  bool result = true ;

  const std::vector<EntityProc>::const_iterator i_end = shares.end();
        std::vector<EntityProc>::const_iterator i     = shares.begin();

  std::vector<unsigned> recv_ordinal ;
  std::vector<uint64_type> recv_relation ;

  while ( result && i != i_end ) {
    const EntityProc & ep = *i ; ++i ;

    Entity & entity = * ep.first ;
    Kernel & kernel = entity.kernel();
    BulkData   & mesh   = kernel.mesh();
    const MetaData & mesh_meta_data = mesh.mesh_meta_data();
    const PartSet & mesh_parts = mesh_meta_data.get_parts();
    const PairIterRelation relation = entity.relations();
    Part * const owns_part = & mesh_meta_data.locally_owned_part();

    const unsigned relation_size = relation.size();
    const entity_key_type key = entity.key();
    const unsigned p_owner  = entity.owner_rank();
    const unsigned p_local  = mesh.parallel_rank();
    const unsigned p_recv   = ep.second ;

    const std::pair<const unsigned * , const unsigned * >
      part_ordinals = kernel.superset_part_ordinals();

    const unsigned parts_size = part_ordinals.second - part_ordinals.first ;

    CommBuffer & b = all.recv_buffer( p_recv );

    entity_key_type recv_key ;   b.unpack<entity_key_type>( recv_key );
    unsigned recv_p_owner ;       b.unpack<unsigned>( recv_p_owner );
    unsigned recv_parts_size ;    b.unpack<unsigned>( recv_parts_size );
    unsigned recv_relation_size ; b.unpack<unsigned>( recv_relation_size );

    recv_ordinal.assign( recv_parts_size , u_zero );
    {
      unsigned * const tmp = & recv_ordinal[0] ;
      b.unpack<unsigned>( tmp , recv_parts_size );
    }

    recv_relation.assign( recv_relation_size * 2 , ul_zero );
    {
      uint64_type * const tmp = & recv_relation[0] ;
      b.unpack<uint64_type>( tmp , recv_relation_size * 2 );
    }

    result = recv_key == key && recv_p_owner == p_owner ;

    if ( result ) {

      unsigned j_this = 0 ;
      unsigned j_recv = 0 ;

      while ( result &&
              j_this < parts_size &&
              j_recv < recv_parts_size ) {

        const int ord_recv = recv_ordinal[ j_recv ];
        if ( ord_recv < 0 || ((int) mesh_parts.size()) <= ord_recv ) {
          // Bad ordinal
          result = false ;
        }
        else {
          Part * const part_this = mesh_parts[ part_ordinals.first[ j_this ] ] ;
          Part * const part_recv = mesh_parts[ ord_recv ];

          if ( owns_part == part_this ) {
            // Local ownership by this processor
            result = p_local == p_owner &&
                     p_local == recv_p_owner ;
            ++j_this ;
          }
          else if ( owns_part == part_recv ) {
            // Remote ownership by sending processor
            result = p_recv == p_owner &&
                     p_recv == recv_p_owner ;
            ++j_recv ;
          }
          else if ( part_this != part_recv ) {
            result = false ;
          }
          else {
            ++j_this ;
            ++j_recv ;
          }
        }
      }
    }

    if ( relation_size != recv_relation_size ) { result = false ; }

    PairIterRelation::iterator ic ;

    ic = relation.begin();

    for ( unsigned k = 0 ; result && k < recv_relation_size ; ++k , ++ic ) {
      const Relation & con = *ic ;
      const unsigned k2 = k * 2 ;
      const relation_attr_type recv_rel_attr = recv_relation[k2] ;
      const entity_key_type    recv_rel_key  = recv_relation[k2+1] ;
      if ( recv_rel_attr != con.attribute() ||
           recv_rel_key  != con.entity()->key() ) { result = false ; }
    }

    if ( ! result ) {
      std::ostringstream os ;

      os << std::endl ;
      os << method << " ERROR :" << std::endl ;
      os << "  P" << p_local << " : " ;
      print_entity_key( os , key );
      os << ".{ Owner = P" << p_owner << " ; Parts =" ;
      for ( const unsigned * j = part_ordinals.first ;
                             j < part_ordinals.second ; ++j ) {
        os << " " << mesh_parts[ *j ]->name();
      }
      os << " ;" << std::endl << "  Relationships =" ;
      for ( ic = relation.begin() ; ic != relation.end() ; ++ic ) {
        Entity & con_e = * ic->entity();

        const bool has_owned = con_e.kernel().has_superset( *owns_part );

        os << std::endl << "    " ;
        os << *ic ;
        os << ".{ Owner = P" << con_e.owner_rank();
        if ( has_owned ) { os << " owned" ; }
        os << " }" ;
      }
      os << " } != " << std::endl ;
      os << "  P" << p_recv << " : " ;
      print_entity_key( os , recv_key );
      os << ".{ Owner = P" << recv_p_owner << " ; Parts =" ;
      for ( unsigned j = 0 ; j < recv_parts_size ; ++j ) {
        const unsigned ord_recv = recv_ordinal[j] ;
        if ( mesh_parts.size() <= ord_recv ) {
          os << " unsigned[" << ord_recv << "]" ;
        }
        else {
          os << " " << mesh_parts[ord_recv]->name();
        }
      }
      os << " ;" << std::endl << "  Relationships =" ;
      for ( unsigned j = 0 ; j < recv_relation_size * 2 ; ) {
        const relation_attr_type attr    = recv_relation[j] ; ++j ;
        const entity_key_type    rel_key = recv_relation[j] ; ++j ;
        os << std::endl << "    " ;
        print_relation( os , attr , rel_key );
      }
      os << " }" << std::endl ;

      error_msg = os.str();
    }
  }

  return result ;
}

bool verify_parallel_attributes(
  BulkData & M , unsigned type , std::string & msg )
{
  bool result = true ;

  const MetaData & S = M.mesh_meta_data();
  const unsigned p_rank = M.parallel_rank();
  Part & uses_part = S.locally_used_part();
  Part & owns_part = S.locally_owned_part();

  const EntitySet & es = M.entities( type );

  std::pair< std::vector<EntityProc>::const_iterator ,
             std::vector<EntityProc>::const_iterator >
    aura_span = span( M.ghost_destination() , type );

  // Iterate everything in identifier ordering.

  const EntitySet::iterator i_end = es.end();
        EntitySet::iterator i     = es.begin();

  while ( result && i != i_end ) {
    Entity & entity = *i ; ++i ;
    Kernel & kernel = entity.kernel();

    const bool uses = kernel.has_superset( uses_part );
    const bool owns = kernel.has_superset( owns_part );
    const unsigned p_owner = entity.owner_rank();

    const bool shares = entity.sharing();

    // Shared is a subset of uses.

    if ( shares && ! uses ) { result = false ; }

    if ( owns ) {
      // Owner is local
      if ( p_owner != p_rank ) { result = false ; }
    }
    else {
      // Owner is remote: if uses then must be shared
      if ( p_owner == p_rank ) { result = false ; }
      if ( uses && ! shares ) { result = false ; }
    }

    bool in_aura = false ;
    unsigned aura_rank = p_owner ;

    // Can appear at most once in the aura parallel relation:

    if ( aura_span.first != aura_span.second &&
         aura_span.first->first == & entity ) {

      in_aura = true ;

      aura_rank = aura_span.first->second ;

      if ( p_owner != aura_rank ) { result = false ; }

      ++aura_span.first ;
    }

    if ( in_aura && uses ) { result = false ; }

    if ( ! result ) {
      std::ostringstream os ;
      os << "P" << p_rank << " : Error with parallel attributes "
         << std::endl ;
      os << "  " ;
      print_entity_key( os , entity.key() );
      os << ".owner(P" << p_owner << ")" ;
      if ( uses )   { os << " , uses" ; }
      if ( owns )   { os << " , owns" ; }
      if ( shares ) { os << " , shares" ; }
      if ( in_aura ) { os << " , in_aura(P" << aura_rank << ")" ; }
      else           { os << " , not_in_aura" ; }
      os << std::endl ;
      msg = os.str();
    }
  }
  return result ;
}

}

bool comm_mesh_verify_parallel_consistency( BulkData & M )
{
  // Exchange the shared entities' identifiers, part mesh ordinals, and
  // owner processor.  Should be fully consistent modulo the owner part.

  const unsigned p_size = M.parallel_size();

  std::string msg ;
  int result = 1 ;

  for ( unsigned i = 0 ; result && i < EntityTypeEnd ; ++i ) {
    result = verify_parallel_attributes( M , i , msg );
  }

  {
    // Verify all shared entities

    const std::vector<EntityProc> & shares = M.shared_entities();

    CommAll all_info( M.parallel() );

    pack_info( all_info , shares );

    all_info.allocate_buffers( p_size / 4 , false );

    pack_info( all_info , shares );

    all_info.communicate();

    if ( ! unpack_info_verify( all_info , shares , msg ) ) { result = false ; }
  }

  // Verify consistency of aura

  if ( ! comm_verify( M.parallel(), M.ghost_source(), M.ghost_destination(), msg) ) {
    result = false ;
  }

  // Global reduction of result flag:

  all_reduce( M.parallel() , Min<1>( & result ) );

  //--------------------------------------------------------------------
  // If an error occured collect the messages on
  // processor zero and output to standard error stream.

  if ( ! result ) { all_write_string( M.parallel() , std::cerr , msg ); }

  return result ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

class SharingComm : public EntityComm {
private:
  SharingComm( const SharingComm & );
  SharingComm & operator = ( const SharingComm & );
public:
  SharingComm() {}
  ~SharingComm() {}

  const char * name() const ;

  void receive_entity(
    CommBuffer & buffer ,
    BulkData & receive_mesh ,
    const unsigned send_source ,
    std::vector<EntityProc> & receive_info ) const ;
};

const char * SharingComm::name() const
{ static const char n[] = "phdmesh::SharingComm" ; return n ; }

void SharingComm::receive_entity(
  CommBuffer & buffer ,
  BulkData & receive_mesh ,
  const unsigned send_source ,
  std::vector<EntityProc> & receive_info ) const
{
  // receive_info is the new_aura_range

  static const char method[] = "phdmesh::SharingComm::receive_entity" ;

  const unsigned local_rank = receive_mesh.parallel_rank();
  const MetaData & S = receive_mesh.mesh_meta_data();
  Part * const uses_part = & S.locally_used_part();
  Part * const owns_part = & S.locally_owned_part();

  entity_key_type       key ;
  unsigned              owner_rank ;
  std::vector<Part*>    add ;
  std::vector<Relation> relations ;
  std::vector<unsigned> send_dest ;

  if ( send_source == local_rank ) {
    std::string msg( method );
    msg.append(" FAILED, Received from local-processor" );
    throw std::logic_error( msg );
  }

  unpack_entity( buffer , receive_mesh ,
                 key , owner_rank ,
                 add , relations , send_dest );

  // Must have been sent by the owner.

  if ( send_source != owner_rank ) {
    std::string msg( method );
    msg.append( " FAILED, Received from non-owner processor" );
    throw std::logic_error( msg );
  }

  // Remove owns part

  {
    std::vector<Part*>::iterator ip = add.begin() ;
    for ( ; ip != add.end() && owns_part != *ip ; ++ip );
    if ( ip != add.end() ) { add.erase( ip ); }
  }

  EntityProc ep ;

  ep.first  = receive_mesh.get_entity( key );
  ep.second = send_source ;

  // If entity exists it must be currently an aura, i.e. not a uses_part

  if ( NULL != ep.first ) {
    std::vector<EntityProc>::iterator k = lower_bound( receive_info , *ep.first );

    if ( ep.first->kernel().has_superset( *uses_part ) ||
         k == receive_info.end() || k->first != ep.first ) {
      std::string msg( method );
      msg.append( "FAILED, entity exists but is not an aura" );
      throw std::logic_error( msg );
    }
    receive_info.erase( k );
  }

  ep.first = & receive_mesh.declare_entity( key , add , owner_rank );

  receive_mesh.declare_relation( *ep.first , relations );
}

bool not_member( const std::vector<EntityProc> & s , const EntityProc & m )
{
  const std::vector<EntityProc>::const_iterator i = lower_bound( s , m );
  return i == s.end() || m != *i ;
}

// Owner packs sharing

void pack_sharing( CommAll & all , const std::vector<EntityProc> & sharing )
{
  const unsigned p_rank = all.parallel_rank();

  std::vector<EntityProc>::const_iterator i , j ;

  for ( i = sharing.begin() ; i != sharing.end() ; ) {
    const std::vector<EntityProc>::const_iterator ib = i ;
    for ( ; i != sharing.end() && i->first == ib->first ; ++i );
    const std::vector<EntityProc>::const_iterator ie = i ;

    if ( p_rank == ib->first->owner_rank() ) {
      const unsigned n = std::distance( ib , ie ) - 1 ;
      const entity_key_type key = ib->first->key();
      for ( i = ib ; i != ie ; ++i ) {
        CommBuffer & buf = all.send_buffer( i->second );
        buf.pack<entity_key_type>( key );
        buf.pack<unsigned>( n );
        if ( n ) {
          for ( j = ib ; j != ie ; ++j ) {
            if ( i != j ) { buf.pack<unsigned>( j->second ); }
          }
        }
      }
    }
  }
}

void unpack_sharing( CommAll & all , BulkData & M ,
                     std::vector<EntityProc> & sharing ,
                     const char * method )
{
  const unsigned p_size = all.parallel_size();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    while ( buf.remaining() ) {
      entity_key_type key ; buf.unpack<entity_key_type>( key );
      unsigned n ;        buf.unpack<unsigned>( n );
      EntityProc tmp ;
      tmp.first = M.get_entity( key , method );
      tmp.second = p ; // The owner sent it
      sharing.push_back( tmp );
      for ( unsigned i = 0 ; i < n ; ++i ) {
        buf.unpack<unsigned>( tmp.second );
        sharing.push_back( tmp );
      }
    }
  }
}


}

void comm_mesh_add_sharing( BulkData & M , const std::vector<EntityProc> & add )
{
  static const char method[] = "phdmesh::comm_mesh_add_sharing" ;

  const SharingComm mgr ;
  const unsigned p_rank = M.parallel_rank();
  const unsigned p_size = M.parallel_size();

  const std::vector<EntityProc> & old_shares = M.shared_entities();

  std::vector<EntityProc> new_shares ;

  //--------------------------------------------------------------------
  // Have the owners send the to-be-shared entities to the
  // specified processors.  If an aura relationship exists
  // then remove it and replace it with a shared relationship.
  {
    std::vector<EntityProc> new_aura_domain( M.ghost_source() );
    std::vector<EntityProc> new_aura_range(  M.ghost_destination() );

    {
      CommAll all( M.parallel() );

      for ( std::vector<EntityProc>::const_iterator
            j = add.begin() ; j != add.end() ; ++j ) {

        const unsigned p_owner = j->first->owner_rank();

        if ( p_owner != j->second ) {
          if ( p_owner != p_rank ) {
            all.send_buffer( p_owner ).skip<entity_key_type>(2);
          }
          else if ( not_member( old_shares , *j ) ) { 
            new_shares.push_back( *j );
          }
        }
      }

      all.allocate_buffers( p_size / 4 , false /* not symmetric */ );

      for ( std::vector<EntityProc>::const_iterator
            j = add.begin() ; j != add.end() ; ++j ) {

        const unsigned p_owner = j->first->owner_rank();

        if ( p_owner != j->second && p_owner != p_rank ) {
          entity_key_type data[2] ;
          data[0] = j->first->key();
          data[1] = j->second ;
          all.send_buffer( p_owner ).pack<entity_key_type>(data,2);
        }
      }

      all.communicate();

      for ( unsigned j = 0 ; j < p_size ; ++j ) {
        CommBuffer & buf = all.recv_buffer(j);
        while ( buf.remaining() ) {
          entity_key_type data[2] ;
          buf.unpack<entity_key_type>( data , 2 );
          EntityProc tmp ;
          tmp.first = M.get_entity( data[0] , method );
          tmp.second = data[1] ;
          if ( not_member( old_shares , tmp ) ) {
            new_shares.push_back( tmp );
          }
        }
      }
    }

    sort_unique( new_shares );

    // Change the newly shared mesh entities on the owner processor
    // Remove from the aura_domain.
    for ( std::vector<EntityProc>::iterator
          i = new_shares.begin() ; i != new_shares.end() ; ++i ) {

      const std::vector<EntityProc>::iterator k = lower_bound( new_aura_domain , *i );

      if ( k != new_aura_domain.end() && *k == *i ) {
        new_aura_domain.erase( k );
      }
    }

    // Communicate new-shared mesh entities from owner to sharing processors.
    // Remove from the aura_range, if present.

    communicate_entities( mgr, M, M, new_shares, new_aura_range, false );

    // Some aura entities have become shared entities.

    M.set_ghosting( new_aura_domain , new_aura_range );
  }
  // Shared mesh entities are now in place.
  //--------------------------------------------------------------------
  // Owners inform the sharing processors of the extent of the sharing.

  new_shares.insert( new_shares.end() , old_shares.begin() ,
                                        old_shares.end() );

  sort_unique( new_shares );

  // Now owners send to sharing the extent of sharing
  {
    CommAll all( M.parallel() );

    pack_sharing( all , new_shares );

    all.allocate_buffers( p_size / 4 , false );

    pack_sharing( all , new_shares );

    all.communicate();

    unpack_sharing( all , M , new_shares , method );
  }

  sort_unique( new_shares );

  M.set_shares( new_shares );

  // Regenerate aura

  comm_mesh_regenerate_aura( M );
}


}

