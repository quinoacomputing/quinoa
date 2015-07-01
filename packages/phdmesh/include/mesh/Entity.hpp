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

#ifndef phdmesh_Entity_hpp
#define phdmesh_Entity_hpp

#include <limits>
#include <iosfwd>

#include <util/PairIter.hpp>
#include <mesh/Types.hpp>
#include <mesh/Kernel.hpp>
#include <mesh/Relation.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
// Functions for encoding / decoding entity keys.

EntityType      entity_type( entity_key_type );
entity_id_type  entity_id(   entity_key_type );
entity_key_type entity_key( EntityType , entity_id_type );

//----------------------------------------------------------------------

/** Span of a sorted relations for a given domain entity.
 *  Members are sorted by
 *  (1) range entity type,
 *  (2) relation type,
 *  (3) identifier, and
 *  (4) range entity identifier.
 */
typedef PairIter< std::vector<Relation>::const_iterator > PairIterRelation ;

//----------------------------------------------------------------------
/** A mesh entity has an entity type, identifier, relations, and
 *  resides within a mesh kernel.  The mesh kernel holds its field data.
 */
class Entity : public SetvMember< entity_key_type > {
private:

  std::vector<Relation> m_relation ;   // Relationships
  KernelSet::iterator   m_kernel ;     // Containing kernel
  unsigned              m_kernel_ord ; // Ordinal in the kernel
  unsigned              m_owner_rank ; // Parallel owner rank
  PairIterEntityProc        m_sharing ;

public:

  EntityType entity_type() const ;

  entity_id_type identifier() const ;

  /** Kernel in which this mesh entity resides */
  Kernel & kernel() const { return * m_kernel ; }

  /** Kernel in which this mesh entity resides */
  unsigned kernel_ordinal() const { return m_kernel_ord ; }

  //------------------------------------
  /** All relations */
  PairIterRelation relations() const { return PairIterRelation( m_relation ); }

  /** Relations with entities of a given entity type */
  PairIterRelation relations( EntityType type , unsigned kind = 0 ) const ;

  //------------------------------------
  /** Owning parallel processor rank */
  unsigned owner_rank() const { return m_owner_rank ; }

  /** Sharing processor information */
  const PairIterEntityProc & sharing() const { return m_sharing ; }

  //------------------------------------

  Entity();
  ~Entity();

private:

  Entity( const Entity & );
  Entity & operator = ( const Entity & );

  friend class BulkData ;
};

typedef Setv<Entity> EntitySet ;

//----------------------------------------------------------------------

/** Print identifier and relations */
std::ostream &
print_entity( std::ostream & , const std::string & lead , const Entity & );

std::ostream &
print_entity_key( std::ostream & , EntityType type , entity_id_type id );

std::ostream &
print_entity_key( std::ostream & , entity_key_type key );

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

enum {
  entity_key_type_shift  =
    std::numeric_limits<entity_key_type>::digits - entity_key_type_digits ,
  entity_id_digits =
     entity_key_type_shift < std::numeric_limits<entity_id_type>::digits ?
     entity_key_type_shift : std::numeric_limits<entity_id_type>::digits
};

enum { entity_id_end = ((entity_key_type) 1) << entity_id_digits };

inline
entity_key_type entity_key( EntityType type , entity_id_type id )
{
  enum { mask = ~((entity_id_type) 0 ) >>
           ( std::numeric_limits<entity_id_type>::digits - entity_id_digits ) };
  return ( ((entity_key_type) type) << entity_key_type_shift ) | ( id & mask );
}

inline
EntityType entity_type( entity_key_type key )
{ return EntityType( key >> entity_key_type_shift ); }

inline
entity_id_type entity_id( entity_key_type key )
{
  enum { mask = ~((entity_id_type) 0 ) >>
           ( std::numeric_limits<entity_id_type>::digits - entity_id_digits ) };
  return ((entity_id_type) key ) & mask ;
}

inline
EntityType Entity::entity_type() const
{ return phdmesh::entity_type( key() ); }

inline
entity_id_type Entity::identifier() const
{ return phdmesh::entity_id( key() ); }

} // namespace phdmesh

#endif

