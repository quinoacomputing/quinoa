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

#ifndef phdmesh_Relation_hpp
#define phdmesh_Relation_hpp

#include <iosfwd>
#include <limits>

#include <mesh/Types.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** Relationship to mesh entities with attributes.
 *
 *  Design is for the DomainEntity to own a subset of the relation,
 *    DomainEntity -> { ( RelationAttribute , RangeEntity ) }
 */
class Relation {
private:
  relation_attr_type m_attr ;
  Entity           * m_entity ;
public:

  ~Relation() {}

  Relation() : m_attr(0), m_entity(NULL) {}

  Relation( const Relation & r )
    : m_attr( r.m_attr ), m_entity(r.m_entity) {}

  Relation & operator = ( const Relation & r )
    { m_attr = r.m_attr ; m_entity = r.m_entity ; return *this ; }

  Relation( relation_attr_type arg_attr , Entity & arg_entity );

  Relation( Entity & arg_entity ,
            unsigned arg_identifier ,
            unsigned arg_kind = 0 ,
            bool     arg_converse = false );

  relation_attr_type attribute() const { return m_attr ; }

  unsigned   kind() const ;
  EntityType entity_type() const ;
  unsigned   identifier() const ;
  bool       forward() const ;
  bool       converse() const ;

  Entity * entity() const { return m_entity ; }

  bool operator == ( const Relation & r ) const
    { return m_attr == r.m_attr && m_entity == r.m_entity ; }

  bool operator != ( const Relation & r ) const
    { return m_attr != r.m_attr || m_entity != r.m_entity ; }

  bool operator < ( const Relation & r ) const ;
};

//----------------------------------------------------------------------
/** Relation attribute type is an encoding of the
 *  a) kind of relationship,
 *  b) entity type,
 *  c) direction of the relationship ( forward or converse ), and then
 *  d) local relationship identifier.
 *
 *  Attribute values are ordered by relation kind, entity type,
 *  direction (forward then converse), and then identifier.
 */ 
relation_attr_type
relation_attr( EntityType entity_type ,
               unsigned   identifier ,
               unsigned   relation_kind = 0 ,
               bool       converse = false );

unsigned   relation_kind(        relation_attr_type );
EntityType relation_entity_type( relation_attr_type );
unsigned   relation_identifier(  relation_attr_type );
bool       relation_forward(     relation_attr_type );
bool       relation_converse(    relation_attr_type );

//----------------------------------------------------------------------

std::ostream &
print_relation( std::ostream & , relation_attr_type );

std::ostream &
print_relation( std::ostream & , relation_attr_type , entity_key_type );

std::ostream & operator << ( std::ostream & , const Relation & );

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

enum {
  relation_attr_kind_digits = 4 ,

  relation_attr_kind_shift =
    std::numeric_limits<relation_attr_type>::digits -
      relation_attr_kind_digits ,

  relation_attr_entity_type_shift =
    relation_attr_kind_shift - entity_key_type_digits ,

  relation_attr_converse_shift =
    relation_attr_entity_type_shift - 1 ,

  relation_attr_identifier_digits =
    relation_attr_converse_shift < std::numeric_limits<unsigned>::digits ?
    relation_attr_converse_shift : std::numeric_limits<unsigned>::digits
};

enum {
  relation_identifer_end =
    ((relation_attr_type) 1) << relation_attr_identifier_digits
};

inline
EntityType relation_entity_type( relation_attr_type attr )
{
  enum { relation_attr_type_mask =
    (~(0u)) >> ( std::numeric_limits<unsigned>::digits -
                 entity_key_type_digits ) };

  return EntityType(
    relation_attr_type_mask & ( attr >> relation_attr_entity_type_shift ) );
}

inline
unsigned relation_identifier( relation_attr_type attr )
{
  enum { relation_attr_identifier_mask =
    (~(0u)) >> ( std::numeric_limits<unsigned>::digits -
                 relation_attr_identifier_digits ) };

  return relation_attr_identifier_mask & attr ;
}

inline
unsigned relation_kind( relation_attr_type attr )
{
  enum { relation_attr_kind_mask =
    (~(0u)) >> ( std::numeric_limits<unsigned>::digits -
                 relation_attr_kind_digits ) };

  return relation_attr_kind_mask & ( attr >> relation_attr_kind_shift );
}

inline
bool relation_converse( relation_attr_type attr )
{ return ( attr >> relation_attr_converse_shift ) & 01 ; }

inline
bool relation_forward( relation_attr_type attr )
{ return ! relation_converse( attr ); }

inline
unsigned Relation::identifier() const
  { return relation_identifier( m_attr ); }

inline
bool Relation::forward() const
  { return relation_forward( m_attr ); }

inline
bool Relation::converse() const
  { return relation_converse( m_attr ); }

inline
unsigned Relation::kind() const
  { return relation_kind( m_attr ); }

inline
EntityType Relation::entity_type() const
  { return relation_entity_type( m_attr ); }

}

#endif

