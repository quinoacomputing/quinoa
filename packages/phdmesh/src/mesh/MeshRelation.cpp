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

#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/Entity.hpp>
#include <mesh/Kernel.hpp>
#include <mesh/Relation.hpp>
#include <mesh/FieldData.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

std::ostream &
print_relation( std::ostream & s , relation_attr_type attr )
{
  s << entity_type_name( relation_entity_type( attr ) );
  s << "[" ;
  s << relation_kind( attr );
  s << "." ;
  s << relation_identifier( attr );
  s << "]" ;
  {
    const char fwd[] = "->" ;
    const char con[] = "<-" ;
    s << ( relation_converse( attr ) ? con : fwd );
  }
  return s ;
}

std::ostream &
print_relation( std::ostream & s , relation_attr_type attr ,
                                   entity_key_type key )
{
  print_relation( s , attr );
  print_entity_key( s , key );

  return s ;
}

std::ostream &
operator << ( std::ostream & s , const Relation & con )
{
  Entity * const e = con.entity();

  print_relation( s , con.attribute() );

  if ( e ) { print_entity_key( s , e->key() ); }
  else     { s << "?" ; }

  return s ;
}

//----------------------------------------------------------------------

relation_attr_type
relation_attr( EntityType arg_entity_type ,
               unsigned   arg_identifier ,
               unsigned   arg_kind ,
               bool       arg_converse )
{
  // Substract one from the id_end to quiet the GNU-4.2 compiler warning
  // when using 64bit relation attribute and 32bit identifier

  enum {
    id_end   = (((relation_attr_type)1) << relation_attr_identifier_digits)-1,
    kind_end = ( (relation_attr_type)1) << relation_attr_kind_digits
  };

  const bool bad_type = (unsigned) EntityTypeEnd <= (unsigned) arg_entity_type ;
  const bool bad_id   = id_end   <= (relation_attr_type) arg_identifier ;
  const bool bad_kind = kind_end <= (relation_attr_type) arg_kind ;

  if ( bad_type || bad_id || bad_kind ) {
    std::ostringstream msg ;
    msg << "phdmesh::relation_attr( " ;
    if ( bad_type ) { msg << "BAD ARGUMENT = " ; }
    msg << entity_type_name( arg_entity_type );
    msg << " , " ;
    if ( bad_id ) { msg << "BAD ARGUMENT = " ; }
    msg << arg_identifier ;
    msg << " , " ;
    if ( bad_kind ) { msg << "BAD ARGUMENT = " ; }
    msg << arg_kind << " , " ;
    if ( arg_converse ) { msg << "true" ; }
    else                { msg << "false" ; }
    msg << " )" ;
    throw std::invalid_argument( msg.str() );
  }

  const relation_attr_type
    attr_converse   = ( arg_converse ? 1 : 0 ) ,
    attr_kind       = arg_kind ,
    attr_rank       = arg_entity_type ,
    attr_identifier = arg_identifier ;

  const relation_attr_type attr =
    ( attr_converse << relation_attr_converse_shift ) |
    ( attr_kind     << relation_attr_kind_shift ) |
    ( attr_rank     << relation_attr_entity_type_shift ) |
    ( attr_identifier );

  return attr ;
}

Relation::Relation( Entity & arg_entity ,
                    unsigned arg_identifier ,
                    unsigned arg_kind ,
                    bool     arg_converse )
  : m_attr( relation_attr( arg_entity.entity_type() ,
                           arg_identifier , arg_kind , arg_converse ) ),
    m_entity( & arg_entity )
{}

Relation::Relation( relation_attr_type arg_attr , Entity & arg_entity )
  : m_attr( arg_attr ), m_entity( & arg_entity )
{
  if ( relation_entity_type( arg_attr ) != arg_entity.entity_type() ) {
    std::ostringstream msg ;
    msg << "phdmesh::Relation::Relation( "  ;
    print_relation( msg , arg_attr );
    msg << " , " ;
    print_entity_key( msg , arg_entity.key() );
    msg << " ) INCOMPATIBLE ARGUMENTS" ;
    throw std::invalid_argument( msg.str() );
  }
}

bool Relation::operator < ( const Relation & r ) const
{
  bool result ;

  if ( m_attr != r.m_attr ) {
    result = m_attr < r.m_attr ;
  }
  else {
    const entity_key_type lhs = m_entity   ? m_entity->key()   : 0 ;
    const entity_key_type rhs = r.m_entity ? r.m_entity->key() : 0 ;
    result = lhs < rhs ;
  }
  return result ;
}

//----------------------------------------------------------------------

namespace {

struct LessRelation {

  inline
  bool operator()( const Relation & lhs , const Relation & rhs ) const
    { return lhs.operator < ( rhs ); }

  inline
  bool operator()( const Relation & lhs , const relation_attr_type rhs ) const
    { return lhs.attribute() < rhs ; }
};

struct LessEntityPointer {
  inline
  bool operator()( const Entity * const lhs , const Entity * const rhs ) const
    {
      const entity_key_type lhs_key = lhs ? lhs->key() : 0 ;
      const entity_key_type rhs_key = rhs ? rhs->key() : 0 ;
      return lhs_key < rhs_key ;
    }
};

}

PairIterRelation con_span( const std::vector<Relation> & con ,
                       const relation_attr_type lo_attr ,
                       const relation_attr_type hi_attr )
{
  std::vector<Relation>::const_iterator i = con.begin();
  std::vector<Relation>::const_iterator e = con.end();

  i = std::lower_bound( i , e , lo_attr , LessRelation() );
  e = std::lower_bound( i , e , hi_attr , LessRelation() );

  return PairIterRelation( i , e );
}

PairIterRelation Entity::relations( EntityType et , unsigned kind ) const
{
  const EntityType et_next = EntityType( et + 1 );
  const relation_attr_type
    lo_attr = relation_attr( et ,      0 , kind ),
    hi_attr = relation_attr( et_next , 0 , kind );

  std::vector<Relation>::const_iterator i = m_relation.begin();
  std::vector<Relation>::const_iterator e = m_relation.end();

  i = std::lower_bound( i , e , lo_attr , LessRelation() );
  e = std::lower_bound( i , e , hi_attr , LessRelation() );

  return PairIterRelation( i , e );
}

void closure( const Entity & e_from , std::vector<Entity*> & eset )
{
  PairIterRelation rel = e_from.relations();
  for ( ; rel ; ++rel ) {
    if ( rel->forward() ) {
      Entity * const e = rel->entity();
      std::vector<Entity*>::iterator i = eset.begin();
      std::vector<Entity*>::iterator j = eset.end();
      i = std::lower_bound( i , j , e , LessEntityPointer() );
      if ( i == j || e != *i ) {
        eset.push_back( e );
        closure( *e , eset );
      }
    }
  }
}

bool in_closure( const Entity & e_from , const Entity & e )
{
  PairIterRelation rel = e_from.relations();

  bool not_in_closure = & e_from != & e ;

  for ( ; not_in_closure && rel ; ++rel ) {
    if ( rel->forward() ) {
      not_in_closure = ! in_closure( * rel->entity() , e );
    }
  }

  return ! not_in_closure ;
}


//----------------------------------------------------------------------

namespace {

// What are this entity's part memberships that can be deduced from
// this entity's relationship.

void deduce_part_relations( const Entity & e_from ,
                            const Entity & e_to ,
                            const unsigned ident ,
                            const unsigned kind ,
                            PartSet & to_parts )
{
  const EntityType t_to   = e_to.entity_type();
  const EntityType t_from = e_from.entity_type();
  const Kernel   & k_to   = e_to.kernel();
  const Kernel   & k_from = e_from.kernel();
  const BulkData     & mesh   = k_to.mesh();
  const MetaData   & mesh_meta_data = mesh.mesh_meta_data();

  const unsigned univ_part_ord = mesh_meta_data.universal_part().mesh_meta_data_ordinal();
  const unsigned uses_part_ord = mesh_meta_data.locally_used_part().mesh_meta_data_ordinal();
  const unsigned owns_part_ord = mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal();

  const std::vector<PartRelation> & part_rel = mesh_meta_data.get_part_relations();

  const std::pair<const unsigned *, const unsigned *>
    kernel_superset_ordinals = k_from.superset_part_ordinals();

  if ( to_parts.empty() ) {
    to_parts.reserve( kernel_superset_ordinals.second -
                      kernel_superset_ordinals.first );

    for ( const unsigned * i = kernel_superset_ordinals.first ;
                           i < kernel_superset_ordinals.second ; ++i ) {
      if ( univ_part_ord != *i && uses_part_ord != *i && owns_part_ord != *i ){
        Part * const p = & mesh_meta_data.get_part( *i );
        to_parts.push_back( p );
      }
    }
  }
  else {

    for ( const unsigned * i = kernel_superset_ordinals.first ;
                           i < kernel_superset_ordinals.second ; ++i ) {
      if ( univ_part_ord != *i && uses_part_ord != *i && owns_part_ord != *i ){
        insert( to_parts , mesh_meta_data.get_part( *i ) );
      }
    }
  }

  // Now the part-relation based contributions:

  for ( std::vector<PartRelation>::const_iterator
        j = part_rel.begin() ; j != part_rel.end() ; ++j ) {

    const PartRelation & stencil = *j ;

    for ( const unsigned * i = kernel_superset_ordinals.first ;
                           i < kernel_superset_ordinals.second ; ++i ) {
      if ( *i == stencil.m_root->mesh_meta_data_ordinal() &&
            0 <= (*stencil.m_function)( t_from , t_to , ident , kind ) ) {
          insert( to_parts , * stencil.m_target );
      }
    }
  }
}

//----------------------------------------------------------------------

void deduce_part_relations( const Entity & e_to , PartSet & to_parts )
{
  PairIterRelation rel = e_to.relations();

  for ( ; rel ; ++rel ) if ( rel->converse() ) {
    deduce_part_relations( * rel->entity() , e_to ,
                             rel->identifier() , rel->kind() , to_parts );
  }
}

//----------------------------------------------------------------------

void set_field_relations( Entity & e_from ,
                          Entity & e_to ,
                          const unsigned ident ,
                          const unsigned kind )
{
  const std::vector<FieldRelation> & field_rels =
    e_from.kernel().mesh().mesh_meta_data().get_field_relations();

  for ( std::vector<FieldRelation>::const_iterator 
        j = field_rels.begin() ; j != field_rels.end() ; ++j ) {

    const FieldRelation & fr = *j ;

    void ** const ptr = (void**) field_data( * fr.m_root , e_from );

    if ( ptr ) {

      void * const src = field_data( * fr.m_target , e_to );

      const int number =
        field_data_size(*fr.m_root,e_from) / sizeof(void*);

      const int offset =
         (*fr.m_function)( e_from.entity_type() , 
                           e_to.entity_type() , ident , kind );

      if ( 0 <= offset && offset < number ) {
        ptr[ offset ] = src ;
      }
    }
  }
}

//----------------------------------------------------------------------

void print_declare_relation( std::ostream & msg ,
                             const char * method ,
                             Entity & e_from ,
                             Entity & e_to ,
                             const unsigned identifier ,
                             const unsigned kind )
{
  msg << method ;
  msg << "( " ;
  print_entity_key( msg , e_from.key() );
  msg << " , " ;
  print_entity_key( msg , e_to.key() );
  msg << " , " ;
  msg << identifier ;
  msg << " , " ;
  msg << kind ;
  msg << " )" ;
}

}

void BulkData::declare_relation( Entity & e_from ,
                                     Entity & e_to ,
                                     const unsigned identifier ,
                                     const unsigned kind )
{
  static const char method[] = "phdmesh::BulkData::declare_relation" ;

  if ( in_closure( e_to , e_from ) ) {
    std::ostringstream msg ;
    print_declare_relation( msg , method , e_from , e_to , identifier , kind );
    msg << " FAILED DUE TO CIRCULAR CLOSURE." ;
    throw std::runtime_error( msg.str() );
  }

  {
    const Relation forward( e_to , identifier , kind , false );
    const std::vector<Relation>::iterator e = e_from.m_relation.end();
          std::vector<Relation>::iterator i = e_from.m_relation.begin();

    i = std::lower_bound( i , e , forward , LessRelation() );

    if ( e == i || forward != *i ) { // Not a duplicate

      if ( e != i && forward.attribute() == i->attribute() ) {
        std::ostringstream msg ;
        print_declare_relation( msg, method, e_from, e_to, identifier, kind );
        msg << " FAILED, ALREADY HAS THIS RELATION TO " ;
        print_entity_key( msg , i->entity()->key() );
        throw std::runtime_error(msg.str());
      }

      e_from.m_relation.insert( i , forward );
    }
  }

  {
    const Relation converse( e_from , identifier , kind , true );
    const std::vector<Relation>::iterator e = e_to.m_relation.end();
          std::vector<Relation>::iterator i = e_to.m_relation.begin();

    i = std::lower_bound( i , e , converse , LessRelation() );

    if ( e == i || converse != *i ) { // Not a duplicate
      e_to.m_relation.insert( i , converse );
    }
  }

  {
    PartSet add , del ;

    deduce_part_relations( e_from , e_to , identifier , kind , add );

    internal_change_entity_parts( e_to , add , del );
  }

  set_field_relations( e_from , e_to , identifier , kind );
}

void BulkData::declare_relation( Entity & entity ,
                                     const std::vector<Relation> & rel )
{
  std::vector<Relation>::const_iterator i ;
  for ( i = rel.begin() ; i != rel.end() ; ++i ) {
    Entity & e = * i->entity();
    unsigned k = i->kind();
    unsigned n = i->identifier();
    if ( i->forward() ) {
      declare_relation( entity , e , n , k );
    }
    else {
      declare_relation( e , entity , n , k );
    }
  }
}

//----------------------------------------------------------------------

namespace {

void clear_field_relations( Entity & e_from ,
                            const EntityType type ,
                            const unsigned ident ,
                            const unsigned kind )
                          
{
  const std::vector<FieldRelation> & field_rels =
    e_from.kernel().mesh().mesh_meta_data().get_field_relations();

  for ( std::vector<FieldRelation>::const_iterator 
        j = field_rels.begin() ; j != field_rels.end() ; ++j ) {

    const FieldRelation & fr = *j ;

    void ** const ptr = (void**) field_data( * fr.m_root , e_from );

    if ( ptr ) {

      const int number =
        field_data_size(*fr.m_root,e_from) / sizeof(void*);

      const int offset =
        (*fr.m_function)( e_from.entity_type() , type , ident , kind );

      if ( 0 <= offset && offset < number ) {
        ptr[ offset ] = NULL ;
      }
    }
  }
}

}

void BulkData::destroy_relation( Entity & e1 , Entity & e2 , unsigned kind )
{
  // When removing a relationship may need to
  // remove part membership and set field relation pointer to NULL

  PartSet del ;
  bool to_e1 = false ;

  std::vector<Relation>::iterator i ;

  for ( i = e1.m_relation.end() ; i != e1.m_relation.begin() ; ) {
    --i ;
    if ( i->entity() == & e2 && i->kind() == kind ) {
      if ( i->forward() ) {
        clear_field_relations( e1 , i->entity_type() ,
                                    i->identifier() ,
                                    i->kind() );
        deduce_part_relations( e1, e2, i->identifier(), i->kind(), del );
      }
      i = e1.m_relation.erase( i );
    }
  }

  for ( i = e2.m_relation.end() ; i != e2.m_relation.begin() ; ) {
    --i ;
    if ( i->entity() == & e1 && i->kind() == kind ) {
      if ( i->forward() ) {
        clear_field_relations( e2 , i->entity_type() ,
                                    i->identifier() ,
                                    i->kind() );
        deduce_part_relations( e2, e1, i->identifier(), i->kind(), del );
        to_e1 = true ;
      }
      i = e2.m_relation.erase( i );
    }
  }

  Entity & e_to = to_e1 ? e1 : e2 ;

  {
    PartSet keep ;

    deduce_part_relations( e_to , keep );

    if ( ! keep.empty() ) {
      // Eliminate the 'keep' from the accumulated 'del'

      for ( PartSet::iterator j = del.end() ; j != del.begin() ; ) {
        --j ;
        if ( contain( keep , **j ) ) { j = del.erase( j ); }
      }
    }
  }

  if ( ! del.empty() ) {
    PartSet add ;
    internal_change_entity_parts( e_to , add , del );
  }
}

//----------------------------------------------------------------------
// Deduce propagation of changes to a part to the related 'to' entities

void BulkData::internal_propagate_part_changes(
  Entity        & entity ,
  const PartSet & removed )
{
  const EntityType etype = entity.entity_type();
  Part * const owns_part = & m_mesh_meta_data.locally_owned_part();
  Part * const uses_part = & m_mesh_meta_data.locally_used_part();

  PairIterRelation rel = entity.relations();

  for ( ; rel ; ++rel ) {
    const unsigned rel_ident = rel->identifier();
    const unsigned rel_kind  = rel->kind();

    if ( rel->forward() ) {

      Entity & e_to = * rel->entity();

      PartSet to_del ;
      PartSet to_add ;

      if ( ! removed.empty() ) {

        const EntityType t_to = e_to.entity_type();

        // Deduce parts for 'e_to' from all upward relations.
        // Any non-parallel part that I removed that is not deduced for
        // 'e_to' must be removed from 'e_to'

        deduce_part_relations( e_to , to_add );

        to_del.reserve( removed.size() );

        for ( PartSet::const_iterator
              j = removed.begin() ; j != removed.end() ; ++j ) {
          Part * const p = *j ;

          if ( p != owns_part && p != uses_part && ! contain( to_add , *p ) ) {

            to_del.push_back( p );

            // What if removing a part with a part-relation ?

            const std::vector<PartRelation> & part_rel =
              m_mesh_meta_data.get_part_relations();

            for ( std::vector<PartRelation>::const_iterator
                  k = part_rel.begin() ; k != part_rel.end() ; ++k ) {

              const PartRelation & stencil = *k ;

              if ( p == stencil.m_root &&
                   0 <= (*stencil.m_function)(etype,t_to,rel_ident,rel_kind) &&
                   ! contain( to_add , * stencil.m_target ) ) {
              }
            }
          }
        }
      }
      else {
        deduce_part_relations( entity , e_to , rel_ident , rel_kind , to_add );
      }

      internal_change_entity_parts( e_to , to_add , to_del );

      set_field_relations( entity, e_to, rel_ident , rel_kind );
    }
    else {
      Entity & e_from = * rel->entity();

      set_field_relations( e_from, entity, rel_ident, rel_kind );
    }
  }
}

void BulkData::internal_propagate_relocation( Entity & entity )
{
  PairIterRelation rel = entity.relations();

  for ( ; rel ; ++rel ) {
    if ( rel->forward() ) {
      Entity & e_to = * rel->entity();

      set_field_relations( entity, e_to, rel->identifier(), rel->kind() );
    }
    else {
      Entity & e_from = * rel->entity();

      set_field_relations( e_from, entity, rel->identifier(), rel->kind() );
    }
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

