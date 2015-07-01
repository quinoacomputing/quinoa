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

#include <strings.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <util/SimpleArrayOps.hpp>
#include <mesh/Types.hpp>
#include <mesh/Field.hpp>
#include <mesh/Part.hpp>
#include <mesh/MetaData.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

const EntityDimension & EntityDimension::tag()
{ static const EntityDimension self ; return self ; }

const char * EntityDimension::name() const
{ static const char n[] = "EntityDimension" ; return n ; }

//----------------------------------------------------------------------

const char * field_state_name( FieldState s )
{
  static const char * name_list[] = {
    "StateNew" ,
    "StateOld" ,
    "StateNM1" ,
    "StateNM2" ,
    "StateNM3" ,
    "StateNM4" ,
    "ERROR" };

  unsigned i = s ;
  if ( StateNM4 < i ) { i = MaximumFieldStates ; }
  return name_list[i] ;
}

//----------------------------------------------------------------------

namespace {

struct DimLess {
  bool operator()( const FieldBase::Restriction & lhs ,
                   const unsigned rhs ) const
    { return lhs < rhs ; }
};

std::vector<FieldBase::Restriction>::const_iterator
find( const std::vector<FieldBase::Restriction> & v , unsigned key )
{
  std::vector<FieldBase::Restriction>::const_iterator
    i = std::lower_bound( v.begin() , v.end() , key , DimLess() );

  if ( i != v.end() && i->key != key ) { i = v.end(); }

  return i ;
}

void insert( std::vector<FieldBase::Restriction> & v ,
             const FieldBase::Restriction & d )
{
  std::vector<FieldBase::Restriction>::iterator
    i = std::lower_bound( v.begin() , v.end() , d );

  if ( i == v.end() || i->key != d.key ) { v.insert( i , d ); }
}

}

//----------------------------------------------------------------------

FieldBase::~Field()
{ }

FieldBase::Field(
  MetaData &            arg_mesh_meta_data ,
  const std::string & arg_name ,
  unsigned scalar_type ,
  unsigned rank ,
  const ArrayDimTag * const * tags ,
  unsigned number_of_states ,
  FieldState state )
: m_cset() ,
  m_name( arg_name ),
  m_mesh_meta_data( arg_mesh_meta_data ),
  m_mesh_meta_data_ordinal(0),
  m_scalar_type( scalar_type ),
  m_rank( rank ),
  m_num_states( number_of_states ),
  m_this_state( state ),
  m_dim_map()
{
  FieldBase * const pzero = NULL ;
  Copy<MaximumFieldStates>( m_field_states , pzero );
  unsigned i = 0 ;
  for ( ; i < rank ; ++i ) { m_dim_tags[i] = tags[i] ; }
  for ( ; i < MaximumFieldDimension ; ++i ) { m_dim_tags[i] = NULL ; }
}

const std::vector<FieldBase::Restriction> & FieldBase::restrictions() const
{ return m_field_states[0]->m_dim_map ; }

std::vector<FieldBase::Restriction> & FieldBase::restrictions()
{ return m_field_states[0]->m_dim_map ; }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

bool equal_field_name( const FieldBase * lhs , const std::string & rhs )
{
  const char * const l_name = lhs->name().c_str();
  const char * const r_name = rhs.c_str();
  return 0 == strcasecmp( l_name , r_name );
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void
print_field_type( std::ostream          & arg_msg ,
                  unsigned                arg_scalar_type ,
                  unsigned                arg_rank ,
                  const ArrayDimTag * const * arg_tags )
{
  arg_msg << "Field< " ;
  arg_msg << NumericEnum<>::name( arg_scalar_type );
  for ( unsigned i = 0 ; i < arg_rank ; ++i ) {
    arg_msg << " , " << arg_tags[i]->name();
  }
  arg_msg << " >" ;
}

}

FieldBase *
MetaData::get_field_base(
  const std::string & arg_name ,
  unsigned          arg_scalar_type ,
  unsigned          arg_rank ,
  const ArrayDimTag * const * arg_tags ,
  int arg_number_states ,
  const char * arg_required_by ) const
{
  static const char declare_method[] = "phdmesh::MetaData::declare_field" ;
  static const char get_method[]     = "phdmesh::MetaData::get_field" ;

  // Potential error conditions:

  bool not_found         = false ;
  bool bad_scalar_type   = false ;
  bool bad_dimension     = false ;
  bool bad_number_states = MaximumFieldStates < arg_number_states ;

  // Find the field by name:

  const std::vector< FieldBase * >::const_iterator e = m_fields.end();
        std::vector< FieldBase * >::const_iterator j = m_fields.begin();

  for ( ; j != e && ! equal_field_name( *j , arg_name ) ; ++j );

  FieldBase * const field = ( j != e ) ? *j : NULL ;

  // If found check for compatibility:

  if ( field != NULL ) {

    bad_scalar_type = arg_scalar_type != field->numeric_type_ordinal();

    bad_dimension = arg_rank != field->m_rank ;

    for ( unsigned i = 0 ; i < arg_rank && ! bad_dimension ; ++i ) {
      bad_dimension = arg_tags[i] != field->m_dim_tags[i] ;
    }

    bad_number_states =
      ( 0 <= arg_number_states ) &&
      ( arg_number_states != (int) field->m_num_states );
  }
  else if ( arg_required_by ) {
    not_found = true ;
  }

  if ( bad_scalar_type || bad_dimension || bad_number_states || not_found ) {

    const char * const method =
      arg_number_states < 0 ? get_method : declare_method ;

    std::ostringstream msg ;

    msg << method << "< " ;
    print_field_type( msg , arg_scalar_type , arg_rank , arg_tags );
    msg << " >( \"" << arg_name << "\"" ;

    if ( 0 <= arg_number_states ) {
      msg << " , " << arg_number_states ;
      if ( bad_number_states ) { msg << " <= NUMBER_OF_STATES IS BAD" ; }
    }

    if ( arg_required_by ) {
      msg << " , " << arg_required_by << " <= REQUIRED_BY" ;
      if ( not_found ) { msg << " BUT NOT FOUND" ; }
    }

    msg << " )" ;
 
    if ( field != NULL ) {
      msg << " FOUND INCOMPATIBLE " ;
      print_field_type( msg , field->numeric_type_ordinal() ,
                              field->m_rank ,
                              field->m_dim_tags );
      msg << "( " << field->m_num_states << " )" ;
    }
    throw std::runtime_error( msg.str() );
  }

  return field ;
}

//----------------------------------------------------------------------

FieldBase &
MetaData::declare_field_base(
  const std::string & arg_name ,
  unsigned            arg_scalar_type ,
  unsigned            arg_rank ,
  const ArrayDimTag * const * arg_tags ,
  unsigned            arg_num_states )
{
  static const char method[] = "phdmesh::MetaData::declare_field" ;

  static const char reserved_state_suffix[6][5] = {
    "_OLD" , "_N" , "_NM1" , "_NM2" , "_NM3" , "_NM4" };

  assert_not_committed( method );

  // Check that the name does not have a reserved suffix

  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    const int len = arg_name.size() - 4 ;
    if ( 0 <= len &&
         ! strcasecmp( arg_name.c_str() + len , reserved_state_suffix[i] ) ) {
      std::ostringstream msg ;
      msg << method << "< " ;
      print_field_type( msg , arg_scalar_type , arg_rank , arg_tags );
      msg << " >( \"" << arg_name ;
      msg << "\" <= HAS RESERVED STATE SUFFIX \"" ;
      msg << reserved_state_suffix[i] ;
      msg << " )" ;
      throw std::runtime_error( msg.str() );
    }
  }

  // Get field and if found verify compatibility

  FieldBase * field = get_field_base( arg_name ,
                                      arg_scalar_type ,
                                      arg_rank , arg_tags ,
                                      arg_num_states , NULL );

  if ( field == NULL ) {

    FieldBase * f[ MaximumFieldStates ] ;

    std::string field_names[ MaximumFieldStates ] ;

    field_names[0] = arg_name ;

    if ( 2 == arg_num_states ) {
      field_names[1] = arg_name ;
      field_names[1].append( reserved_state_suffix[0] );
    }
    else {
      for ( unsigned i = 1 ; i < arg_num_states ; ++i ) {
        field_names[i] = arg_name ;
        field_names[i].append( reserved_state_suffix[i] );
      }
    }

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {

      const std::vector< FieldBase * >::iterator e = m_fields.end();
            std::vector< FieldBase * >::iterator j = m_fields.begin();

      for ( ; j != e && ! equal_field_name( *j , field_names[i] ) ; ++j );

      if ( j != e ) {
        std::string msg( method );
        msg.append(" CATASTROPHIC INTERNAL LOGIC ERROR FOR FIELD NAME \"" );
        msg.append( field_names[i] );
        msg.append( "\"");
        throw std::logic_error( msg );
      }

      f[i] = new FieldBase( *this ,
                            field_names[i] ,
                            arg_scalar_type,
                            arg_rank , arg_tags ,
                            arg_num_states , (FieldState) i );

      f[i]->m_mesh_meta_data_ordinal = m_fields.size();

      m_fields.push_back( f[i] );
    }

    field = f[0] ;

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {
      FieldBase & tmp = * f[i] ;
      for ( unsigned k = 0 ; k < arg_num_states ; ++k ) {
        tmp.m_field_states[k] = f[k] ;
      }
    }
  }

  return *field ;
}

//----------------------------------------------------------------------

void MetaData::declare_field_relation(
  FieldBase & pointer_field ,
  relation_stencil_ptr stencil ,
  FieldBase & referenced_field )
{
  static const char method[] = "phdmesh::MetaData::declare_field_relation" ;

  static const int offset = NumericEnum<void*>::value -
                            NumericEnum<void>::value ;

  if ( referenced_field.numeric_type_ordinal() + offset !=
       pointer_field.numeric_type_ordinal() ) {

    std::string msg( method );
    msg.append( ": FAILED, " );
    msg.append( pointer_field.name() );
    msg.append( " and " );
    msg.append( referenced_field.name() );
    msg.append( " have incompatible numerical types." );
    throw std::invalid_argument( msg );
  }

  FieldRelation tmp ;
  tmp.m_root   = & pointer_field ;
  tmp.m_target = & referenced_field ;
  tmp.m_function = stencil ;

  m_field_relations.push_back( tmp );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void print_field_dim( std::ostream & msg ,
                      const FieldBase & f ,
                      const FieldBase::Restriction & d )
{
  const Part   & p = f.mesh_meta_data().get_part( d.ordinal() );
  const unsigned n = f.rank();

  msg << "( " << entity_type_name( d.type() );
  msg << " , " ;
  msg << p.name();
  msg << " , " << d.stride[0] ;
  for ( unsigned j = 1 ; j < n ; ++j ) {
    if ( d.stride[j-1] ) {
      msg << " , " << ( d.stride[j] / d.stride[j-1] );
    }
    else {
      msg << " , ERROR " ;
    }
  }
  msg << ")" ;
}

void assert_field_dimension_compatible(
  const char * const method ,
  const FieldBase & field ,
  const FieldBase::Restriction & a ,
  const FieldBase::Restriction & b )
{
  if ( NotEqual< MaximumFieldDimension >( a.stride , b.stride ) ) {
    std::ostringstream msg ;
    msg << method << " FOUND INCOMPATIBLE SIZES FOR " ;
    print_field_type( msg , field.numeric_type_ordinal() ,
                            field.rank(), field.dimension_tags() );
    msg << "[" << field.name() << "] " ;

    print_field_dim( msg , field , a );
    msg << " INCOMPATIBLE WITH " ;
    print_field_dim( msg , field , b );
    throw std::runtime_error( msg.str() );
  }
}

}

//----------------------------------------------------------------------
// Setting the dimension for one field sets the dimension
// for the corresponding fields of the FieldState array.
// If subset exists then replace it.
// If exists or superset exists then do nothing.

void MetaData::declare_field_restriction(
  FieldBase    & arg_field ,
  EntityType     arg_entity_type ,
  const Part   & arg_part ,
  const unsigned * arg_stride )
{
  static const char method[] =
    "phdmesh::MetaData::declare_field_restriction" ;

  assert_not_committed( method );
  assert_same_mesh_meta_data( method , arg_field.m_mesh_meta_data );
  assert_same_mesh_meta_data( method , arg_part.m_mesh_meta_data );

  FieldBase::Restriction tmp( arg_entity_type ,
                              arg_part.mesh_meta_data_ordinal() );

  if ( arg_field.m_rank ) {
    Copy<7>( tmp.stride , arg_stride );
  }
  else {
    tmp.stride[0] = 1 ; // For scalar fields
  }

  // If a dimension is defined for a superset of the part
  // then that dimension must be compatible and this
  // declaration is redundant.

  bool redundant  = false ;

  std::vector<FieldBase::Restriction> & dim_map =
    arg_field.m_field_states[0]->m_dim_map ;

  std::vector<FieldBase::Restriction>::iterator i ;

  for ( i = dim_map.begin() ; i != dim_map.end() ; ) {
    bool remove = false ;

    if ( arg_entity_type == i->type() ) {

      const Part & part = get_part( i->ordinal() );

      if ( ( arg_part == part ) || ( contain( arg_part.subsets() , part ) ) ) {
        redundant = true ;
      }
      else if ( contain( arg_part.supersets() , part ) ) {
        remove = true ;
      }

      if ( redundant || remove ) {
        assert_field_dimension_compatible( method , arg_field , tmp , *i );
      }
    }

    if ( remove ) {
      i = dim_map.erase( i );
    }
    else {
      ++i ;
    }
  }

  if ( ! redundant ) { insert( dim_map , tmp ); }
}

//----------------------------------------------------------------------
// If a part and one of its subset parts show up in the dimension map
// verify compatibility of dimensions and delete the subset part.

void MetaData::clean_field_restrictions()
{
  static const char method[] = "phdmesh::MetaData::clean_field_restrictions" ;
  const int zero = 0 ;

  for ( std::vector<FieldBase *>::iterator
        f = m_fields.begin() ; f != m_fields.end() ; ++f ) {
    FieldBase & field = **f ;

    std::vector<FieldBase::Restriction> & dim_map = field.m_dim_map ;
    std::vector<int> flag( dim_map.size() , zero );

    for ( unsigned i = 0 ; i < dim_map.size() ; ++i ) {
      const FieldBase::Restriction & dim = dim_map[i] ;
      const EntityType type = dim.type();
      const Part     & part = get_part( dim.ordinal() );
      const PartSet  & sub = part.subsets();

      for ( unsigned j = 0 ; j < dim_map.size() ; ++j ) {
        const Part     & p = get_part( dim_map[j].ordinal() );
        const EntityType t = dim_map[j].type();
        if ( i != j && type == t && contain( sub , p ) ) {
          assert_field_dimension_compatible( method, field, dim, dim_map[j] );
          flag[j] = 1 ;
        }
      }
    }

    for ( unsigned i = dim_map.size() ; i-- ; ) {
      if ( flag[i] ) {
        std::vector<FieldBase::Restriction>::iterator j = dim_map.begin();
        std::advance( j , i );
        dim_map.erase( j );
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// This part or any superset of this part

const FieldBase::Restriction &
FieldBase::restriction( EntityType etype , const Part & part ) const
{
  static const FieldBase::Restriction empty ;

  const std::vector<FieldBase::Restriction> & dim_map = restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator ie = dim_map.end() ;
        std::vector<FieldBase::Restriction>::const_iterator i ;

  const PartSet::const_iterator ipe = part.supersets().end();
        PartSet::const_iterator ip  = part.supersets().begin() ;

  unsigned key = Restriction::key_value(etype,part.mesh_meta_data_ordinal());

  while ( ie == ( i = find( dim_map , key ) ) && ipe != ip ) {
    key = Restriction::key_value( etype , (*ip)->mesh_meta_data_ordinal() );
    ++ip ;
  }

  return ie == i ? empty : *i ;
}

unsigned FieldBase::max_size( EntityType etype ) const
{
  unsigned max = 0 ;

  const std::vector<FieldBase::Restriction> & dim_map = restrictions();
  const std::vector<FieldBase::Restriction>::const_iterator ie= dim_map.end();
        std::vector<FieldBase::Restriction>::const_iterator i = dim_map.begin();

  for ( ; i != ie ; ++i ) {
    if ( i->type() == etype ) {
      const unsigned len = m_rank ? i->stride[ m_rank - 1 ] : 1 ;
      if ( max < len ) { max = len ; }
    }
  }

  return max ;
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const FieldBase & field )
{
  print_field_type( s , field.numeric_type_ordinal() ,
                        field.rank() , field.dimension_tags() );
  s << "[ \"" ;
  s << field.name() ;
  s << " \" , " ;
  s << field.number_of_states();
  s << " ]" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b ,
                      const FieldBase & field )
{
  const std::vector<FieldBase::Restriction> & dim_map = field.restrictions();
  s << field ;
  s << " {" ;
  for ( std::vector<FieldBase::Restriction>::const_iterator
        i = dim_map.begin() ; i != dim_map.end() ; ++i ) {
    s << std::endl << b << "  " ;
    print_field_dim( s , field , *i );
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------

} // namespace phdmesh

