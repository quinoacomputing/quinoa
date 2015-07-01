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

#include <string.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/Comm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

void MetaData::assert_not_committed( const char * method ) const
{
  if ( m_commit ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED: mesh MetaData has been committed." );
    throw std::logic_error( msg );
  }
}

void MetaData::assert_committed( const char * method ) const
{
  if ( ! m_commit ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED: mesh MetaData has not been committed." );
    throw std::logic_error( msg );
  }
}

void MetaData::assert_same_mesh_meta_data( const char * method ,
                                 const MetaData & rhs ) const
{
  if ( this != & rhs ) {
    std::string msg ;
    msg.append( method )
       .append( " FAILED Different mesh_meta_data." );
    throw std::logic_error( msg );
  }
}

//----------------------------------------------------------------------

MetaData::MetaData()
  : m_commit( false ),
    m_universal_part( *this , std::string( "{UNIVERSAL}" ) , 0 ),
    m_uses_part( NULL ),
    m_owns_part( NULL )
{
  // Declare remaining predefined parts
  const std::string uses_part_name( "{USES}" );
  const std::string owns_part_name( "{OWNS}" );

  {
    Part * u = & m_universal_part ;
    m_universal_part.m_subsets.push_back( u );
  }

  m_uses_part = & declare_part( uses_part_name );
  m_owns_part = & declare_part( owns_part_name );

  declare_part_subset( * m_uses_part , * m_owns_part );
}

void MetaData::commit()
{
  static const char method[] = "phdmesh::MetaData::commit" ;

  assert_not_committed( method );

  { // Verify parts
    std::string msg ;
    msg.append( method );
    msg.append( " FAILED : " );
    const std::vector<Part*> & parts = m_universal_part.subsets();
    for ( unsigned i = 0 ; i < parts.size() ; ++i ) {
      Part & p = * parts[i] ;
      if ( ! verify( p , msg ) ) {
        throw std::logic_error( msg );
      }
    }
  }

  clean_field_restrictions();

  m_commit = true ; // Cannot add or change parts or fields now
}

MetaData::~MetaData()
{
  // Destroy the fields, used 'new' to allocate so now use 'delete'

  {
    std::vector<FieldBase * >::iterator j = m_fields.begin();

    for ( ; j != m_fields.end() ; ++j ) { delete *j ; }

    m_fields.clear();
  }

  // Destroy the parts, used 'new' to allocate so now use 'delete'
  {
    std::vector<Part*> & parts = m_universal_part.m_subsets ;

    std::vector< Part * >::iterator j = parts.begin();

    for ( ; j != parts.end() ; ++j ) {
      if ( *j != & m_universal_part ) { delete *j ; }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Verify parallel consistency of fields and parts

namespace {

void pack( CommBuffer & b , const PartSet & pset )
{
  PartSet::const_iterator i , j ;
  for ( i = pset.begin() ; i != pset.end() ; ++i ) {
    const Part & p = **i ;
    const PartSet & subsets   = p.subsets();
    const PartSet & intersect = p.intersection_of();

    const unsigned     name_len = p.name().size() + 1 ;
    const char * const name_ptr = p.name().c_str();

    {
      const unsigned ord = p.mesh_meta_data_ordinal();
      b.pack<unsigned>( ord );
    }

    b.pack<unsigned>( name_len );
    b.pack<char>( name_ptr , name_len );

    const unsigned subset_size = subsets.size();
    b.pack<unsigned>( subset_size );
    for ( j = subsets.begin() ; j != subsets.end() ; ++j ) {
      const Part & s = **j ;
      const unsigned ord = s.mesh_meta_data_ordinal();
      b.pack<unsigned>( ord );
    }
    const unsigned intersect_size = intersect.size();
    b.pack<unsigned>( intersect_size );
    for ( j = intersect.begin() ; j != intersect.end() ; ++j ) {
      const Part & s = **j ;
      const unsigned ord = s.mesh_meta_data_ordinal();
      b.pack<unsigned>( ord );
    }
  }
}

bool unpack_verify( CommBuffer & b , const PartSet & pset )
{
  enum { MAX_TEXT_LEN = 4096 };
  char b_text[ MAX_TEXT_LEN ];
  unsigned b_tmp ;

  bool ok = true ;
  PartSet::const_iterator i , j ;
  for ( i = pset.begin() ; ok && i != pset.end() ; ++i ) {
    const Part & p = **i ;
    const PartSet & subsets   = p.subsets();
    const PartSet & intersect = p.intersection_of();
    const unsigned     name_len = p.name().size() + 1 ;
    const char * const name_ptr = p.name().c_str();

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == p.mesh_meta_data_ordinal();
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == name_len ;
    }
    if ( ok ) {
      b.unpack<char>( b_text , name_len );
      ok = 0 == strcmp( name_ptr , b_text );
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == subsets.size() ;
    }
    for ( j = subsets.begin() ; ok && j != subsets.end() ; ++j ) {
      const Part & s = **j ;
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == s.mesh_meta_data_ordinal();
    }

    if ( ok ) {
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == intersect.size();
    }
    for ( j = intersect.begin() ; ok && j != intersect.end() ; ++j ) {
      const Part & s = **j ;
      b.unpack<unsigned>( b_tmp );
      ok = b_tmp == s.mesh_meta_data_ordinal();
    }
  }
  return ok ;
}

void pack( CommBuffer & ,
           const std::vector< FieldBase * > & )
{
}

bool unpack_verify( CommBuffer & ,
                    const std::vector< FieldBase * > & )
{
  bool ok = true ;
  return ok ;
}

}

//----------------------------------------------------------------------

void verify_parallel_consistency( const MetaData & s , ParallelMachine pm )
{
  static const char method[] = "phdmesh::verify_parallel_consistency(MetaData)" ;

  const unsigned p_rank = parallel_machine_rank( pm );

  const bool is_root = 0 == p_rank ;

  CommBroadcast comm( pm , 0 );

  if ( is_root ) {
    pack( comm.send_buffer() , s.get_parts() );
    pack( comm.send_buffer() , s.get_fields() );
  }

  comm.allocate_buffer();

  if ( is_root ) {
    pack( comm.send_buffer() , s.get_parts() );
    pack( comm.send_buffer() , s.get_fields() );
  }

  comm.communicate();

  int ok[ 2 ];

  ok[0] = unpack_verify( comm.recv_buffer() , s.get_parts() );
  ok[1] = unpack_verify( comm.recv_buffer() , s.get_fields() );

  all_reduce( pm , Min<2>( ok ) );

  if ( ! ok[0] || ! ok[1] ) {
    std::ostringstream msg ;
    msg << "P" << p_rank ;
    msg << ": " << method ;
    msg << " : FAILED for:" ;
    if ( ! ok[0] ) { msg << " Parts" ; }
    if ( ! ok[1] ) { msg << " Fields" ; }
    throw std::logic_error( msg.str() );
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

