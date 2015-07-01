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

#include <ctype.h>
#include <strings.h>

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include <util/NamedValue.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

bool less_nocase( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) < 0 ; }

bool equal_nocase( const std::string & lhs , const std::string & rhs )
{ return strcasecmp( lhs.c_str() , rhs.c_str() ) == 0 ; }


bool value_ios_get_token( std::istream & s , int t , bool required )
{
  static const char method_name[] = "phdmesh::value_ios_get_token" ;

  // while ( s.good() && isspace( s.peek() ) ) { s.get(); }
   while ( isspace( s.peek() ) ) { s.get(); }

  const bool found = s.peek() == t ;

  if ( found ) {
    s.get();
  }
  else if ( required ) {
    const char c[2] = { (char) t , 0 };
    std::string msg ;
    msg.append( method_name )
       .append( " FAILED to find '" )
       .append( c )
       .append( "'" );
    throw std::runtime_error(msg);
  }
  return found ;
}

bool good_c_name( const char * name )
{
  bool result = isalpha( *name );

  while ( result && *++name ) {
    result = isalnum(*name) || *name == '_' ;
  }
  return result ;
}

int peek_non_space( std::istream & s )
{
  while ( s.good() && isspace( s.peek() ) ) { s.get(); }
  return s.peek();
}

std::string read_name( std::istream & s )
{
  std::string name ;
  int c = peek_non_space(s);

  if ( isalpha( c ) ) {
    while ( isalnum( c ) || c == '_' ) {
      const char tmp[2] = { (char) c , 0 };
      name.append( tmp );
      s.get();
      c = s.peek();
    }
  }

  return name ;
}

//----------------------------------------------------------------------

struct less_pset {
  bool operator()( const NamedValue<void> * lhs ,
                   const std::string & rhs ) const
    { return less_nocase( lhs->name , rhs ); }
};

std::vector< NamedValue<void>* >::const_iterator
lower_bound( const std::vector<NamedValue<void>* > & v , const std::string & n )
{
  return std::lower_bound( v.begin(), v.end(), n , less_pset() );
}

std::vector< NamedValue<void>* >::iterator
lower_bound( std::vector<NamedValue<void>*> & v , const std::string & n )
{
  return std::lower_bound( v.begin(), v.end(), n , less_pset() );
}

//----------------------------------------------------------------------

void verify_name_policy( const char * method_name , const std::string & name )
{
  if ( ! good_c_name( name.c_str() ) ) {
    std::string msg ;
    msg.append( method_name )
       .append( "( " )
       .append( name )
       .append( " ) FAILED: bad name" );
    throw std::runtime_error( msg );
  }
}

void verify_no_loop( const char * const method ,
                     const NamedValueSet * const vs ,
                     const NamedValue<void> & p ,
                     std::string path )
{
  if ( p.type == typeid(NamedValueSet) ) {
    const NamedValueSet * ps =
      & static_cast< const NamedValue<NamedValueSet> & >(p).value ;

    if ( vs == ps ) {
      std::string msg ;
      msg.append( method )
         .append(" FAILED, Loop attempted: " )
         .append( p.name )
         .append( "." )
         .append( path )
         .append( p.name );
      throw std::runtime_error( msg );
    }
    else {
      path.append( p.name ).append( "." );

      const std::vector<NamedValue<void>*> & vec = ps->get();

      for ( std::vector<NamedValue<void>*>::const_iterator
            i = vec.begin() ; i != vec.end() ; ++i ) {
        verify_no_loop( method , vs , **i , path );
      }
    }
  }
}

//----------------------------------------------------------------------

void remove_this( std::vector< NamedValueSet *> & v , NamedValueSet * const ps )
{
  if ( ps ) {
    std::vector<NamedValueSet*>::iterator i ;
    for ( i = v.begin() ; i != v.end() && ps != *i ; ++i ) {}
    if ( i != v.end() ) { v.erase( i ); }
  }
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

NamedValue<void>::~NamedValue()
{
  while ( ! m_holders.empty() ) {
    NamedValueSet & ps = *m_holders.back() ;
    ps.remove( this );
  }
}

void NamedValue<void>::get_throw( const std::type_info & t , unsigned i ) const 
{
  std::ostringstream msg ;
  msg << "phdmesh::NamedValue< " ;
  tell( msg );
  msg << " >::get<" << t.name() << ">(" << i << ") FAILED" ;
  throw std::runtime_error( msg.str() );
}

void NamedValue<void>::put_throw( const std::type_info & t , unsigned i ) const 
{
  std::ostringstream msg ;
  msg << "phdmesh::NamedValue< " ;
  tell( msg );
  msg << " >::put<" << t.name()<< ">(" << i << ") FAILED" ;
  throw std::runtime_error( msg.str() );
}

//----------------------------------------------------------------------

NamedValueSet::~NamedValueSet()
{ clear(); }

NamedValueSet::NamedValueSet() : m_members()
{}

void NamedValueSet::clear()
{
  NamedValueSet * const myself = this ;
  while ( ! m_members.empty() ) {
    NamedValue<void> & v = * m_members.back();
    remove_this( v.m_holders , myself );
    m_members.pop_back();
  }
}

//----------------------------------------------------------------------

NamedValue<void> *
NamedValueSet::find( const std::string & n ,
                     const char sep ) const
{
  NamedValue<void> * v = NULL ;

  const std::string::size_type len = n.size();
  const std::string::size_type p   = sep ? n.find( sep ) : len ;

  if ( len == p || std::string::npos == p ) {
    // Local:
    const std::vector<NamedValue<void>*>::const_iterator
      i = lower_bound( m_members , n );

    if ( m_members.end() != i && equal_nocase( n , (*i)->name ) ) { v = *i ; }
  }
  else {
    // Nested:
    std::string local_name  = n.substr( 0 , p ); // before sep
    std::string nested_name = n.substr( p + 1 ); // after  sep

    const std::vector<NamedValue<void>*>::const_iterator
      i = lower_bound( m_members , local_name );

    NamedValueSet & nested =
      static_cast< NamedValue<NamedValueSet> *>( *i )->value ;

    v = nested.find( nested_name , sep );
  }

  return v ;
}

//----------------------------------------------------------------------

NamedValue<void> *
NamedValueSet::insert( NamedValue<void> * v )
{
  static const char method_name[] = "phdmesh::NamedValueSet::insert" ;

  NamedValue<void> * m = NULL ;

  verify_name_policy( method_name , v->name );

  verify_no_loop( method_name , this , *v , std::string() );

  std::vector<NamedValue<void>*>::iterator i =
    lower_bound( m_members , v->name );

  if ( m_members.end() != i && equal_nocase( v->name , (*i)->name ) ) {
    m = *i ;
  }
  else {
    NamedValueSet * const myself = this ;
    v->m_holders.push_back( myself );
    m = v ;
    m_members.insert( i , m );
  }

  return m ;
}

void NamedValueSet::remove( NamedValue<void> * p )
{
  std::vector<NamedValue<void>*>::iterator
    i = lower_bound( m_members , p->name );

  if ( m_members.end() != i && *i == p ) {
    m_members.erase( i );
    remove_this( p->m_holders , this );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

std::istream & operator >> ( std::istream & s , NamedValueSet & v )
{
  static const char method[] =
    "phdmesh::operator >> ( std::istream & , NamedValueSet & )" ;

  if ( value_ios_get_token( s , '{' , false ) ) {

    while ( ! value_ios_get_token( s , '}' , false ) ) {
      std::string ex_msg ;
      std::string name = read_name( s );

      NamedValue<void> * m ;

      if ( name.empty() ) {
        ex_msg.append(": failed to read <name>");
      }
      else if ( ! value_ios_get_token( s , '=' , false ) ) {
        ex_msg.append(": '" );
        ex_msg.append(name);
        ex_msg.append("' failed to read initial '='");
      }
      else if ( NULL == ( m = v.find( name ) ) ) {
        ex_msg.append(": '" );
        ex_msg.append(name);
        ex_msg.append("' is not a member");
      }
      else {
        m->read(s);
        if ( ! value_ios_get_token( s , ';' , false ) ) {
          ex_msg.append(": '" );
          ex_msg.append(name);
          ex_msg.append("' failed to read terminating ';'");
        }
      }

      if ( ex_msg.size() ) {
        std::string msg( method );
        msg.append( ex_msg );
        throw std::runtime_error( msg );
      }
    }
  }

  return s ;
}

std::ostream & operator << ( std::ostream & s , const NamedValueSet & v )
{
  const std::vector< NamedValue<void> *> & members = v.get();

  s << "{" ;

  if ( ! members.empty() ) {
    s << std::endl ;

    for ( std::vector< NamedValue<void> * >::const_iterator
          i = members.begin() ; i != members.end() ; ++i ) {

      s << (*i)->name << " = " ;
      (*i)->write( s );
      s << " ;" << std::endl ;
    }
  }

  s << "}" ;

  return s ;
}

//----------------------------------------------------------------------

void NamedValue< NamedValueSet >::tell( std::ostream & s ) const
{
  const std::vector< NamedValue<void> *> & members = value.get();

  s << name << " = NamedValueSet {" ;

  if ( ! members.empty() ) {

    s << std::endl ;

    for ( std::vector<NamedValue<void>*>::const_iterator
         i = members.begin() ; i != members.end() ; ++i ) {

      (*i)->tell( s );
      s << " ;" << std::endl ;
    }

    s << "}" ;
  }
}


//----------------------------------------------------------------------

NamedValue< NamedValueSet >::~NamedValue() {}

unsigned NamedValue< NamedValueSet >::pack( void * ) const { return 0 ; }
unsigned NamedValue< NamedValueSet >::unpack( void * ) { return 0 ; }


} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------


