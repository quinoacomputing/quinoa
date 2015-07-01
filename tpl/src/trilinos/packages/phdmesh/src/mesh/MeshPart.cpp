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

#include <mesh/Part.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/Entity.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

namespace {

inline
bool equal_no_case( const std::string & lhs ,
                    const std::string & rhs )
{
  const char * const lhs_c_str = lhs.c_str();
  const char * const rhs_c_str = rhs.c_str();
  return 0 == strcasecmp( lhs_c_str , rhs_c_str );
}

struct PartLess {

  inline
  bool operator()( const Part * lhs , const Part & rhs ) const
  {
    const unsigned l = lhs->mesh_meta_data_ordinal();
    const unsigned r = rhs.mesh_meta_data_ordinal();
    return l < r ;
  }

  inline
  bool operator()( const Part * lhs , const Part * rhs ) const
  {
    const unsigned l = lhs->mesh_meta_data_ordinal();
    const unsigned r = rhs->mesh_meta_data_ordinal();
    return l < r ;
  }
};

}

//----------------------------------------------------------------------

Part * find( const PartSet & parts , const std::string & name )
{
  PartSet::const_iterator i = parts.begin();

  while ( i != parts.end() && ! equal_no_case((*i)->name(),name) ) { ++i ; }

  return i != parts.end() ? *i : (Part*) NULL ;
}

//----------------------------------------------------------------------

bool verify( const Part & p , std::string & msg )
{
  bool ok = true ;

  // Superset/subset consistency

  const Part    & universal = p.mesh_meta_data().universal_part();
  const PartSet & supersets = p.supersets();
  const PartSet & subsets   = p.subsets();
  const PartSet & intersection = p.intersection_of();

  std::vector<Part*>::const_iterator i , j ;

  if ( & p == & universal ) {
    if ( ! supersets.empty() ) {
      ok = false ;
      msg.append(" ");
      msg.append( p.name() );
      msg.append( " cannot have supersets ;" );
    }
  }
  else {

    // Unversial superset with symmetry

    if ( ! contain( supersets , universal ) ) {
      ok = false ;
      msg.append(" ");
      msg.append( p.name() );
      msg.append( " not-in " );
      msg.append( universal.name() );
      msg.append( " ;" );
    }

    if ( ! contain( universal.subsets() , p ) ) {
      ok = false ;
      msg.append(" ");
      msg.append( universal.name() );
      msg.append( " not-contain " );
      msg.append( p.name() );
      msg.append( " ;" );
    }

    // Superset symmetry and noncircular

    for ( i = supersets.begin() ; i != supersets.end() ; ++i ) {
      Part & s = **i ;
      if ( ! contain( s.subsets() , p ) ) {
        ok = false ;
        msg.append( " Asymmetry " );
        msg.append( p.name() );
        msg.append( " in " );
        msg.append( s.name() );
        msg.append( " ;" );
      }

      if ( contain( subsets , s ) ) {
        ok = false ;
        msg.append( " Circular " );
        msg.append( p.name() );
        msg.append( " in " );
        msg.append( s.name() );
        msg.append( " ;" );
      }
    }

    // Subset symmetry, noncircular, and transitive

    for ( i = subsets.begin() ; i != subsets.end() ; ++i ) {
      Part & sub = **i ;
      if ( ! contain( sub.supersets() , p ) ) {
        ok = false ;
        msg.append( " Asymmetry " );
        msg.append( p.name() );
        msg.append( " contain " );
        msg.append( sub.name() );
        msg.append( " ;" );
      }
      if ( contain( supersets , sub ) ) {
        ok = false ;
        msg.append( " Circular " );
        msg.append( p.name() );
        msg.append( " contain " );
        msg.append( sub.name() );
        msg.append( " ;" );
      }
      for ( j = supersets.begin() ; j != supersets.end() ; ++j ) {
        Part & s = **j ;
        if ( ! contain( s.subsets() , sub ) ) {
          ok = false ;
          msg.append( " Not-transitive " );
          msg.append( s.name() );
          msg.append( " contain " );
          msg.append( p.name() );
          msg.append( " contain " );
          msg.append( sub.name() );
          msg.append( " ;" );
        }
      }
    }
  }

  // Intersection is subclass of supersets.
  // Intersection members cannot be subset/superset of one another.

  for ( i = intersection.begin() ; i != intersection.end() ; ) {
    Part & a = **i ; ++i ;
    if ( ! contain( supersets , a ) ) {
      ok = false ;
      msg.append( " Intersection-superset " );
      msg.append( p.name() );
      msg.append( " not-in " );
      msg.append( a.name() );
      msg.append( " ;" );
    }

    for ( j = i ; j != intersection.end() ; ) {
      Part & b = **j ; ++j ;
      if ( contain( a.supersets() , b ) ||
           contain( b.supersets() , a ) ) {
        ok = false ;
        msg.append( " Intersection " );
        msg.append( p.name() );
        msg.append( " ( " );
        msg.append( a.name() );
        msg.append( " , " );
        msg.append( b.name() );
        msg.append( " );" );
      }
    }
  }

  return ok ;
}

//----------------------------------------------------------------------

std::ostream &
print( std::ostream & os , const char * const lead , const Part & p )
{
  const PartSet & supersets = p.supersets();
  const PartSet & subsets   = p.subsets();
  const PartSet & intersection = p.intersection_of();

  std::vector<Part*>::const_iterator i ;

  if ( lead != NULL ) { os << lead ; }
  os << "Part[ " ;
  os << p.name() ;
  os << " , " ;
  os << p.mesh_meta_data_ordinal() ;
  os << " ] {" ;
  os << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Supersets {" ;
  for ( i = supersets.begin() ; i != supersets.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " }" << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Intersection_Of {" ;
  for ( i = intersection.begin() ; i != intersection.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " } }" << std::endl ;

  if ( lead != NULL ) { os << lead ; }
  os << "  Subsets {" ;
  for ( i = subsets.begin() ; i != subsets.end() ; ++i ) {
    const std::string & n = (*i)->name() ; os << " " << n ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

void order( PartSet & v )
{
  PartSet::iterator ev = v.end();
  PartSet::iterator iv = v.begin();
  std::sort( iv , ev , PartLess() );
  iv = std::unique( iv , ev );
  v.erase( iv , ev );
}

bool insert( PartSet & v , Part & part )
{
  const PartSet::iterator e = v.end();
        PartSet::iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  const bool new_member = i == e || *i != & part ;

  if ( new_member ) { Part * const tmp = & part ; v.insert( i , tmp ); }
  return new_member ;
}

bool contain( const PartSet & v , const Part & part )
{
  const PartSet::const_iterator e = v.end();
        PartSet::const_iterator i = v.begin();

  i = std::lower_bound( i , e , part , PartLess() );

  return i != e && *i == & part ;
}

bool contain( const PartSet & super , const PartSet & sub )
{
  bool result = ( ! sub.empty() ) && ( sub.size() <= super.size() );

  if ( result ) {
    PartLess comp ;

    const PartSet::const_iterator ev = super.end();
          PartSet::const_iterator iv = super.begin();

    const PartSet::const_iterator ep = sub.end();
          PartSet::const_iterator ip = sub.begin();

    while ( result && ip != ep ) {
      Part * const q = *ip ; ++ip ;
      iv = std::lower_bound( iv , ev , q , comp );
      result = iv != ev && *iv == q ; 
    }
  }

  return result ;
}

unsigned intersect( const PartSet & v , const PartSet & p )
{
  // Both lists must be sorted, assume v.size() > p.size()

  const PartSet::const_iterator ev = v.end();
        PartSet::const_iterator iv = v.begin();

  const PartSet::const_iterator ep = p.end();
        PartSet::const_iterator ip = p.begin();

  unsigned count = 0 ;

  for ( ; ip != ep && iv != ev ; ++ip ) {
    Part * const q = *ip ;
    iv = std::lower_bound( iv , ev , q , PartLess() );
    if ( iv != ev && *iv == q ) { ++count ; }
  }

  return count ;
}

unsigned intersect( const PartSet & v , const PartSet & p , PartSet & r )
{
  // Both lists must be sorted, assume v.size() > p.size()

  const PartSet::const_iterator ev = v.end();
        PartSet::const_iterator iv = v.begin();

  const PartSet::const_iterator ep = p.end();
        PartSet::const_iterator ip = p.begin();

  for ( ; ip != ep && iv != ev ; ++ip ) {
    Part * const q = *ip ;
    iv = std::lower_bound( iv , ev , q , PartLess() );
    if ( iv != ev && *iv == q ) { r.push_back( q ); }
  }

  return r.size() ;
}

bool intersect( const Part & a , const Part & b )
{
  const PartSet & a_sub = a.subsets();
  const PartSet & b_sub = b.subsets();
  return contain( a_sub , b ) ||
         contain( b_sub , a ) ||
         intersect( b_sub , a_sub );
}

//----------------------------------------------------------------------

Part::~Part()
{}

// The constructor must only be called by
// the 'MetaData::declare_part( const std::string & )' method.

Part::Part( MetaData & m , const std::string & n , unsigned ordinal )
  : m_name( n ),
    m_cset(),
    m_subsets() , m_supersets() , m_intersect() ,
    m_mesh_meta_data( m ) ,
    m_mesh_meta_data_ordinal( ordinal ),
    m_relation_target( false )
{}

Part & MetaData::declare_part( const std::string & p_name )
{
  static const char method[] = "phdmesh::BulkData::declare_part" ;

  assert_not_committed( method );

  Part * const u = & m_universal_part ;

  Part * p = find( u->m_subsets , p_name );

  if ( p == NULL ) {
    const unsigned ord = m_universal_part.m_subsets.size();

    p = new Part( *this , p_name , ord );

    u->m_subsets  .push_back( p );
    p->m_supersets.push_back( u );
  }

  return *p ;
}

//----------------------------------------------------------------------

namespace {

void assert_not_relation_target( const char * const method ,
                                 const char * const what ,
                                 const Part & part )
{
  if ( part.is_relation_target() ) {
    std::string msg ;
    msg.append( method );
    msg.append(" : FAILED, " );
    msg.append( part.name() );
    msg.append(" is a part relation target and cannot be a " );
    msg.append( what );
    throw std::invalid_argument(msg);
  }
}

}

//----------------------------------------------------------------------

namespace {

void clean_intersection( const char * const method ,
                         PartSet     & pset ,
                         std::string & name )
{
  static const char separator[] = "^" ;

  order( pset );

  PartSet::iterator i ;

  for ( i = pset.begin() ; i != pset.end() ; ) {
    // If a subset of 'i' is contained then 'i' is redundant
    if ( intersect( (*i)->subsets() , pset ) ) {
      i = pset.erase( i );
    }
    else {
      ++i ;
    }
  }

  if ( pset.size() < 2 ) {
    std::string msg ;
    msg.append(method);
    msg.append(" : FAILED, Cannot intersect fewer than two unique parts." );
    msg.append(" Input {" );
    PartSet::iterator j ;
    for ( j = pset.begin() ; j != pset.end() ; ++j ) {
      msg.append(" ");
      msg.append( (*j)->name() );
    }
    msg.append(" } Clean {" );
    for ( i = pset.begin() ; i != pset.end() ; ) {
      msg.append(" ");
      msg.append( (*i)->name() );
    }
    throw std::invalid_argument(msg);
  }

  name.assign("{");
  for ( i = pset.begin() ; i != pset.end() ; ++i ) {
    if ( i != pset.begin() ) { name.append( separator ); }
    name.append( (*i)->name() );
  }
  name.append("}");
}

}

Part & MetaData::declare_part( const PartSet & pset )
{
  static const char method[] = "phdmesh::BulkData::declare_part" ;

  assert_not_committed( method );

  for ( PartSet::const_iterator i = pset.begin() ; i != pset.end() ; ++i ) {
    assert_not_relation_target( method , "superset" , **i );
  }

  PartSet pset_clean( pset );

  std::string p_name ;

  clean_intersection( method , pset_clean , p_name );

  Part * p = find( m_universal_part.m_subsets , p_name );

  if ( p == NULL ) {
    Part & p_new = declare_part( p_name );
    p_new.m_intersect = pset_clean ; // Copy
    PartSet::iterator i ;
    for ( i = pset_clean.begin() ; i != pset_clean.end() ; ++i ) {
      declare_part_subset( **i , p_new );
    }
    p = & p_new ;
  }
  else if ( pset_clean != p->m_intersect ) {
    std::string msg ;
    msg.append(method);
    msg.append(" : FAILED, Redundant incompatible intersection declaration.");
    throw std::invalid_argument(msg);
  }

  return *p ;
}

//----------------------------------------------------------------------

void MetaData::declare_part_subset( Part & superset , Part & subset )
{
  static const char method[] = "phdmesh::MetaData::declare_part_subset" ;

  assert_not_relation_target( method , "superset" , superset );
  assert_not_relation_target( method , "subset" , subset );

  if ( ! contain( superset.m_subsets , subset ) ) {

    assert_not_committed( method );
    assert_same_mesh_meta_data(   method , superset.mesh_meta_data() );
    assert_same_mesh_meta_data(   method , subset.mesh_meta_data() );

    if ( & m_universal_part == & subset ||
         & superset         == & subset ||
         contain( superset.m_supersets , subset ) ) {
      std::string msg ;
      msg.append( method )
         .append( "[ " )
         .append( superset.name() )
         .append( " ] ( " )
         .append( subset.name() )
         .append( " ) FAILED, IS CIRCULAR" );
      throw std::invalid_argument( msg );
    }

    // Symmetry:

    insert( subset.m_supersets , superset );
    insert( superset.m_subsets , subset );

    PartSet::iterator i ;

    // Transitive, is also a subset of superset's supersets:

    for ( i =  superset.m_supersets.begin() ;
          i != superset.m_supersets.end() ; ++i ) {
      declare_part_subset( **i , subset );
    }

    // Intersection-part membership check.

    for ( i =  superset.m_subsets.begin() ;
          i != superset.m_subsets.end() ; ++i ) {
      Part & p_sub = **i ;
      if ( & p_sub != & subset ) {

        // If subset is fully contained in a superset->subset
        // intersection-part then is also a subset of that
        // intersection-part.

        // If subset is an intersection-part and one of the superset's
        // subsets is fully contained in the subset's intersection then
        // that superset-subset is a subset of this subset.

        if ( contain( subset.m_supersets , p_sub.m_intersect ) ) {
          declare_part_subset( p_sub , subset );
        }
        else if ( contain( p_sub.m_supersets , subset.m_intersect ) ) {
          declare_part_subset( subset , p_sub );
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void MetaData::declare_part_relation(
  Part & root_part ,
  relation_stencil_ptr stencil ,
  Part & target_part )
{
  static const char method[] = "phdmesh::MetaData::declare_part_relation" ;

  assert_not_relation_target( method , "root of part relation" , root_part );

  if ( 0 != target_part.subsets().size() ||
       0 != target_part.intersection_of().size() ||
       1 != target_part.supersets().size() ) {
    std::string msg ;
    msg.append( method );
    msg.append( ": FAILED, target " );
    msg.append( target_part.name() );
    msg.append( " cannot be a superset or subset" );
    throw std::runtime_error( msg );
  }

  target_part.m_relation_target = true ;

  PartRelation tmp ;
  tmp.m_root = & root_part ;
  tmp.m_target = & target_part ;
  tmp.m_function = stencil ;

  m_part_relations.push_back( tmp );
}

//----------------------------------------------------------------------

Part * MetaData::get_part( const std::string & p_name ,
                         const char * required_by ) const
{
  const PartSet & all_parts = m_universal_part.m_subsets ;

  Part * const p = find( all_parts , p_name );

  if ( required_by && NULL == p ) { // ERROR
    static const char method[] = "phdmesh::BulkData::get_part" ;
    std::string msg ;
    msg.append( method )
       .append( "( " )
       .append( p_name )
       .append( " , " )
       .append( required_by )
       .append( " ) FAILED to find part" );
    throw std::runtime_error( msg );
  }

  return p ;
}

//----------------------------------------------------------------------

} // namespace phdmesh

