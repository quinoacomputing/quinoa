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

#ifndef phdmesh_Part_hpp
#define phdmesh_Part_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>
#include <vector>

#include <util/CSet.hpp>
#include <mesh/Types.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

/** Supersets and subsets of parts are of the same mesh_meta_data and
 *  maintained in ordinal order.
 */
typedef std::vector<Part*> PartSet ;

/** A Part defines a subset of a to-be-discretized problem domain.
 *  A Part-defined subset can correspond to geometric subdomains,
 *  parallel distributed subdomains, types of discretization 
 *  entities (e.g., hexahedrons, tetrahedrons, ...), or any other
 *  subsetting need.
 */
class Part {
public:

  /** MetaData in which this part resides */
  MetaData & mesh_meta_data() const { return m_mesh_meta_data ; }

  /** Text name of this mesh part */
  const std::string & name() const { return m_name ; }

  /** Ordinal of this part within the mesh mesh_meta_data */
  unsigned mesh_meta_data_ordinal() const { return m_mesh_meta_data_ordinal ; }

  /** Parts that are supersets of this part.  */
  const PartSet & supersets() const { return m_supersets ; }

  /** Parts that are subsets of this part. */
  const PartSet & subsets() const { return m_subsets ; }

  /** Parts for which this part is defined to be the intersection */
  const PartSet & intersection_of() const { return m_intersect ; }

  /** A subset of entities may be deduced through relationships
   *  with other entities.
   */
  bool is_relation_target() const { return m_relation_target ; }

  bool operator == ( const Part & rhs ) const { return this == & rhs ; }
  bool operator != ( const Part & rhs ) const { return this != & rhs ; }

  template<class A>
  const A * attribute() const { return m_cset.template get<A>(); }

private:

  /* Only a MetaData can create and delete parts */
  friend class MetaData ;

  /** Construct a part within a given mesh */
  Part( MetaData & , const std::string & , unsigned );

  ~Part();
  Part();
  Part( const Part & );
  Part & operator = ( const Part & );

  const std::string m_name ;
  CSet              m_cset ;
  PartSet           m_subsets ;
  PartSet           m_supersets ;
  PartSet           m_intersect ;
  MetaData          & m_mesh_meta_data ;
  const unsigned    m_mesh_meta_data_ordinal ;
  bool              m_relation_target ;
};

//----------------------------------------------------------------------
/** Order a collection of parts: invoke sort and then unique */
void order( PartSet & );

/** Insert a part into a properly ordered collection of parts.
 *  Return if this is a new insertion.
 */
bool insert( PartSet & , Part & );

/** Find a part by name in a properly ordered collection of parts. */
Part * find( const PartSet & , const std::string & );

/** Query containment for properly ordered PartSet */
bool contain( const PartSet & , const Part & );
bool contain( const PartSet & , const PartSet & );

/** Query cardinality of intersection of two PartSets */
unsigned intersect( const PartSet & , const PartSet & );

/** Generate intersection of two PartSets */
unsigned intersect( const PartSet & , const PartSet & , PartSet & );

/** Query if two parts intersect:
 *  If one is a subset of the other or they share a common subset
 */
bool intersect( const Part & , const Part & );

//----------------------------------------------------------------------
/** Verify consistency of supersets, subsets, and partitions */
bool verify( const Part & , std::string & );

/** Print a part.  Each line starts with the given leader string */
std::ostream & print( std::ostream & , const char * const , const Part & );

//----------------------------------------------------------------------
/** A part relation defined as follows:
 *  if Entity 'e1' is a member of part 'm_root' and
 *     there exists a relation from Entity 'e1' to Entity 'e2' that
 *     is in the domain of the relation stencil 'm_function'
 *  then Entity 'e2' is a member of part 'm_target'.
 *
 *  This data structure is used internally and should never need to be
 *  used by a user of the phdBulkData package.
 */
struct PartRelation {
  Part               * m_root ;
  Part               * m_target ;
  relation_stencil_ptr m_function ;

  PartRelation() : m_root( NULL ), m_target( NULL ), m_function( NULL ) {}

  PartRelation( const PartRelation & rhs )
    : m_root( rhs.m_root ),
      m_target( rhs.m_target ),
      m_function( rhs.m_function ) {}

  PartRelation & operator = ( const PartRelation & rhs )
    {
      m_root = rhs.m_root ;
      m_target = rhs.m_target ;
      m_function = rhs.m_function ;
      return *this ;
    }
};

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

