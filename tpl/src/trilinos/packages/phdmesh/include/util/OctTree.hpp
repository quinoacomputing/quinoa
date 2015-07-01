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
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   January 2007
 */

#ifndef pddgeom_OctTree_hpp
#define pddgeom_OctTree_hpp

#include <limits>
#include <iosfwd>
#include <util/Basics.hpp>
#include <util/SimpleArrayOps.hpp>

namespace phdmesh {

class OctTreeKey ;

}

namespace std {

ostream & operator << ( ostream & , const phdmesh::OctTreeKey & );

}

namespace phdmesh {

class OctTreeKey {
public:
  typedef unsigned value_type ;
  enum { MaxDepth     = 16 };
  enum { MaskIndex    = 0x0f };
  enum { BitsPerIndex = 4 };
  enum { BitsPerWord  = std::numeric_limits<value_type>::digits };
  enum { IndexPerWord = BitsPerWord / BitsPerIndex };
  enum { NWord        = MaxDepth / IndexPerWord };
  enum { OKBits = StaticAssert< 0 == BitsPerWord % BitsPerIndex >::OK };
  enum { OKWord = StaticAssert< 0 == MaxDepth    % IndexPerWord >::OK };
private:
  value_type m_value[ NWord ];
public:

  OctTreeKey()
    { Copy<NWord>( m_value , 0u ); }

  OctTreeKey( const OctTreeKey & k )
    { Copy<NWord>( m_value , k.m_value ); }

  OctTreeKey & operator = ( const OctTreeKey & k )
    { Copy<NWord>( m_value , k.m_value ); return *this ; }

  bool operator == ( const OctTreeKey & k ) const
    { return Equal<NWord>( m_value , k.m_value ); }

  bool operator != ( const OctTreeKey & k ) const
    { return Equal<NWord>( m_value , k.m_value ); }

  bool operator < ( const OctTreeKey & k ) const
    { return Less<NWord>( m_value , k.m_value ); }

  /** Query depth of this key */
  unsigned depth() const ;

  /** Index of the key at the depth [1..8]
   *  A zero value indicates it is not defined at that depth.
   */
  template<unsigned Depth> unsigned index() const ;

  /** Index of the key at depth */
  unsigned index( const unsigned Depth ) const ;

  /** Clear index at depth */
  template<unsigned D> OctTreeKey & clear_index() ;

  /** Clear index at depth */
  OctTreeKey & clear_index( const unsigned Depth );

  /** Set index at depth, range is [1..8] */
  template<unsigned D> OctTreeKey & set_index( const unsigned );

  /** Set index at depth */
  OctTreeKey & set_index( const unsigned Depth , const unsigned );

  /** Clear the value (root node) */
  OctTreeKey & clear();

  /** Set to the maximum value */
  OctTreeKey & set_maximum();

  /** Query raw value */
  const value_type * value() const { return m_value ; }

  /** Set raw value */
  OctTreeKey & set_value( const value_type * );

  /** Intersects if either key contains the other. */
  bool intersect( const OctTreeKey & k ) const ;

private:
  void throw_depth( const unsigned ) const ;
  void throw_index( const unsigned , const unsigned ) const ;
};

//----------------------------------------------------------------------

/** Generate a 3D Hilbert space filling curve oct-tree key. */

OctTreeKey hsfc3d( const unsigned Depth , const unsigned * const coord );

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<unsigned Depth> struct OctTreeSize ;

template<>
struct OctTreeSize<0>
{ enum { value = 1 }; };

template<unsigned Depth>
struct OctTreeSize
{
  enum { MaxDepth = 10 , N = Depth }; // Size representable by an unsigned int

  enum { OK = StaticAssert< N <= MaxDepth >::OK };

  enum { value = 1 + 8 * OctTreeSize<Depth-1>::value };
};

unsigned oct_tree_size( const unsigned Depth );

/** Offset of a oct-tree node in a dense tree of a given depth. */
unsigned oct_tree_offset( const unsigned Depth , const OctTreeKey & );

template<unsigned Depth>
inline
unsigned oct_tree_offset( const OctTreeKey & k )
{ return oct_tree_offset( Depth , k ); }

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

template<unsigned Depth>
inline
unsigned OctTreeKey::index() const
{
  enum { OK = StaticAssert< 0 < Depth && Depth <= MaxDepth >::OK };
  enum { which = ( Depth - 1 ) / IndexPerWord };
  enum { shift = BitsPerWord - BitsPerIndex * ( Depth % IndexPerWord ) };

  return ( m_value[ which ] >> shift ) & MaskIndex ;
}

template<unsigned Depth>
inline
OctTreeKey & OctTreeKey::clear_index()
{
  enum { OK = StaticAssert< 0 < Depth && Depth <= MaxDepth >::OK };
  enum { which = ( Depth - 1 ) / IndexPerWord };
  enum { shift = BitsPerWord - BitsPerIndex * ( Depth % IndexPerWord ) };

  const value_type m = MaskIndex ;

  m_value[ which ] &= ~( m << shift );

  return *this ;
}

template<unsigned Depth>
inline
OctTreeKey & OctTreeKey::set_index( const unsigned Index )
{
  enum { OK = StaticAssert< 0 < Depth && Depth <= MaxDepth >::OK };
  enum { which = ( Depth - 1 ) / IndexPerWord };
  enum { shift = BitsPerWord - BitsPerIndex * ( Depth % IndexPerWord ) };

  if ( 8 < Index ) { throw_index( Depth , Index ); }

  const value_type m = MaskIndex ;

  ( m_value[which] &= ~( m << shift ) ) |= Index << shift ;

  return *this ;
}

}

#endif

