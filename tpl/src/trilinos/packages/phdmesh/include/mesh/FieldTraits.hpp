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

#ifndef phdmesh_FieldTraits_hpp
#define phdmesh_FieldTraits_hpp

//----------------------------------------------------------------------

#include <util/Array.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** Dimension traits for cartesian vector. */

struct Cartesian : public ArrayDimTag {

  /** Default to 3D */
  enum { Size = 3 };

  enum { X = 0 , Y = 1 , Z = 2 };

  const char * name() const ;
  std::string to_string( unsigned size , unsigned index ) const ;
  unsigned    to_index(  unsigned size , const std::string & ) const ;
  static const Cartesian & tag();

protected:
  Cartesian() {}
private:
  Cartesian( const Cartesian & );
  Cartesian & operator = ( const Cartesian & );
};

//----------------------------------------------------------------------
/** Dimension traits for cylindrical vector */

struct Cylindrical : public ArrayDimTag {

  enum { Radius = 0 , R = 0 ,
         Angle = 1 ,  A = 1 ,
         Z = 2 };

  const char * name() const ;
  std::string to_string( unsigned size , unsigned index ) const ;
  unsigned    to_index(  unsigned size , const std::string & ) const ;
  static const Cylindrical & tag();

private:
  Cylindrical() {}
  Cylindrical( const Cylindrical & );
  Cylindrical & operator = ( const Cylindrical & );
};

//----------------------------------------------------------------------

struct FullTensor : public ArrayDimTag {

  enum { Size = 9 };
  enum { XX = 0 , XY = 3 , XZ = 6 ,
         YX = 1 , YY = 4 , YZ = 7 ,
         ZX = 2 , ZY = 5 , ZZ = 8 };

  const char * name() const ;
  std::string to_string( unsigned, unsigned) const  ;
  unsigned    to_index(  unsigned, const std::string & ) const  ;
  static const FullTensor & tag();

private:
  FullTensor() {}
  FullTensor( const FullTensor & );
  FullTensor & operator = ( const FullTensor & );
};

//----------------------------------------------------------------------

struct SymmetricTensor : public ArrayDimTag {

  enum { Size = 6 };
  enum { XX = 0 , YY = 1 , ZZ = 2, XY = 3, YZ = 4, XZ = 5 };

  const char * name() const  ;
  std::string to_string( unsigned, unsigned) const  ;
  unsigned    to_index(  unsigned , const std::string & ) const ;
  static const SymmetricTensor & tag();

private:
  SymmetricTensor() {}
  SymmetricTensor( const SymmetricTensor & );
  SymmetricTensor & operator = ( const SymmetricTensor & );
};

//----------------------------------------------------------------------
/** Dimension traits for array of entities */

struct EntityDimension : public ArrayDimTag {

  const char * name() const ;

  static const EntityDimension & tag();

private:
  EntityDimension() {}
  EntityDimension( const EntityDimension & );
  EntityDimension & operator = ( const EntityDimension & );
};

//----------------------------------------------------------------------

}

#endif

