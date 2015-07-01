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

#ifndef util_TypeName_hpp
#define util_TypeName_hpp

#include <complex>
#include <vector>
#include <string>
#include <typeinfo>

namespace phdmesh {

template< typename T >
struct TypeName {
  static std::string value() { return std::string( typeid(T).name() ); }
};

//----------------------------------------------------------------------

#define GENERATE_SIMPLE_TYPE_NAME_VALUE( T ) \
  template<> struct TypeName< T > { \
    static std::string value() { return std::string( # T ); } \
  }

GENERATE_SIMPLE_TYPE_NAME_VALUE( void );
GENERATE_SIMPLE_TYPE_NAME_VALUE( char );
GENERATE_SIMPLE_TYPE_NAME_VALUE( unsigned char );
GENERATE_SIMPLE_TYPE_NAME_VALUE( short );
GENERATE_SIMPLE_TYPE_NAME_VALUE( unsigned short );
GENERATE_SIMPLE_TYPE_NAME_VALUE( int );
GENERATE_SIMPLE_TYPE_NAME_VALUE( unsigned int );
GENERATE_SIMPLE_TYPE_NAME_VALUE( long );
GENERATE_SIMPLE_TYPE_NAME_VALUE( unsigned long );
GENERATE_SIMPLE_TYPE_NAME_VALUE( float );
GENERATE_SIMPLE_TYPE_NAME_VALUE( double );
GENERATE_SIMPLE_TYPE_NAME_VALUE( std::complex<float> );
GENERATE_SIMPLE_TYPE_NAME_VALUE( std::complex<double> );
GENERATE_SIMPLE_TYPE_NAME_VALUE( std::string );

//----------------------------------------------------------------------

template< typename T >
struct TypeName< const T >
{
  static std::string value()
    { return std::string( "const ").append( TypeName<T>::value() ); }
};

template< typename T >
struct TypeName< T * >
{
  static std::string value()
    { return std::string( TypeName<T>::value() ).append( " *" ); }
};

template< typename T >
struct TypeName< T & >
{
  static std::string value()
    { return std::string( TypeName<T>::value() ).append( " &" ); }
};

//----------------------------------------------------------------------

std::string type_name_array( const std::string & , unsigned );

template< typename T , unsigned N >
struct TypeName< T[N] >
{
  static std::string value()
    { return type_name_array( TypeName<T>::value() , N ); }
};

template< typename T >
std::string type_name_array( unsigned n )
{ return type_name_array( TypeName<T>::value() , n ); }

//----------------------------------------------------------------------

std::string type_name_vector( const std::string & , unsigned = 0 );

template< typename T >
struct TypeName< std::vector<T> >
{
  static std::string value()
    { return type_name_vector( TypeName<T>::value() ); }
};

template< typename T >
std::string type_name_vector( unsigned n = 0 )
{ return type_name_vector( TypeName<T>::value() , n ); }

//----------------------------------------------------------------------

} // namespace phdmesh

#endif /* util_TypeName_hpp */

