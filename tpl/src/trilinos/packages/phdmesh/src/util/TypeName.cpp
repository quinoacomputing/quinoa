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
 * @date   October 2007
 */

#include <sstream>
#include <util/TypeName.hpp>

namespace phdmesh {

std::string type_name_array( const std::string & type , unsigned n )
{
  std::ostringstream tmp ;
  tmp << type << "[" << n << "]" ;
  return tmp.str();
}

std::string type_name_vector( const std::string & type , unsigned n )
{
  std::ostringstream tmp ;
  tmp << "std::vector< " << type << " >" ;
  if ( n ) { tmp << "(" << n << ")" ; }
  return tmp.str();
}

}

#if 0

#define TYPE_NAME_VALUE( T ) \
const char * TypeName< T >::value() { static const char v[] = # T ; return v ; }

namespace phdmesh {

TYPE_NAME_VALUE( void )
TYPE_NAME_VALUE( char )
TYPE_NAME_VALUE( unsigned char )
TYPE_NAME_VALUE( short )
TYPE_NAME_VALUE( unsigned short )
TYPE_NAME_VALUE( int )
TYPE_NAME_VALUE( unsigned int )
TYPE_NAME_VALUE( long )
TYPE_NAME_VALUE( unsigned long )

TYPE_NAME_VALUE( float )
TYPE_NAME_VALUE( double )
TYPE_NAME_VALUE( std::complex<float> )
TYPE_NAME_VALUE( std::complex<double> )

TYPE_NAME_VALUE( void * )
TYPE_NAME_VALUE( char * )
TYPE_NAME_VALUE( unsigned char * )
TYPE_NAME_VALUE( short * )
TYPE_NAME_VALUE( unsigned short * )
TYPE_NAME_VALUE( int * )
TYPE_NAME_VALUE( unsigned int * )
TYPE_NAME_VALUE( long * )
TYPE_NAME_VALUE( unsigned long * )

TYPE_NAME_VALUE( float * )
TYPE_NAME_VALUE( double * )
TYPE_NAME_VALUE( std::complex<float> * )
TYPE_NAME_VALUE( std::complex<double> * )

TYPE_NAME_VALUE( std::string )

} // namespace phdmesh

#endif

