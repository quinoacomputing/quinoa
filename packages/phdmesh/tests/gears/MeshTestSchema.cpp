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

#include <util/TPI.h>
#include <mesh/MetaData.hpp>

using namespace phdmesh ;

namespace {

void require( bool condition , const char * const method ,
                               const char * const detail )
{
  if ( ! condition ) {
    std::string msg ;
    msg.append( method );
    msg.append( " : FAILED " );
    msg.append( detail );
    std::cerr << msg << std::endl ;
    std::cerr.flush();
    throw std::runtime_error( msg );
  }
}

void require_get( MetaData & mesh_meta_data , const Part * part )
{
  static const char method[] = "phdmesh::test::schema_parts::require_get" ;

  require( part == mesh_meta_data.get_part( part->name() ) ,
           method , "part = get_part( part->name()" );
}

}

void test_schema_parts( ParallelMachine comm , std::istream & )
{
  static const char method[] = "phdmesh::test::schema_parts" ;

  if ( 0 < parallel_machine_rank( comm ) ) { return ; }

  MetaData mesh_meta_data ;

  // Should have exactly five predefined parts

  require( mesh_meta_data.get_parts().size() == 3 , method , "predefined parts" );

  const Part & universal = mesh_meta_data.universal_part();
  const Part & uses      = mesh_meta_data.locally_used_part();
  const Part & owns      = mesh_meta_data.locally_owned_part();

  require_get( mesh_meta_data , & universal );
  require_get( mesh_meta_data , & uses );
  require_get( mesh_meta_data , & owns );

  //--------------------------------------------------------------------
  // Parts and subsets

  Part & a     = mesh_meta_data.declare_part( std::string("a") );
  Part & a_a   = mesh_meta_data.declare_part( std::string("a_a") );
  Part & a_b   = mesh_meta_data.declare_part( std::string("a_b") );
  Part & a_b_c = mesh_meta_data.declare_part( std::string("a_b_c") );
  Part & b     = mesh_meta_data.declare_part( std::string("b") );
  Part & b_a   = mesh_meta_data.declare_part( std::string("b_a") );
  Part & b_b   = mesh_meta_data.declare_part( std::string("b_b") );
  Part & b_b_c = mesh_meta_data.declare_part( std::string("b_b_c") );

  mesh_meta_data.declare_part_subset( a , a_a );
  mesh_meta_data.declare_part_subset( a , a_b );
  mesh_meta_data.declare_part_subset( a_b , a_b_c );

  mesh_meta_data.declare_part_subset( b , b_a );
  mesh_meta_data.declare_part_subset( b , b_b );
  mesh_meta_data.declare_part_subset( b_b , b_b_c );

  PartSet x_intersect ;
  {
    Part * tmp ;
    tmp = & a_b ;   x_intersect.push_back( tmp );
    tmp = & a_b_c ; x_intersect.push_back( tmp );
    tmp = & b_b ;   x_intersect.push_back( tmp );
    tmp = & b_b_c ; x_intersect.push_back( tmp );
  }

  Part & x = mesh_meta_data.declare_part( x_intersect );

  Part & y = mesh_meta_data.declare_part( std::string("y") );

  mesh_meta_data.declare_part_subset( a_b_c , y );
  mesh_meta_data.declare_part_subset( b_b_c , y );

  //--------------------------------------------------------------------
  // Test expectations:

/*
  require( contain( universal , a ) ,     method , "universal > a" );
  require( contain( universal , a_b ) ,   method , "universal > a_b" );
  require( contain( universal , a_a ) ,   method , "universal > a_a" );
  require( contain( universal , a_b_c ) , method , "universal > a_b_c" );

  require( contain( universal , b ) ,     method , "universal > b" );
  require( contain( universal , b_b ) ,   method , "universal > b_b" );
  require( contain( universal , b_a ) ,   method , "universal > b_a" );
  require( contain( universal , b_b_c ) , method , "universal > b_b_c" );
*/

  require( contain( a.subsets() , a_a ) ,   method , "a > a_a" );
  require( contain( a.subsets() , a_b ) ,   method , "a > a_b" );
  require( contain( a.subsets() , a_b_c ) , method , "a > a_b_c" );
  require( contain( a_b.subsets() , a_b_c ) , method , "a_b > a_b_c" );

  require( ! intersect( a_a , b_a ) , method , "a_a ^ b_a" );
  require( ! intersect( b_a , a_a ) , method , "b_a ^ a_a" );
  require( ! intersect( a , b_a ) , method , "a ^ b_a" );
  require( ! intersect( b , a_a ) , method , "b ^ a_a" );

  require( contain( a.subsets() , x ) ,     method , "a > x" );
  require( contain( a_b.subsets() , x ) ,   method , "a_b > x" );
  require( contain( a_b_c.subsets() , x ) , method , "a_b_c > x" );
  require( contain( b.subsets() , x ) ,     method , "b > x" );
  require( contain( b_b.subsets() , x ) ,   method , "b_b > x" );
  require( contain( b_b_c.subsets() , x ) , method , "b_b_c > x" );
  require( ! intersect( a_a , x ) ,  method , "a_a ^ x" );
  require( ! intersect( b_a , x ) ,  method , "b_a ^ x" );

  require( contain(x.intersection_of(),a_b_c) , method , "x < a_b_c" );
  require( contain(x.intersection_of(),b_b_c) , method , "x < b_b_c" );
  require( ! contain(x.intersection_of(),a_b) , method , "x ! a_b" );
  require( ! contain(x.intersection_of(),b_b) , method , "x ! b_b" );

  require( contain( x.subsets() , y ) , method , "x > y" );

  //--------------------------------------------------------------------
  // Test errors:

  try {
    mesh_meta_data.declare_part_subset( a_b , a );
    require( false , method , "accepted circular a_b > a" );
  }
  catch(...) {
    ;
  }

  try {
    mesh_meta_data.declare_part_subset( a_b_c , a );
    require( false , method , "accepted circular a_b_c > a" );
  }
  catch(...) {
    ;
  }

  //--------------------------------------------------------------------

  mesh_meta_data.commit();

  std::cout << std::endl ;
  print( std::cout , NULL , a );
  print( std::cout , NULL , a_a );
  print( std::cout , NULL , a_b );
  print( std::cout , NULL , a_b_c );
  print( std::cout , NULL , b );
  print( std::cout , NULL , b_a );
  print( std::cout , NULL , b_b );
  print( std::cout , NULL , b_b_c );
  print( std::cout , NULL , x );
  print( std::cout , NULL , y );
  std::cout << std::endl ;
}

