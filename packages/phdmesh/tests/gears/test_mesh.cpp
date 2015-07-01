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
#include <stdlib.h>
#include <string.h>
#include <map>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include <util/TestDriver.hpp>
#include <util/Parallel.hpp>
#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <util/ParallelInputStream.hpp>

using namespace phdmesh ;

void test_simple_mesh( ParallelMachine , std::istream & );
void test_gears( ParallelMachine , std::istream & );
void test_schema_parts( ParallelMachine , std::istream & );

//----------------------------------------------------------------------
//----------------------------------------------------------------------

int main( int argc , char ** argv )
{
  TestDriverMap test_map ;

  test_map[ std::string( "schema_parts" ) ] = & test_schema_parts ;
  test_map[ std::string( "simple_mesh" ) ]  = & test_simple_mesh ;
  test_map[ std::string( "gears" ) ]        = & test_gears ;

  //----------------------------------

  phdmesh::ParallelMachine comm = phdmesh::parallel_machine_init(&argc,&argv);

  const unsigned p_rank = phdmesh::parallel_machine_rank( comm );

  int result = -1 ;

  try {
    result = phdmesh::test_driver( comm , test_map , argc , argv );
  }
  catch( const std::exception & x ) {
    std::cout << "P" << p_rank << ": " << x.what() << std::endl ;
  }

  phdmesh::parallel_machine_finalize();

  //----------------------------------

  return result ;
}

