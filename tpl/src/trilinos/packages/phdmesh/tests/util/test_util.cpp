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

#include <iostream>
#include <util/Parallel.hpp>
#include <util/ParallelInputStream.hpp>
#include <util/TestDriver.hpp>

//----------------------------------------------------------------------

void test_array( phdmesh::ParallelMachine , std::istream & );
void test_containers( phdmesh::ParallelMachine , std::istream & );
void test_comm_bounds( phdmesh::ParallelMachine , std::istream & );
void test_comm_sparse( phdmesh::ParallelMachine , std::istream & );
void test_comm_dense(  phdmesh::ParallelMachine , std::istream & );
void test_comm_all(  phdmesh::ParallelMachine , std::istream & );
void test_global_box( phdmesh::ParallelMachine comm , std::istream & );
void test_oct_tree( phdmesh::ParallelMachine , std::istream & );
void test_oct_tree_part_course( phdmesh::ParallelMachine , std::istream & );

void test_oct_tree_comm_part( phdmesh::ParallelMachine , std::istream & );

void test_oct_tree_global_search( phdmesh::ParallelMachine , std::istream & );

void test_oct_tree_global_search_time(
  phdmesh::ParallelMachine , std::istream & );

//----------------------------------------------------------------------

int main( int argc , char **argv )
{
  phdmesh::TestDriverMap test_map ;

  test_map[ std::string("array") ] = & test_array ;
  test_map[ std::string("containers") ] = & test_containers ;
  test_map[ std::string("bounds") ] = & test_comm_bounds ;
  test_map[ std::string("sparse") ] = & test_comm_sparse ;
  test_map[ std::string("comm_all") ] = & test_comm_all ;
  test_map[ std::string("dense") ] = & test_comm_dense ;
  test_map[ std::string("global_box") ] = & test_global_box ;
  test_map[ std::string("oct_tree") ] = & test_oct_tree ;
  test_map[ std::string("oct_tree_part_course") ] = & test_oct_tree_part_course ;
  test_map[ std::string("oct_tree_comm_part") ] = & test_oct_tree_comm_part ;
  test_map[ std::string("oct_tree_global_search") ] = & test_oct_tree_global_search ;
  test_map[ std::string("oct_tree_global_search_time") ] = & test_oct_tree_global_search_time ;

  phdmesh::ParallelMachine comm = phdmesh::parallel_machine_init(&argc,&argv);

  const unsigned p_rank = phdmesh::parallel_machine_rank( comm );

  /* Last command line argument to ignore the 'sierra' script junk */

  const char * const file_name = 1 < argc ? argv[argc-1] : NULL ;

  phdmesh::ParallelInputStream
    is( comm , ( 0 == p_rank ? file_name : NULL ) );

  const int result = phdmesh::test_driver( comm , is , test_map );

  phdmesh::parallel_machine_finalize();

  return result ;
}

