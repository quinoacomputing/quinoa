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
#include <util/TPI.h>
#include <util/Parallel.hpp>
#include <util/ParallelInputStream.hpp>
#include <util/TestDriver.hpp>

//----------------------------------------------------------------------

void test_reduce( phdmesh::ParallelMachine , std::istream & );
void test_accuracy( phdmesh::ParallelMachine, std::istream & );
void test_timing_blas1( phdmesh::ParallelMachine , std::istream & );
void test_timing_mxv( phdmesh::ParallelMachine, std::istream & );

//----------------------------------------------------------------------

int main( int argc , char **argv )
{
  phdmesh::TestDriverMap test_map ;

  test_map[ std::string("reduce") ] = & test_reduce ;
  test_map[ std::string("accuracy") ] = & test_accuracy ;
  test_map[ std::string("timing_mxv") ] = & test_timing_mxv ;
  test_map[ std::string("timing_blas1") ] = & test_timing_blas1 ;

  phdmesh::ParallelMachine comm = phdmesh::parallel_machine_init(&argc,&argv);

  const unsigned p_rank = phdmesh::parallel_machine_rank( comm );

  phdmesh::ParallelInputStream
    is( comm , ( 0 == p_rank && argc == 2 ? argv[1] : NULL ) );

  const int result = phdmesh::test_driver( comm , is , test_map );

  phdmesh::parallel_machine_finalize();

  return result ;
}

