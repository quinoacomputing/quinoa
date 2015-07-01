/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
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

#define ARRAY_BOUNDS_CHECKING

#include <iostream>
#include <util/Parallel.hpp>
#include <util/Array.hpp>

namespace {

//----------------------------------------------------------------------

struct TagA : public phdmesh::ArrayDimTag {
  TagA() {}
  const char * name() const ;
  static const TagA & tag();
};

const char * TagA::name() const { static const char n[] = "A" ; return n ; }

const TagA & TagA::tag() { static const TagA self ; return self ; }

//----------------------------------------------------------------------

struct TagB : public phdmesh::ArrayDimTag {
  TagB() {}
  const char * name() const ;
  static const TagB & tag();
};

const char * TagB::name() const { static const char n[] = "B" ; return n ; }

const TagB & TagB::tag() { static const TagB self ; return self ; }

//----------------------------------------------------------------------

struct TagC : public phdmesh::ArrayDimTag {
  TagC() {}
  const char * name() const ;
  static const TagC & tag();
};

const char * TagC::name() const { static const char n[] = "C" ; return n ; }

const TagC & TagC::tag() { static const TagC self ; return self ; }

//----------------------------------------------------------------------

struct TagD : public phdmesh::ArrayDimTag {
  TagD() {}
  const char * name() const ;
  static const TagD & tag();
};

const char * TagD::name() const { static const char n[] = "D" ; return n ; }

const TagD & TagD::tag() { static const TagD self ; return self ; }

//----------------------------------------------------------------------

using namespace phdmesh ;

void myfortranfunc( const Array<double,FortranOrder> xf )
{
  std::cout << "myfortranfunc( Array<double,FortranOrder" ;

  if ( xf.rank() && NULL != xf.tag(0) ) {
    for ( unsigned i = 0 ; i < xf.rank() ; ++i ) {
      std::cout << "," << xf.tag(i)->name();
    }
  }

  std::cout << ">(" ;
  std::cout << (void*) xf.contiguous_data();
  for ( unsigned i = 0 ; i < xf.rank() ; ++i ) {
    std::cout << "," << xf.dimension(i);
  }
  std::cout << ") )" << std::endl ;
}

void mynaturalfunc( const Array<double,NaturalOrder> xf )
{
  std::cout << "mynaturalfunc( Array<double,NaturalOrder" ;

  if ( xf.rank() && NULL != xf.tag(0) ) {
    for ( unsigned i = 0 ; i < xf.rank() ; ++i ) {
      std::cout << "," << xf.tag(i)->name();
    }
  }

  std::cout << ">(" ;
  std::cout << (void*) xf.contiguous_data();
  for ( unsigned i = 0 ; i < xf.rank() ; ++i ) {
    std::cout << "," << xf.dimension(i);
  }
  std::cout << ") )" << std::endl ;
}

//----------------------------------------------------------------------

void myfortranA( const Array<double,FortranOrder,TagA> )
{}

void myfortranAB( const Array<double,FortranOrder,TagA,TagB> )
{}

void myfortranABC( const Array<double,FortranOrder,TagA,TagB,TagC> )
{}

void myfortranABCD( const Array<double,FortranOrder,TagA,TagB,TagC,TagD> )
{}

//----------------------------------------------------------------------

void local_test_array()
{
  double storage[100000];

  Array<double,FortranOrder,TagA>
    af1( storage , 2 );

  Array<double,FortranOrder,TagA,TagB>
    af2( storage , 2 , 3 );

  Array<double,FortranOrder,TagA,TagB,TagC>
    af3( storage , 2 , 3 , 4 );

  Array<double,FortranOrder,TagA,TagB,TagC,TagD>
    af4( storage , 2 , 3 , 4 , 5 );

  Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA>
    af5( storage , 2 , 3 , 4 , 5 , 6 );

  Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB>
    af6( storage , 2 , 3 , 4 , 5 , 6 , 7 );

  Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC>
    af7( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 );

  Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC,TagD>
    af8( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 );

  Array<double,NaturalOrder,TagA>
    an1( storage , 2 );

  Array<double,NaturalOrder,TagA,TagB>
    an2( storage , 2 , 3 );

  Array<double,NaturalOrder,TagA,TagB,TagC>
    an3( storage , 2 , 3 , 4 );

  Array<double,NaturalOrder,TagA,TagB,TagC,TagD>
    an4( storage , 2 , 3 , 4 , 5 );

  Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA>
    an5( storage , 2 , 3 , 4 , 5 , 6 );

  Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA,TagB>
    an6( storage , 2 , 3 , 4 , 5 , 6 , 7 );

  Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC>
    an7( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 );

  Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC,TagD>
    an8( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 );

  std::cout << std::endl << "FORTRAN ARRAYS:" << std::endl ;

  myfortranfunc( af1 );
  myfortranfunc( af2 );
  myfortranfunc( af3 );
  myfortranfunc( af4 );
  myfortranfunc( af5 );
  myfortranfunc( af6 );
  myfortranfunc( af7 );
  myfortranfunc( af8 );

  mynaturalfunc( af1 );
  mynaturalfunc( af2 );
  mynaturalfunc( af3 );
  mynaturalfunc( af4 );
  mynaturalfunc( af5 );
  mynaturalfunc( af6 );
  mynaturalfunc( af7 );
  mynaturalfunc( af8 );

  myfortranfunc( af8.truncate(0) );

  std::cout << std::endl << "NATURAL ARRAYS:" << std::endl ;

  mynaturalfunc( an1 );
  mynaturalfunc( an2 );
  mynaturalfunc( an3 );
  mynaturalfunc( an4 );
  mynaturalfunc( an5 );
  mynaturalfunc( an6 );
  mynaturalfunc( an7 );
  mynaturalfunc( an8 );

  myfortranfunc( an1 );
  myfortranfunc( an2 );
  myfortranfunc( an3 );
  myfortranfunc( an4 );
  myfortranfunc( an5 );
  myfortranfunc( an6 );
  myfortranfunc( an7 );
  myfortranfunc( an8 );

  mynaturalfunc( an8.truncate(0) );

  //------------------------------

  myfortranA( af1 );
  myfortranA( an1 ); // Implicit conversion-construction is good

  // myfortranA( af2 ); // Compile error catches correctly
  // myfortranA( af3 ); // Compile error catches correctly
  // myfortranA( af4 ); // Compile error catches correctly

  myfortranAB( af2 );

  // myfortranAB( an2 ); // Compile error catches correctly
  // myfortranAB( af3 ); // Compile error catches correctly
  // myfortranAB( af4 ); // Compile error catches correctly

  myfortranABC( af3 );
  myfortranABCD( af4 );

  //------------------------------

  try {
    for ( unsigned i = 0 ; i < 100 ; ++i ) { af1(i); }
  }
  catch( const std::exception & x ) {
    std::cout << "Correctly caught exception:"
              << std::endl << x.what() << std::endl ;
  }
}

}

void test_array( phdmesh::ParallelMachine , std::istream & )
{
  try {
    local_test_array();
    std::cout << "TEST_ARRAY PASSED" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "TEST_ARRAY FAILED: " << x.what() << std::endl ;
  }
  catch( ... ) {
    std::cout << "TEST_ARRAY FAILED: <unknown>" << std::endl ;
  }
}

