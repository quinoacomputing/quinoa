/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ElemSampleQMTest.cpp
 *  \brief Unit tests for ElementSampleQM class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ElemSampleQM.hpp"
#include <cppunit/extensions/HelperMacros.h>

using namespace Mesquite;

class ElemSampleQMTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(ElemSampleQMTest);
  CPPUNIT_TEST (test_handle_from_sample);
  CPPUNIT_TEST_SUITE_END();

public:
  void test_handle_from_sample();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ElemSampleQMTest, "ElemSampleQMTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ElemSampleQMTest, "Unit");

void test_handle( Sample sample, size_t elem )
{
  size_t handle = ElemSampleQM::handle( sample, elem );
  CPPUNIT_ASSERT_EQUAL( sample.dimension, ElemSampleQM::sample( handle ).dimension );
  CPPUNIT_ASSERT_EQUAL( sample.number,    ElemSampleQM::sample( handle ).number );
  CPPUNIT_ASSERT_EQUAL( elem,             (size_t)ElemSampleQM::elem( handle ) );
}

void ElemSampleQMTest::test_handle_from_sample( )
{
  Sample min_sample_no(0,0);
  Sample max_sample_no( Sample::SIDE_DIMENSON_MASK, Sample::SIDE_NUMBER_MASK );
  size_t min_elem_no = 0;
  size_t max_elem_no = ElemSampleQM::MAX_ELEM_PER_PATCH - 1;
  
  test_handle( min_sample_no, min_elem_no );
  test_handle( min_sample_no, max_elem_no );
  test_handle( max_sample_no, min_elem_no );
  test_handle( max_sample_no, max_elem_no );
}
