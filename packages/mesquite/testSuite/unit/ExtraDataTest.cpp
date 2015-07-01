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


/** \file ExtraDataTest.cpp
 *  \brief Test the ExtraData functionality of the PatchData class
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "PatchDataInstances.hpp"
#include "ExtraData.hpp"

using namespace Mesquite;

#include <iostream>

class ExtraDataTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(ExtraDataTest);
  CPPUNIT_TEST (test_initialize);
  CPPUNIT_TEST (test_finalize);
  CPPUNIT_TEST (test_notify_destroyed);
  CPPUNIT_TEST (test_notify_subpatch);
  CPPUNIT_TEST (test_notify_new_patch_fill);
  CPPUNIT_TEST (test_notify_new_patch_sub);
  CPPUNIT_TEST (test_multiple_data);
  CPPUNIT_TEST_SUITE_END();

public:
  void test_initialize();
  void test_finalize();
  void test_notify_destroyed();
  void test_notify_subpatch();
  void test_notify_new_patch_fill();
  void test_notify_new_patch_sub();
  void test_multiple_data();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ExtraDataTest, "ExtraDataTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(ExtraDataTest, "Unit");

class TestExtraData : public Mesquite::ExtraData
{
  public:
    enum EventType {
      NONE = 0,
      DESTROYED,
      NEW,
      SUB
    };
    
  private:
    EventType lastEvent;
    PatchData* subPatchPtr;
  
  public:
    TestExtraData( PatchData& patch )
      : ExtraData( patch ) { clear(); }
    
    void clear()
      { lastEvent = NONE; subPatchPtr = 0; }
    EventType get_last_event( ) const
      { return lastEvent; }
    PatchData* get_sub_patch_ptr() const
      { return subPatchPtr; }
    
    virtual void notify_patch_destroyed();
    virtual void notify_new_patch( );
    virtual void notify_sub_patch( PatchData& sub_patch, 
                                   const size_t* vertex_map,
                                   const size_t* element_map,
                                   MsqError& err );
                                   
    std::vector<size_t> vertexMap, elementMap;
};

void TestExtraData::notify_patch_destroyed() { lastEvent = DESTROYED; }
void TestExtraData::notify_new_patch( ) { lastEvent = NEW; }
void TestExtraData::notify_sub_patch( PatchData& p, 
                                      const size_t* vertex_map,
                                      const size_t* element_map,
                                      MsqError& ) 
{ 
  lastEvent = SUB;
  subPatchPtr = &p;
  vertexMap.resize( p.num_nodes() );
  elementMap.resize( p.num_elements() );
  std::copy( vertex_map, vertex_map+p.num_nodes(), vertexMap.begin() );
  std::copy( element_map, element_map+p.num_elements(), elementMap.begin() );
}

void ExtraDataTest::test_initialize()
{
  PatchData patch;
  TestExtraData data1( patch );
  TestExtraData data2( patch );
  CPPUNIT_ASSERT_EQUAL( data1.get_patch_data(), &patch );
  CPPUNIT_ASSERT_EQUAL( data2.get_patch_data(), &patch );
  CPPUNIT_ASSERT_EQUAL( data1.get_last_event(), TestExtraData::NONE );
  CPPUNIT_ASSERT_EQUAL( data2.get_last_event(), TestExtraData::NONE );
}

void ExtraDataTest::test_finalize()
{
  PatchData patch;
  
  // PatchData doesn't expose any public method to
  // verify that the ExtraData object was removed when it
  // was destroyed!.  Do some funky stuff to be able to
  // call the ExtraData destructor w/out invalidating the
  // memory of the actual object.
  unsigned char* mem = new unsigned char[sizeof(TestExtraData)];
  TestExtraData* data = new (mem) TestExtraData( patch );
  data->~TestExtraData();
  CPPUNIT_ASSERT_EQUAL( data->get_patch_data(), (PatchData*)0 );
  delete [] mem;
}

void ExtraDataTest::test_notify_destroyed()
{
  PatchData* patch = new PatchData;
  TestExtraData data1(*patch);
  TestExtraData data2(*patch);

  delete patch;
  CPPUNIT_ASSERT_EQUAL( data1.get_last_event(), TestExtraData::DESTROYED );
  CPPUNIT_ASSERT_EQUAL( data2.get_last_event(), TestExtraData::DESTROYED );
  CPPUNIT_ASSERT_EQUAL( data1.get_patch_data(), (PatchData*)0 );
  CPPUNIT_ASSERT_EQUAL( data2.get_patch_data(), (PatchData*)0 );
}

void ExtraDataTest::test_notify_subpatch()
{
  MsqPrintError err(std::cerr);
  PatchData patch, subpatch;
  create_four_quads_patch( patch, err ); CPPUNIT_ASSERT(!err);
  
  TestExtraData data1(patch);
  TestExtraData data2(patch);
  patch.get_subpatch( 0, 1, subpatch, err ); CPPUNIT_ASSERT(!err);
  
  CPPUNIT_ASSERT_EQUAL( data1.get_last_event(), TestExtraData::SUB );
  CPPUNIT_ASSERT_EQUAL( data1.get_sub_patch_ptr(), &subpatch );
  CPPUNIT_ASSERT_EQUAL( data1.get_patch_data(), &patch );
  
  size_t i;
  for (i = 0; i < subpatch.num_nodes(); ++i)
    CPPUNIT_ASSERT_EQUAL( subpatch.get_vertex_handles_array()[i],
                          patch.get_vertex_handles_array()[data1.vertexMap[i]] );
  for (i = 0; i < subpatch.num_elements(); ++i)
    CPPUNIT_ASSERT_EQUAL( subpatch.get_element_handles_array()[i],
                          patch.get_element_handles_array()[data1.elementMap[i]] );
  
  CPPUNIT_ASSERT_EQUAL( data2.get_last_event(), TestExtraData::SUB );
  CPPUNIT_ASSERT_EQUAL( data2.get_sub_patch_ptr(), &subpatch );
  CPPUNIT_ASSERT_EQUAL( data2.get_patch_data(), &patch );
  
  for (i = 0; i < subpatch.num_nodes(); ++i)
    CPPUNIT_ASSERT_EQUAL( subpatch.get_vertex_handles_array()[i],
                          patch.get_vertex_handles_array()[data2.vertexMap[i]] );
  for (i = 0; i < subpatch.num_elements(); ++i)
    CPPUNIT_ASSERT_EQUAL( subpatch.get_element_handles_array()[i],
                          patch.get_element_handles_array()[data2.elementMap[i]] );
}

void ExtraDataTest::test_notify_new_patch_fill()
{
  MsqPrintError err(std::cerr);
  PatchData patch;
  TestExtraData data1(patch);
  TestExtraData data2(patch);
  
  create_four_quads_patch( patch, err ); CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( data1.get_last_event(), TestExtraData::NEW );
  CPPUNIT_ASSERT_EQUAL( data2.get_last_event(), TestExtraData::NEW );
  
  data1.clear();
  data2.clear();
  
  create_four_quads_patch( patch, err ); CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( data1.get_last_event(), TestExtraData::NEW );
  CPPUNIT_ASSERT_EQUAL( data2.get_last_event(), TestExtraData::NEW );
}


void ExtraDataTest::test_notify_new_patch_sub()
{
  MsqPrintError err(std::cerr);
  PatchData patch, subpatch;
  TestExtraData data1(subpatch);
  TestExtraData data2(subpatch);
  
  create_four_quads_patch( patch, err ); CPPUNIT_ASSERT(!err);
  patch.get_subpatch( 0, 1, subpatch, err ); CPPUNIT_ASSERT(!err);
  CPPUNIT_ASSERT_EQUAL( data1.get_last_event(), TestExtraData::NEW );
  CPPUNIT_ASSERT_EQUAL( data2.get_last_event(), TestExtraData::NEW );
}


void ExtraDataTest::test_multiple_data()
{
  PatchData* patch;
  TestExtraData *data1, *data2, *data3, *data4;
  
  patch = new PatchData;
  data1 = new TestExtraData( *patch );
  data2 = new TestExtraData( *patch );
  data3 = new TestExtraData( *patch );
  data4 = new TestExtraData( *patch );
  
  delete data2;
  delete data4;
  CPPUNIT_ASSERT_EQUAL( data1->get_patch_data(), patch );
  CPPUNIT_ASSERT_EQUAL( data3->get_patch_data(), patch );
  
  delete patch;
  CPPUNIT_ASSERT_EQUAL( data1->get_last_event(), TestExtraData::DESTROYED );
  CPPUNIT_ASSERT_EQUAL( data3->get_last_event(), TestExtraData::DESTROYED );
  CPPUNIT_ASSERT_EQUAL( data1->get_patch_data(), (PatchData*)0 );
  CPPUNIT_ASSERT_EQUAL( data3->get_patch_data(), (PatchData*)0 );

  delete data1;
  delete data3;
}

