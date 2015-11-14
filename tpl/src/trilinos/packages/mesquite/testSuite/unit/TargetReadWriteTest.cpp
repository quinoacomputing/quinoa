/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retian certain rights to this software.

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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TargetReadWriteTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_TargetWriter.hpp"
#include "Mesquite_TargetReader.hpp"
#include "Mesquite_WeightReader.hpp"
#include "Mesquite_MeshImpl.hpp"
#include "Mesquite_PatchData.hpp"
#include "Mesquite_Settings.hpp"
#include "Mesquite_ElemSampleQM.hpp"

#include "Mesquite_cppunit/extensions/HelperMacros.h"
#include "UnitUtil.hpp"

#include <iostream>

using namespace Mesquite;

class TargetReadWriteTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE( TargetReadWriteTest );
  CPPUNIT_TEST( read_write_3D_targets );
  CPPUNIT_TEST( read_write_2D_targets );
  CPPUNIT_TEST( read_write_surface_targets );
  CPPUNIT_TEST( read_write_weights );
  CPPUNIT_TEST_SUITE_END();
  
  MeshImpl myMesh;   // mesh data
  PatchData myPatch; // global patch for mesh data
  Settings linearMaps; 
public:
  
  void setUp();
  void tearDown();
  void read_write_3D_targets();
  void read_write_2D_targets();
  void read_write_surface_targets();
  void read_write_weights();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetReadWriteTest, "TargetReadWriteTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(TargetReadWriteTest, "Unit");

static const char vtk_file_data[] =
"# vtk DataFile Version 2.0\n"
"Mesquite Mesh\n"
"ASCII\n"
"DATASET UNSTRUCTURED_GRID\n"
"POINTS 11 float\n"
"-1 -1 -1\n"
" 1 -1 -1\n"
" 1  1 -1\n"
"-1  1 -1\n"
"-1 -1  1\n"
" 1 -1  1\n"
" 1  1  1\n"
"-1  1  1\n"
" 0  0  2\n"
"-2  0  2\n"
"-2  0 -1\n"
"CELLS 6 36\n"
"3 0 1 10\n"
"4 0 1 2 3\n"
"4 4 5 8 9\n"
"5 4 5 6 7 8\n"
"6 4 5 9 1 0 10\n"
"8 0 1 2 3 4 5 6 7\n"
"CELL_TYPES 6\n"
"5 9 10 14 13 12\n"
"\n";


void TargetReadWriteTest::setUp()
{
    // create input file
  const char filename[] = "target_reader_test.vtk";
  FILE* file = fopen( filename, "w" );
  CPPUNIT_ASSERT(file);
  int rval = fputs( vtk_file_data, file );
  fclose(file);
  CPPUNIT_ASSERT( rval != EOF );
  
    // read input file
  MsqError err;
  myMesh.read_vtk( filename, err );
  remove(filename);
  if (err) std::cout << err << std::endl;
  CPPUNIT_ASSERT( !err );
  
    // Construct global patch
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::VertexHandle> verts;
  myMesh.get_all_elements( elems, err ); CPPUNIT_ASSERT(!err);
  myMesh.get_all_vertices( verts, err ); CPPUNIT_ASSERT(!err);
  myPatch.set_mesh( &myMesh );
  myPatch.attach_settings( &linearMaps );
  myPatch.set_mesh_entities( elems, verts, err );
  CPPUNIT_ASSERT(!err);
}

void TargetReadWriteTest::tearDown()
{
  myMesh.clear();
}

class FakeTargetCalc : public TargetCalculator, public WeightCalculator
{
  bool surfOrient;

public:
  FakeTargetCalc(bool surface_orient = true) : surfOrient(surface_orient) {}

  ~FakeTargetCalc() {}
  
  bool get_3D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,3>& W_out,
                      MsqError& err );

  bool get_2D_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<2,2>& W_out,
                      MsqError& err );

  bool get_surface_target( PatchData& pd, 
                      size_t element,
                      Sample sample,
                      MsqMatrix<3,2>& W_out,
                      MsqError& err );

  double get_weight( PatchData& pd, 
                     size_t element,
                     Sample sample,
                     MsqError& err );
  
  bool have_surface_orient() const { return surfOrient; }
                     
  unsigned long make_value( Mesh::ElementHandle elem, Sample sample, unsigned idx );
                                   
};

bool FakeTargetCalc::get_3D_target( PatchData& pd, size_t elem, 
                                    Sample sample,
                                    MsqMatrix<3,3>& W_out, MsqError& )
{
  CPPUNIT_ASSERT_EQUAL( 3u, TopologyInfo::dimension( pd.element_by_index(elem).get_element_type() ) );
  unsigned i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < i; ++j)
      W_out(i,j) = 0.0;
    for (j = i; j < 3; ++j)
      W_out(i,j) = make_value( pd.get_element_handles_array()[elem], sample, 3*i+j+1 );
  }
  return true;
}

bool FakeTargetCalc::get_surface_target( PatchData& pd, size_t elem, 
                                    Sample sample,
                                    MsqMatrix<3,2>& W_out, MsqError& )
{
  CPPUNIT_ASSERT_EQUAL( 2u, TopologyInfo::dimension( pd.element_by_index(elem).get_element_type() ) );
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 2; ++j)
      W_out(i,j) = make_value( pd.get_element_handles_array()[elem], sample, 2*i+j );
  return true;
}

bool FakeTargetCalc::get_2D_target( PatchData& pd, size_t elem, 
                                    Sample sample,
                                    MsqMatrix<2,2>& W_out, MsqError& )
{
  CPPUNIT_ASSERT_EQUAL( 2u, TopologyInfo::dimension( pd.element_by_index(elem).get_element_type() ) );
  for (unsigned i = 0; i < 2; ++i)
    for (unsigned j = 0; j < 2; ++j)
      W_out(i,j) = make_value( pd.get_element_handles_array()[elem], sample, (2-i)*(2-j) );
  return true;
}

double FakeTargetCalc::get_weight( PatchData& pd, size_t elem, 
                                   Sample sample,
                                   MsqError& )
{
  return make_value( pd.get_element_handles_array()[elem], sample, 0 );
}

unsigned long FakeTargetCalc::make_value( Mesh::ElementHandle elem, Sample sample, unsigned idx )
{
  const unsigned index_bits = 4;
  CPPUNIT_ASSERT( idx < (1<<index_bits) );
  unsigned long result = (unsigned long)elem;
  result = (result << Sample::SIDE_DIMENSION_BITS) | sample.dimension;
  result = (result << Sample::SIDE_NUMBER_BITS) | sample.number;
  result = (result << index_bits ) | idx;
  return result;
}

void TargetReadWriteTest::read_write_3D_targets()
{
  const bool oriented = true; // doesn't matter which value for 3D

  MsqPrintError err( std::cout );
  FakeTargetCalc tc(oriented);
  
    // Write the targets
  TargetWriter writer( &tc );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, 0);
  writer.loop_over_mesh( &mesh_and_domain, &linearMaps, err );
  CPPUNIT_ASSERT(!err);
  
    // Compare all target matrices
  bool checked_something = false; // make sure mesh actually contains volume elements
  TargetReader reader(oriented);
  for (size_t i = 0; i < myPatch.num_elements(); ++i) {
    const unsigned d = TopologyInfo::dimension( myPatch.element_by_index(i).get_element_type() );
    if (d != 3)
      continue;
    
    checked_something = true;
    std::vector<Sample> samples;
    myPatch.get_samples( i, samples, err ); ASSERT_NO_ERROR(err);
    for (size_t j = 0; j < samples.size(); ++j) {
      MsqMatrix<3,3> expected, read;
      tc.get_3D_target( myPatch, i, samples[j], expected, err );
      CPPUNIT_ASSERT(!err);
      reader.get_3D_target( myPatch, i, samples[j], read, err );
      CPPUNIT_ASSERT(!err);
      ASSERT_MATRICES_EQUAL( expected, read, 1e-12 );
    }
  }
  
  CPPUNIT_ASSERT(checked_something);
}

void TargetReadWriteTest::read_write_surface_targets()
{
  const bool oriented = true;

  MsqPrintError err( std::cout );
  FakeTargetCalc tc(oriented);
  
    // Write the targets
  TargetWriter writer( &tc );
    MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, 0);
  writer.loop_over_mesh( &mesh_and_domain, &linearMaps, err );
  CPPUNIT_ASSERT(!err);
  
    // Compare all target matrices
  bool checked_something = false; // make sure mesh actually contains surface elements
  TargetReader reader(oriented);
  for (size_t i = 0; i < myPatch.num_elements(); ++i) {
    const unsigned d = TopologyInfo::dimension( myPatch.element_by_index(i).get_element_type() );
    if (d != 2)
      continue;
    
    checked_something = true;
    std::vector<Sample> samples;
    myPatch.get_samples( i, samples, err ); ASSERT_NO_ERROR(err);
    for (size_t j = 0; j < samples.size(); ++j) {
      MsqMatrix<3,2> expected, read;
      tc.get_surface_target( myPatch, i, samples[j], expected, err );
      CPPUNIT_ASSERT(!err);
      reader.get_surface_target( myPatch, i, samples[j], read, err );
      CPPUNIT_ASSERT(!err);
      ASSERT_MATRICES_EQUAL( expected, read, 1e-6 );
    }
  }
  
  CPPUNIT_ASSERT(checked_something);
}

void TargetReadWriteTest::read_write_2D_targets()
{
  const bool oriented = false;

  MsqPrintError err( std::cout );
  FakeTargetCalc tc(oriented);
  
    // Write the targets
  TargetWriter writer( &tc );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, 0);
  writer.loop_over_mesh( &mesh_and_domain, &linearMaps, err );
  CPPUNIT_ASSERT(!err);
  
    // Compare all target matrices
  bool checked_something = false; // make sure mesh actually contains surface elements
  TargetReader reader(oriented);
  for (size_t i = 0; i < myPatch.num_elements(); ++i) {
    const unsigned d = TopologyInfo::dimension( myPatch.element_by_index(i).get_element_type() );
    if (d != 2)
      continue;
    
    checked_something = true;
    std::vector<Sample> samples;
    myPatch.get_samples( i, samples, err ); ASSERT_NO_ERROR(err);
    for (size_t j = 0; j < samples.size(); ++j) {
      MsqMatrix<2,2> expected, read;
      tc.get_2D_target( myPatch, i, samples[j], expected, err );
      CPPUNIT_ASSERT(!err);
      reader.get_2D_target( myPatch, i, samples[j], read, err );
      CPPUNIT_ASSERT(!err);
      ASSERT_MATRICES_EQUAL( expected, read, 1e-6 );
    }
  }
  
  CPPUNIT_ASSERT(checked_something);
}

void TargetReadWriteTest::read_write_weights()
{
  MsqPrintError err( std::cout );
  FakeTargetCalc tc;
  
    // Write the targets
  TargetWriter writer( 0, &tc );
  MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&myMesh, 0);
  writer.loop_over_mesh( &mesh_and_domain, &linearMaps, err );
  CPPUNIT_ASSERT(!err);
  
    // Compare all target matrices
  WeightReader reader;
  for (size_t i = 0; i < myPatch.num_elements(); ++i) {
    std::vector<Sample> samples;
    myPatch.get_samples( i, samples, err ); ASSERT_NO_ERROR(err);
    for (size_t j = 0; j < samples.size(); ++j) {
      double expected = tc.get_weight( myPatch, i, samples[j], err );
      CPPUNIT_ASSERT(!err);
      double read = reader.get_weight( myPatch, i, samples[j], err );
      CPPUNIT_ASSERT(!err);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, read, 1e-12 );
    }
  }
}


