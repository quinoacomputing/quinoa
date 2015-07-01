/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 12-Nov-02 at 18:05:56
//  LAST-MOD:  5-May-03 at 15:59:29 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqFreeVertexIndexIteratorTest.cpp

Unit testing of various functions in the MsqFreeVertexIndexIterator class. 

 */
// DESCRIP-END.
//


#include "MsqFreeVertexIndexIterator.hpp"
#include "PatchDataInstances.hpp"

#include <math.h>
#include <iostream>

#include "cppunit/extensions/HelperMacros.h"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class MsqFreeVertexIndexIteratorTest : public CppUnit::TestFixture
{

private:
  CPPUNIT_TEST_SUITE(MsqFreeVertexIndexIteratorTest);
  CPPUNIT_TEST (test_hard_fixed_flags);
  CPPUNIT_TEST (test_soft_fixed_flags);
  CPPUNIT_TEST_SUITE_END();

private:
  PatchData pd;

public:
  void setUp()
  {
    MsqPrintError err(cout);

    /*      7____6____5___11
            |    |    |    |
            | 2  |  3 | 5  |
            8-_  |  _-4---10       vertex 1 is at (0,0)
            |  -_0_-  |    |       vertex 11 is at (3,2)
            | 0  |  1 | 4  |
            1----2----3----9
    */
    create_six_quads_patch_with_domain(pd, err); 
  }

  void tearDown()
  {
    destroy_patch_with_domain(pd);
  }
  
public:
  MsqFreeVertexIndexIteratorTest()
    {}
  
  void test_hard_fixed_flags()
  {   
     MsqPrintError err(cout);
     int indices[10];
     int i=0;
     MsqFreeVertexIndexIterator ind(pd, err);
     ind.reset();
     while (ind.next()) {
        indices[i] = ind.value();
//         cout << i << "th free vertex value: " << ind.value() << endl; 
        ++i;
     } 

     CPPUNIT_ASSERT(i==2); // number of free vertices.
     CPPUNIT_ASSERT(pd.vertex_by_index(indices[0]).is_free_vertex());
     CPPUNIT_ASSERT(pd.vertex_by_index(indices[1]).is_free_vertex());
  }

  void test_soft_fixed_flags()
  {   
     MsqPrintError err(cout);
     pd.set_vertex_culled( 0 );

     int indices[10];
     int i=0;
     MsqFreeVertexIndexIterator ind(pd, err);
     ind.reset();
     while (ind.next()) {
        indices[i] = ind.value();
        ++i;
     } 

     CPPUNIT_ASSERT(i==1); // number of free vertices.
     CPPUNIT_ASSERT(pd.vertex_by_index(indices[0]).is_free_vertex());
  }


};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqFreeVertexIndexIteratorTest, "MsqFreeVertexIndexIteratorTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqFreeVertexIndexIteratorTest, "Unit");
