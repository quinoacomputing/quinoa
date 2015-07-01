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
// ORIG-DATE: 13-Jan-03 at 09:05:56
//  LAST-MOD:  9-Apr-03 at 17:17:26 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqHessianTest.cpp

Unit testing of the MsqHessian class. 

*/
// DESCRIP-END.
//


#include "PatchDataInstances.hpp"
#include "MsqHessian.hpp"

#include <math.h>

#include "cppunit/extensions/HelperMacros.h"

using namespace Mesquite;
using std::cout;
using std::cerr;
using std::endl;

class MsqHessianTest
  : public CppUnit::TestFixture,
    public Mesquite::MsqHessian
{

private:
  CPPUNIT_TEST_SUITE(MsqHessianTest);
  CPPUNIT_TEST (test_initialize);
  CPPUNIT_TEST (test_axpy);
  CPPUNIT_TEST (test_cg_solver);
  CPPUNIT_TEST (test_cholesky_preconditioner);
  CPPUNIT_TEST_SUITE_END();

private:
   
  PatchData twoTriangles; 

public:
  void setUp()
  {
    MsqPrintError err(cout);
    create_two_tri_patch(twoTriangles, err); CPPUNIT_ASSERT(!err);
  }

  void tearDown()
  {
    destroy_patch_with_domain( twoTriangles );
  }

public:
  MsqHessianTest()
  {}

  void accumulate_entries( const PatchData &pd, 
                           size_t elem_index,
                           const Matrix3D* mat3d_array, 
                           MsqError &err );

  void test_initialize()
  {
    MsqPrintError err(cout);
    
    MsqHessian::initialize(twoTriangles, err); CPPUNIT_ASSERT(!err);

    // Checks values of mRowStart are correct.
    size_t i, row_start[] = {0, 4, 7, 8, 9};
    for (i=0; i<5; ++i) 
      CPPUNIT_ASSERT( mRowStart[i] == row_start[i] );

    // Checks values of mColIndex are correct.
    size_t col_index[] = {0,1,2,3,1,2,3,2,3};
    for (i=0; i<9; ++i) 
      CPPUNIT_ASSERT( mColIndex[i] == col_index[i] );
  }

  void test_axpy()
  {
    size_t i;
	MsqPrintError err(cout);
    
    MsqHessian::initialize(twoTriangles, err); CPPUNIT_ASSERT(!err);

    size_t hs = MsqHessian::size();

    Vector3D* res = new Vector3D[hs];
    Vector3D* x = new Vector3D[hs];
    Vector3D* y = new Vector3D[hs];
    Vector3D* ans = new Vector3D[hs];

    Matrix3D blocks[6]; // 6 blocks correspond to a triangular element (n+1)n/2 .

    blocks[0] = "4 4 7   4 5 7   7 7 3 ";
    blocks[1] = "4 8 7   3 5 7   1 2 3 ";
    blocks[2] = "4 4 7   6 5 9   1 8 5 ";
    blocks[3] = "4 4 2   4 5 3   2 3 3 ";
    blocks[4] = "2 4 7   3 2 7   1 4 3 ";
    blocks[5] = "8 4 9   4 5 7   9 7 3 ";

    accumulate_entries(twoTriangles, 0, blocks, err); CPPUNIT_ASSERT(!err);

    blocks[2] += blocks[5];
    blocks[5] = blocks[3];
    
    accumulate_entries(twoTriangles, 1, blocks, err); CPPUNIT_ASSERT(!err);

    Matrix3D entries_6_ans("2 3 1   4 2 4   7 7 3");
    CPPUNIT_ASSERT( mEntries[6] == entries_6_ans );

    x[0].set(4, 5, 6);
    x[1].set(2, 5, 9);
    x[2].set(1, 2, 6);
    x[3].set(1, 5, 9);

    y[0].set(0, 0, 0);
    y[1].set(0, 0, 0);
    y[2].set(0, 0, 0);
    y[3].set(0, 0, 0);

    axpy(res, hs, *this, x, hs, y, hs, err); CPPUNIT_ASSERT(!err);

    ans[0].set(636, 635, 453);
    ans[1].set(365, 460, 461);
    ans[2].set(150, 199, 220);
    ans[3].set(166, 204, 174);

    for (i=0; i<hs; ++i)
      CPPUNIT_ASSERT( res[i] == ans[i] );

    
    y[0].set(3, 2, 6);
    y[1].set(1, 2, 4);
    y[2].set(3, 6, 9);
    y[3].set(2, 4, 4);

    ans[0].set(639, 637, 459);
    ans[1].set(366, 462, 465);
    ans[2].set(153, 205, 229);
    ans[3].set(168, 208, 178);

    axpy(res, hs, *this, x, hs, y, hs, err); CPPUNIT_ASSERT(!err);

    for (i=0; i<hs; ++i)
      CPPUNIT_ASSERT( res[i] == ans[i] );

    delete[] res;
    delete[] x;
    delete[] y;
    delete[] ans;
  }

  
  void test_cg_solver()
  {
    size_t i;
    MsqPrintError err(cout);
    
    MsqHessian::initialize(twoTriangles, err); CPPUNIT_ASSERT(!err);

    size_t hs = MsqHessian::size();

    Vector3D* res = new Vector3D[hs];
    Vector3D* x = new Vector3D[hs];
    Vector3D* y = new Vector3D[hs];
    Vector3D* ans = new Vector3D[hs];

    Matrix3D blocks[6]; // 6 blocks correspond to a triangular element (n+!)n/2 .

    blocks[0] = "14 4 7   4 15 7   7 7 13 ";
    blocks[1] = "4 8 7   3 5 7   1 2 3 ";
    blocks[2] = "4 4 7   6 5 9   1 8 5 ";
    blocks[3] = "14 4 2   4 15 3   2 3 13 ";
    blocks[4] = "2 4 7   3 2 7   1 4 3 ";
    blocks[5] = "18 4 9   4 15 7   9 7 13 ";

    accumulate_entries(twoTriangles, 0, blocks, err); CPPUNIT_ASSERT(!err);

    blocks[2] -= blocks[5];
    blocks[5] = blocks[3];
    
    accumulate_entries(twoTriangles, 1, blocks, err); CPPUNIT_ASSERT(!err);
        
    y[0].set(3, 2, 6);
    y[1].set(1, 2, 4);
    y[2].set(3, 6, 9);
    y[3].set(2, 4, 4);

    cg_solver(x, y, err); CPPUNIT_ASSERT(!err);

//     for (int i=0; i<4; ++i) 
//       cout << x[i];
    
    ans[0].set(3.2338, 7.6431, -7.0735);
    ans[1].set(-1.0068, 4.5520, -4.6628);
    ans[2].set(-0.4361, 4.7640, -8.1006);
    ans[3].set(-0.1218, -0.8817, -4.5571);

    for (i=0; i<hs; ++i)
      for (short j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(x[i][j], ans[i][j], 10e-1);

    delete[] res;
    delete[] x;
    delete[] y;
    delete[] ans;    
          
  }

    void test_cholesky_preconditioner()
  {
    MsqPrintError err(cout);

    MsqHessian::initialize(twoTriangles, err); CPPUNIT_ASSERT(!err);
    MsqHessian::zero_out();

    Matrix3D blocks[6]; // 6 blocks correspond to a triangular element (n+1)n/2 .

    blocks[0] = " 2  1    1 "
                " 1  2.5  0.5 "
                " 1  0.5  2.5 ";
    blocks[1] = "0 0 0   0 0 0   0 0 0 ";
    blocks[2] = "0 0 0   0 0 0   0 0 0 ";
    blocks[3] = "0 0 0   0 0 0   0 0 0 ";
    blocks[4] = "0 0 0   0 0 0   0 0 0 ";
    blocks[5] = "0 0 0   0 0 0   0 0 0 ";

    accumulate_entries(twoTriangles, 0, blocks, err); CPPUNIT_ASSERT(!err);

    MsqHessian::compute_preconditioner(err); CPPUNIT_ASSERT(!err);
    Matrix3D block_0 = mPreconditioner[0];

    Matrix3D correct(" 0.5  0.5  0.5 "
                     " 0    0.5  0  "
                     " 0    0    0.5 ");

    for (short i=0; i<3; ++i)
      for (short j=0; j<3; ++j)
        CPPUNIT_ASSERT_DOUBLES_EQUAL( block_0[i][j], correct[i][j], 10e-10);

//     cout << "block 0: \n" << block_0 << endl;
//     cout << "correct: \n" << correct << endl;
    
  }
  
};

void MsqHessianTest::accumulate_entries( const PatchData &pd, 
                                         size_t elem_index,
                                         const Matrix3D* mat3d_array, 
                                         MsqError &err )
{
  const MsqMeshEntity& elem = pd.element_by_index( elem_index );
  const size_t nv = elem.vertex_count();
  const size_t* v = elem.get_vertex_index_array();
  for (size_t r = 0; r < nv; ++r) {
    for (size_t c = r; c < nv; ++c) {
      add( v[r], v[c], *mat3d_array, err );
      CPPUNIT_ASSERT(!err);
      ++mat3d_array;
    }
  }
}


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqHessianTest, "MsqHessianTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(MsqHessianTest, "Unit");
