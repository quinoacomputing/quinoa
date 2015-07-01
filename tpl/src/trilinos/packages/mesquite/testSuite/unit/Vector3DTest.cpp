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
#ifndef __VECTOR3DTEST_H__
#define __VECTOR3DTEST_H__

#include "Vector3D.hpp"
#include "cppunit/extensions/HelperMacros.h"
#include <iostream>

class Vector3DTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Vector3DTest);
  CPPUNIT_TEST (test_default_constructor);
  CPPUNIT_TEST (test_double_constructor);
  CPPUNIT_TEST (test_copy_constructor);
//  CPPUNIT_TEST_EXCEPTION (throw_exception, CppUnit::SignalException);
  CPPUNIT_TEST_SUITE_END();

public:
  Vector3DTest()
    {}
  
  void test_default_constructor()
    {
      Mesquite::Vector3D v;
      CPPUNIT_ASSERT(v.x() == 0);
      CPPUNIT_ASSERT(v.y() == 0);
      CPPUNIT_ASSERT(v.z() == 0);
    }
  void test_double_constructor()
    {
      Mesquite::Vector3D v(3, 2, 1);
      CPPUNIT_ASSERT(v.x() == 3);
      CPPUNIT_ASSERT(v.y() == 2);
      CPPUNIT_ASSERT(v.z() == 1);
    }
  void test_copy_constructor()
    {
      Mesquite::Vector3D v2(3, 2, 1);
      Mesquite::Vector3D v(v2);
      
      CPPUNIT_ASSERT(v.x() == 3);
      CPPUNIT_ASSERT(v.y() == 2);
      CPPUNIT_ASSERT(v.z() == 1);
    }
  void throw_exception()
    {
      Mesquite::Vector3D* v = NULL;
      double d = v->x();
      std::cout << d << std::endl;
    }
};

class Vector3DTest2 : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Vector3DTest);
  CPPUNIT_TEST_SUITE_END();

public:
  Vector3DTest2()
    {}
  
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Vector3DTest, "Unit");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Vector3DTest, "Misc");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(Vector3DTest, "Vector3DTest");

#endif
