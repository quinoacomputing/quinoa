/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TAbsQualityMetricTest.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_AWQualityMetric.hpp"
#include "TMPQualityMetricTest.hpp"
#include "Mesquite_TShapeNB1.hpp"
#include "Mesquite_AWMetric.hpp"

using namespace Mesquite;

class FauxAbsShapeMetric : public AWMetric {
  TShapeNB1 mMetric;
public:
  std::string get_name() const { return mMetric.get_name(); }
  bool evaluate( const MsqMatrix<2,2>& A,
                 const MsqMatrix<2,2>& W,
                 double& result,
                 MsqError& err )
    { return mMetric.evaluate( A * inverse(W), result, err ); }
  bool evaluate( const MsqMatrix<3,3>& A,
                 const MsqMatrix<3,3>& W,
                 double& result,
                 MsqError& err )
    { return mMetric.evaluate( A * inverse(W), result, err ); }
};

template <> class TMPTypes<AWQualityMetric> {
public:
  typedef AWMetric MetricType;
  typedef FauxAbsShapeMetric TestType;
};

class AWQualityMetricTest : public TMPQualityMetricTest<AWQualityMetric>
{
  CPPUNIT_TEST_SUITE(AWQualityMetricTest);
  
  REGISTER_TMP_TESTS
  
  CPPUNIT_TEST_SUITE_END();
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AWQualityMetricTest, "AWQualityMetricTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(AWQualityMetricTest, "Unit");
