/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */


#include "PySundanceFIATQuadratureAdapter.hpp"
#include <stack>

using namespace std;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

//#ifdef FOOBAR_FIAT_QUADRATURE

namespace Sundance 
{
  FIATQuadratureAdapter::FIATQuadratureAdapter( PyObject *py_quad_factory ,
						int order ) :
    QuadratureFamilyBase( order ) , pts_(3) , wts_(3)
  {
    
    for (i=0;i<3;i++) {
      stack<PyObject *> to_decref;
      // call the factory function with shape and order
      // to get the FIAT quadrature rule for that shape
      
      PyObject *py_quad_rule = 
	PyObject_CallFunction( py_quad_factory , "(ii)" , i+1 , order );
      TEUCHOS_TEST_FOR_EXCEPTION( !py_quad_rule , std::runtime_error , 
			  "Unable to construct quadrature rule" );
      to_decref.push( py_quad_rule );
      
      PyObject *py_quad_points =
	PyObject_CallMethod( py_quad_rule , "get_points" , NULL );
      
      TEUCHOS_TEST_FOR_EXCEPTION( !py_quad_points , std::runtime_error ,
			  "Unable to extract quadrature points" );
      
      to_decref.push( py_quad_points );
      
      PyObject *py_quad_weights =
	PyObject_CallMethod( py_quad_rule , "get_weights" , NULL );
      
      
      TEUCHOS_TEST_FOR_EXCEPTION( !py_quad_weights , std::runtime_error , 
			  "Unable to extract quadrature weights" );
      
      int num_points = PyObject_Length( py_quad_points );
      
      TEUCHOS_TEST_FOR_EXCEPTION( (num_points < 1) , std::runtime_error ,
			  "Not enough points" );
      
      TEUCHOS_TEST_FOR_EXCEPTION( (num_points != PyObject_Length( py_quad_weights ) ) ,
			  std::runtime_error ,
			  "Number of points and weights don't match" );
      
      to_decref.push( py_quad_weights );
      
      
      // get ready to copy data
      pts_[i].resize( num_points );
      wts_[i].resize( num_points );
      
      for (int j=0;j<num_points;j++) {
	pts_[i][j].resize(i+1);				
      }
      
      // copy data for the points
      // lookup functions into list and tuple are borrowed reference
      // and hence don't need to be decrefed
      for (int j=0;j<num_points;j++) {
	// Get a tuple that is the current point
	PyObject *py_pt_cur = PyList_GetItem( py_quad_points , j );
	for (k=0;k<i+1;k++) {
	  PyObject *py_pt_cur_k = PyTuple_GetItem( py_pt_cur , k );
	  pts_[i][j][k] = PyFloat_AsDouble( py_pt_cur_k );
	}
      }
      
      // copy data for weights
      // since they are a Numeric array, I go through the
      // object protocol and have to decref.
      for (int j=0;j<num_points;j++) {
	PyObject *py_j = PyInt_FromLong( (long) j );
	PyObject *py_wt_cur = PyObject_GetItem( py_quad_weights , py_j );
	wts_[i][j] = PyFloat_AsDouble( py_wt_cur );
	to_decref.push( py_j );
	to_decref.push( py_wt_cur );
      }
      
      while( !to_decref.empty() ) {
	PyObject *foo = to_decref.top();
	Py_DECREF( foo );
	to_decref.pop();
      }
    }
  }
  
  
  void FIATQuadratureAdapter::getPoints( const CellType& cellType ,
					 Array<Point>& quadPoints,
					 Array<double>& quadWeights ) const
  {
    int cellDim = dimension( cellType );
    const Array<Point>& pts = pts_[cellDim-1];
    const Array<double>& wts = wts_[cellDim-1];
    quadPoints.resize( pts.size() );
    quadWeights.resize( wts.size() );
    for (i=0;i<pts.size();i++) {
      quadPoints[i].resize( pts[i].dim() );
      for (int j=0;j<pts[i].dim();j++) {
	quadPoints[i][j] = pts[i][j];
      }
      quadWeights[i] = wts[i];
    }
    
    return;
  }
  
  XMLObject FIATQuadratureAdapter::toXML() const 
  {
    XMLObject rtn("FIATQuadrature");
    rtn.addAttribute("order", Teuchos::toString(order()));
    return rtn;
  }
   
}
//#endif
