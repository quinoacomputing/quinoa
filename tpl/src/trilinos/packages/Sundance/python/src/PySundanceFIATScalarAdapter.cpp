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


#include "PySundanceFIATScalarAdapter.hpp"
#include <stack>
#include <iostream>
using namespace std;
using namespace Sundance;
using namespace Sundance;


static int line_sdvert_to_fvert[] = {0,1};
static int line_sdline_to_fline[] = {0};
static int *line_sd_to_fiat[] = {line_sdvert_to_fvert,line_sdline_to_fline};
static int tri_sdvert_to_fvert[] = {0,1,2};
static int tri_sdline_to_fline[] = {2,0,1};
static int tri_sdtri_to_ftri[] = {0};
static int *tri_sd_to_fiat[] = {tri_sdvert_to_fvert,
                                tri_sdline_to_fline,
                                tri_sdtri_to_ftri};
static int tet_sdvert_to_fvert[] = {0,1,2,3};
static int tet_sdline_to_fline[] = {2,0,1,3,4,5};
static int tet_sdtri_to_ftri[] = {0,1,2,3};
static int tet_sdtet_to_ftet[] = {0};

static int *tet_sd_to_fiat[] = {tet_sdvert_to_fvert,
                                tet_sdline_to_fline,
                                tet_sdtri_to_ftri,
                                tet_sdtet_to_ftet};

static int **sd_to_fiat[] = {NULL,
                             line_sd_to_fiat,
                             tri_sd_to_fiat,
                             tet_sd_to_fiat};

static CellType sdim_to_cellType[] = 
{PointCell,LineCell,TriangleCell,TetCell};

#define HAVE_PY_FIAT
#ifdef HAVE_PY_FIAT

namespace Sundance
{
FIATScalarAdapter::FIATScalarAdapter( PyObject *pyfamilyclass ,
  int order ) :
  order_( order )
{
  stack<PyObject *> to_decref;

  /* instantiate basis for each shape */
  bases_.resize( 3 );
  for (int i=0;i<3;i++) {
    PyObject *arglist = Py_BuildValue( "(ii)" , i+1 , order_ );
    if (!arglist) cout << "barf" << std::endl;
    bases_[i] = PyObject_CallObject( pyfamilyclass , arglist  );
    if (!bases_[i]) cout << "puke" << std::endl;
    to_decref.push( arglist );
  }

  /* pretabulate all the dof in the constructor */
  /* one dof array for each cellType I instantiate */
  dof_.resize( 3 );
  for (int cd=1;cd<=3;cd++) {
    dof_[cd-1].resize( cd + 1 );
    PyObject *basis = bases_[cd-1];
    PyObject *dual_basis = PyObject_CallMethod( basis , "dual_basis" ,
      NULL );
    to_decref.push( dual_basis );
    for (int i=0;i<=cd;i++) {
      PyObject *nodes_per_dim =
        PyObject_CallMethod( dual_basis , "getNodeIDs" ,
          "i" , i );
      to_decref.push( nodes_per_dim );
      int num_facets = PyObject_Length( nodes_per_dim );
      dof_[cd-1][i].resize( num_facets );
      for (int j=0;j<num_facets;j++) {
        //cout << "facet number " << j << std::endl;
        PyObject *pyj = PyInt_FromLong( (long) j );
        to_decref.push( pyj );
        PyObject *nodes_per_facet =
          PyObject_GetItem( nodes_per_dim , pyj );
        to_decref.push( nodes_per_facet );
        int num_nodes_this_facet = PyObject_Length( nodes_per_facet );
	  
        dof_[cd-1][i][j].resize( num_nodes_this_facet );
        for (int k=0;k<num_nodes_this_facet;k++) {
          PyObject *pyk = PyInt_FromLong( (long) k );
          to_decref.push( pyk );
          PyObject *pynodecur = PyObject_GetItem( nodes_per_facet , pyk );
          to_decref.push( pynodecur );
          dof_[cd-1][i][j][k] = (int) PyInt_AsLong( pynodecur );
        }
      }
    }
  }

  while (!to_decref.empty()) {
    PyObject *foo = to_decref.top();
    Py_DECREF( foo );
    to_decref.pop();
  }	

}
  
FIATScalarAdapter::~FIATScalarAdapter()
{
  for (int i=0;i<3;i++) {
    Py_DECREF( bases_[i] );
  }
}




bool FIATScalarAdapter::supportsCellTypePair(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case LineCell:
      switch(cellType)
      {
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    case TriangleCell:
      switch(cellType)
      {
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    case TetCell:
      switch(cellType)
      {
        case TetCell:
        case TriangleCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
    default:
      return false;
  }
}

/* loops over pretabulated data structures and copies into dofs */
void FIATScalarAdapter::getReferenceDOFs(const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const
{
//		cout << "getting local dof" << std::endl;
  if (PointCell == cellType) {
    dofs.resize(1);
    dofs[0] = tuple<Array<int> >(tuple(0));
  }
  else {
    int cd = dimension( cellType );
    const Array<Array<Array<int> > >& dofs_cur = dof_[cd-1];
    dofs.resize( dofs_cur.size() );
    for (i=0;i<dofs_cur.size();i++) {
      dofs[i].resize( dofs_cur[i].size() );
      for (j=0;j<dofs_cur[i].size();j++) {
        dofs[i][j].resize( dofs_cur[i][j].size() );
        for (k=0;k<dofs_cur[i][j].size();k++) {
          dofs[i][j][k] = dofs_cur[i][j][k];
        }
      }
    }
  }
//  		cout << "done getting local dof" << std::endl;
  return;
}

// sum over spatial dimensions up to and including spatialDim
// for the given cellType
int FIATScalarAdapter::nReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const 
{
//		cout << "getting nnodes" << std::endl;
  if (PointCell == cellType) {
//			cout << "done getting nodes" << std::endl;
    return 1;
  }
  else {
    int cellDim = dimension( cellType );
    int nn = 0;
    const Array<Array<Array<int> > >& dofs_cur = dof_[cellDim-1];
    for (int i=0;i<=cellDim;i++) {
      nn += numFacets( cellType , i ) * dofs_cur[i][0].size();
    }

//			cout << "done getting nnodes" << std::endl;

    return nn;
  }
}


/* eval the basis associated with spatialDim,
   extract the values of basis functions associated only
   with facet 0 of cellType and all lower-dimensional
   facets covering it */
void FIATScalarAdapter::refEval(
  const CellType& maximalCellType,
  const CellType& cellType,
  const Array<Point>& pts,
  const MultiIndex& deriv,
  Array<Array<Array<double> > >& result) const
{
  int spatialDim = dimension(maximalCellType);

  result.resize(1);
  if (PointCell == cellType) {
    result[0] = tuple<Array<double> >(tuple(1.0));
  }
  else {
// 			cout << "refevaling" << std::endl;
// 			cout << "spatial dim: " << spatialDim << std::endl;
// 			cout << "cell dim: " << dimension( cellType ) << std::endl;
// 			cout << "pts sizes: " << pts.size() << " " << pts[0].dim() << std::endl;
// 			cout << nNodes( dimension(cellType) , cellType ) << " nodes" << std::endl;
    stack<PyObject *> to_decref;
    int cellDim = dimension( cellType );
 
    result[0].resize( pts.size() );
    int nn = nReferenceDOFs( maximalCellType, cellType );
    for (i=0;i<pts.size();i++) {
      result[0][i].resize(nn);
    }
 
// 			cout << "result sizes" << std::endl;
//  			for (i=0;i<result.size();i++) {
//  				cout << result[0][i].size() << std::endl;
//  			}
// 			cout << "making points" << std::endl;
    // This is the list of points converted into a Python list
    PyObject *py_list_of_points = PyList_New(pts.size());
    TEUCHOS_TEST_FOR_EXCEPTION( !py_list_of_points , std::runtime_error, 
      "Unable to create list" );
    to_decref.push( py_list_of_points );
    for (int i=0;i<pts.size();i++) {
      // Create a Python tuple for the point, converting from
      // Sundance (0,1)-based coordinates to FIAT (-1,1)-based 
      // coordinates
      PyObject *py_pt_cur = PyTuple_New( spatialDim );
      TEUCHOS_TEST_FOR_EXCEPTION( !py_pt_cur , std::runtime_error ,
        "Unable to create tuple" );
      for (int j=0;j<cellDim;j++) {
        double coord_cur = 2.0 * ( pts[i][j] - 0.5 );
        PyObject *py_coord = PyFloat_FromDouble( coord_cur );
        TEUCHOS_TEST_FOR_EXCEPTION( !py_coord , std::runtime_error ,
          "Unable to create PyFloat" );
        int msg = PyTuple_SetItem( py_pt_cur , j , py_coord );
        TEUCHOS_TEST_FOR_EXCEPTION( msg==-1 , std::runtime_error ,
          "Unable to set tuple item" );
        // PyTuple_SetItem steals a reference to py_coord;
        // not added to the decref stack
      }
      for (int j=cellDim;j<spatialDim;j++) {
        PyObject *py_coord = PyFloat_FromDouble( -1.0 );
        int msg = PyTuple_SetItem( py_pt_cur , j , py_coord );
        TEUCHOS_TEST_FOR_EXCEPTION( msg==-1 , std::runtime_error ,
          "Unable to set tuple item" );
        // PyTuple_SetItem steals a reference to py_coord;
        // not added to the decref stack
      }
//				cout << "putting into list" << std::endl;
      int msg = PyList_SetItem( py_list_of_points , i , py_pt_cur );
      // reference to py_pt_cur stolen
      TEUCHOS_TEST_FOR_EXCEPTION( msg==-1 , std::runtime_error ,
        "Unable to set tuple item" );
//				cout << "done putting into list" << std::endl;
    }

//			cout << "done making points" << std::endl;

    // Extract the function space from the basis
    PyObject *py_basis = bases_[spatialDim-1];
    PyObject *py_function_space = 
      PyObject_CallMethod( py_basis , "function_space" , NULL );
    TEUCHOS_TEST_FOR_EXCEPTION( !py_function_space , std::runtime_error, 
      "Could not extract function space" );
    to_decref.push( py_function_space );

			

    // convert the multiindex from Sundance into a Python tuple
    PyObject *py_alpha = PyTuple_New( spatialDim );
    TEUCHOS_TEST_FOR_EXCEPTION( !py_alpha , std::runtime_error, 
      "Unable to create new tuple" );
    to_decref.push( py_alpha );
    for (int i=0;i<spatialDim;i++) {
      PyObject *py_alpha_i = PyInt_FromLong( (long) deriv[i] );
      TEUCHOS_TEST_FOR_EXCEPTION( !py_alpha_i , std::runtime_error,
        "Unable to create PyInt" );
      int msg = PyTuple_SetItem( py_alpha , i , py_alpha_i );
      // reference to py_alpha_i stolen
      TEUCHOS_TEST_FOR_EXCEPTION( msg==-1 , std::runtime_error ,
        "Unable to set tuple item" );
    }

    TEUCHOS_TEST_FOR_EXCEPTION( !PyObject_HasAttrString( py_function_space ,
        "multi_deriv_all" ) ,
      std::runtime_error , "???" );

    // Get the set of all partial derivatives
    PyObject *py_deriv_space =
      PyObject_CallMethod( py_function_space , "multi_deriv_all" , 
        "(O)" , py_alpha );
    TEUCHOS_TEST_FOR_EXCEPTION( !py_deriv_space , std::runtime_error ,
      "Unable to take derivatives" );
    to_decref.push( py_deriv_space );

    // This should be a Numeric.array object.  I can
    // go to the low-level Numeric API and grab the C pointer to
    // the data directly if we have a speed problem.
    // 
    // This function should give a 2d array that is num_pts by
    // the number of basis functions associated with cellType
    // (and things that cover it)
    PyObject *py_tabulation =
      PyObject_CallMethod( py_deriv_space , "tabulate" ,
        "(O)" , py_list_of_points );
    TEUCHOS_TEST_FOR_EXCEPTION( !py_tabulation , std::runtime_error ,
      "Unable to tabulate derivatives" );
    to_decref.push( py_tabulation );

    Array<Array<Array<int> > > dofs;
    getReferenceDOFs(maximalCellType, cellType, dofs );


    double factor = pow( 2.0 , deriv.order() );

    int cur = 0;
    int **sd_to_fiat_spd = sd_to_fiat[spatialDim];
    for (d=0;d<=(unsigned)cellDim;d++) {
//				cout << "copying points for dimension " << d << std::endl;
      int *sd_to_fiat_spd_d = sd_to_fiat_spd[d];
      for (int e=0;e<numFacets(cellType,d);e++) {
//					cout << "\tcopying points for facet " << e << std::endl;
        int fiat_e = sd_to_fiat_spd_d[e];
//					cout << "fiat_e: " << fiat_e << std::endl;
        for (int n=0;n<dofs[d][e].size();n++) {
//						cout << "\t\tcopying points for local bf " << n << std::endl;
//						cout << "\t\tcur " << cur << std::endl;
          for (unsigned p=0;p<pts.size();p++) {
//							cout << "\t\t\tcopying value for point " << p << std::endl;
//							cout << "ij tuple" << std::endl;
            PyObject *py_ij_tuple = 
              Py_BuildValue( "(ii)" , dofs[d][fiat_e][n] , p );
            to_decref.push( py_ij_tuple );
//							cout << "lookup" << std::endl;
            PyObject *py_tab_cur_ij = 
              PyObject_GetItem( py_tabulation , py_ij_tuple );
            to_decref.push( py_tab_cur_ij );
            result[0][p][cur] = factor * PyFloat_AsDouble( py_tab_cur_ij );
          }
          cur++;
        }
      }
    }

//			cout << "done refevaling" << std::endl;


    while (!to_decref.empty()) {
      PyObject *cur = to_decref.top();
      Py_DECREF( cur );
      to_decref.pop();
    }
  }
		
  return;
}
}
#endif
