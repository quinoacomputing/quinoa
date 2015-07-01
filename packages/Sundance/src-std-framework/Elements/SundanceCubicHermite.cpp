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

#include "SundanceCubicHermite.hpp"
#include "SundanceADReal.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;



bool CubicHermite::supportsCellTypePair(
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
//     case TetCell:
//       switch(cellType)
//       {
//         case TetCell:
//         case TriangleCell:
//         case LineCell:
//         case PointCell:
//           return true;
//         default:
//           return false;
//       }
    default:
      return false;
  }
}

void CubicHermite::print(std::ostream& os) const 
{
  os << "CubicHermite";
}

int CubicHermite::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  switch(maximalCellType)
  {
    case LineCell:
      switch(cellType)
      {
        case PointCell:
          return 2;
        case LineCell:
          return 0;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
          return -1;
      }
      break;
    case TriangleCell:
      switch(cellType)
      {
        case PointCell:
          return 3;
        case LineCell:
          return 0;
        case TriangleCell:
          return 1;
        default:
          TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
          return -1;
      }
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
      return -1;
  }
  
}

void CubicHermite::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;
  switch(cellType)
  {
    case PointCell:
    {
      dofs.resize(1);
      if (maximalCellType==LineCell)
        dofs[0] = tuple<Aint>(tuple(0,1));
      else if (maximalCellType==TriangleCell)
        dofs[0] = tuple<Aint>(tuple(0,1,2));
      else TEUCHOS_TEST_FOR_EXCEPT(1);
      return;
    }
    break;
    case LineCell:
    {
      dofs.resize(2);
      dofs[0].resize(2);
      if (maximalCellType==LineCell)
      {
        dofs[0][0].resize(2);
        dofs[0][0][0] = 0;
        dofs[0][0][1] = 1;
        dofs[0][1].resize(2);
        dofs[0][1][0] = 2;
        dofs[0][1][1] = 3;
      }
      else if (maximalCellType==TriangleCell)
      {
        dofs[0][0].resize(3);
        dofs[0][0][0] = 0;
        dofs[0][0][1] = 1;
        dofs[0][0][2] = 2;
        dofs[0][1].resize(3);
        dofs[0][1][0] = 3;
        dofs[0][1][1] = 4;
        dofs[0][1][1] = 5;
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPT(1);
      }
      dofs[1].resize(1);
      dofs[1][0].resize(0);
      return;
    }
    break;
    case TriangleCell:
    {
      dofs.resize(3);
      dofs[0].resize(3);
      dofs[0][0] = tuple(0,1,2);
      dofs[0][1] = tuple(3,4,5);
      dofs[0][2] = tuple(6,7,8);
      dofs[1].resize(3);
      dofs[1][0].resize(0);
      dofs[1][1].resize(0);
      dofs[1][2].resize(0);
      dofs[2].resize(1);
      dofs[2][0].resize(1);
      dofs[2][0][0] = 9;
    }
    break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
        << cellType << " not implemented in CubicHermite basis");
  }
}


void CubicHermite::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    std::runtime_error,
    "cannot evaluate spatial derivative " << sds << " on CubicHermite basis");
  const MultiIndex& deriv = sds.mi();
  typedef Array<double> Adouble;
  result.resize(1);
  result[0].resize(pts.length());

  switch(cellType)
  {
    case PointCell:
      result[0] = tuple<Adouble>(tuple(1.0));
      return;
    case LineCell:
      for (int i=0; i<pts.length(); i++)
      {
        evalOnLine(pts[i], deriv, result[0][i]);
      }
      return;
    case TriangleCell:
      for (int i=0; i<pts.length(); i++)
      {
        evalOnTriangle(pts[i], deriv, result[0][i]);
      }
      return;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "CubicHermite::refEval() unimplemented for cell type "
        << cellType);

  }
}

/* ---------- evaluation on different cell types -------------- */

void CubicHermite::evalOnLine(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  result.resize(4);
  ADReal x = ADReal(pt[0],0,1);
  Array<ADReal> tmp(4);

  tmp[0] = 1 + x * x * ( -3 + 2 * x );
  tmp[1] = x * ( 1 + (x - 2 ) * x );
  tmp[2] = ( 3 - 2*x ) * x * x;
  tmp[3] = (-1+x)*x*x;

  for (int i=0; i<tmp.length(); i++)
  {
    if (deriv.order()==0) 
    {
      result[i] = tmp[i].value();
    }
    else 
    {
      result[i] = tmp[i].gradient()[0];
    }
  }  
  return;
}

void CubicHermite::evalOnTriangle(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const

{
  result.resize(10);
  ADReal x = ADReal(pt[0], 0, 2);
  ADReal y = ADReal(pt[1], 1, 2);
  ADReal one(1.0, 2);
  
  Array<ADReal> tmp(10);

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
    << y.value());
 
  tmp[0] = 1 - 3*x*x + 2 * x*x*x - 13*x*y + 13*x*x*y - 3*y*y + 13 *x*y*y + 2 *y*y*y;
  tmp[1] = x - 2 *x*x + x*x*x - 3*x*y + 3*x*x*y + 2*x*y*y;
  tmp[2] = y - 3 *x *y + 2* x*x* y - 2* y*y + 3* x*y*y + y*y*y;
  tmp[3] = 3 * x*x - 2* x*x*x - 7* x* y + 7* x*x *y + 7* x*y*y;
  tmp[4] = -x*x + x*x*x + 2*x *y - 2* x*x* y - 2* x* y*y;
  tmp[5] = -x* y + 2* x*x* y + x* y*y;
  tmp[6] = -7* x* y + 7* x*x*y + 3* y*y + 7* x* y*y - 2* y*y*y;
  tmp[7] = -x *y + x*x* y + 2* x* y*y;
  tmp[8] = 2 *x *y - 2* x*x* y - y*y - 2* x* y*y + y*y*y;
  tmp[9] = 27* x *y - 27* x*x* y - 27* x* y*y;

  for (int i=0; i<tmp.length(); i++)
  {
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
}

void CubicHermite::evalOnTet(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  ADReal x = ADReal(pt[0], 0, 3);
  ADReal y = ADReal(pt[1], 1, 3);
  ADReal z = ADReal(pt[2], 2, 3);
  ADReal one(1.0, 3);
  
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

void CubicHermite::preApplyTransformation( const CellType &maxCellType ,
					   const Mesh &mesh, 
					   const Array<int> &cellLIDs,
					   const CellJacobianBatch& JVol,
					   RCP<Array<double> >& A ) const
  {
    switch(maxCellType)
      {
      case TriangleCell:
	preApplyTransformationTriangle( mesh , cellLIDs, JVol , A );
	break;
      default:
	TEUCHOS_TEST_FOR_EXCEPT(1);
	break;
      }
    return;
  }

void CubicHermite::postApplyTransformation( const CellType &maxCellType ,
					    const Mesh &mesh, 
					    const Array<int> &cellLIDs,
					    const CellJacobianBatch& JVol,
					    RCP<Array<double> >& A ) const
  {
    switch(maxCellType)
      {
      case TriangleCell:
	postApplyTransformationTriangle( mesh , cellLIDs, JVol , A );
	break;
      default:
	TEUCHOS_TEST_FOR_EXCEPT(1);
	break;
      }
    return;
  }

void CubicHermite::preApplyTransformationTranspose( const CellType &maxCellType ,
						    const Mesh &mesh, 
						    const Array<int> &cellLIDs,
						    const CellJacobianBatch& JVol ,
						    Array<double> & A ) const
{
  switch(maxCellType)
    {
    case TriangleCell:
      preApplyTransformationTransposeTriangle( mesh , cellLIDs, JVol , A );
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(1);
      break;
    }
  return;
}


void CubicHermite::preApplyTransformationTriangle( const Mesh &mesh, 
						   const Array<int> &cellLIDs,
						   const CellJacobianBatch& JVol,
						   RCP<Array<double> >& A) const
{
  Array<double> &Anoptr = *A;
  
  Array<double> cellVertH;
  getVertexH( mesh , cellLIDs , cellVertH );
  
  
  // this applies M from the left on each cell
  // A for each cell has 10 rows because it is Hermite
  // so, this gives the number of columns per cell
  const int numCols = Anoptr.size() / JVol.numCells() / 10;
  
  for (int i=0;i<JVol.numCells();i++) 
    {
      const int cell_start = i * numCols * 10;
      const double *invJ = JVol.jVals(i);
      
      for (int j=0;j<numCols;j++) 
	{
	  const int col_start = cell_start + 10 * j;
	  const double a1 = Anoptr[col_start + 1];
	  const double a2 = Anoptr[col_start + 2];
	  const double a4 = Anoptr[col_start + 4];
	  const double a5 = Anoptr[col_start + 5];
	  const double a7 = Anoptr[col_start + 7];
	  const double a8 = Anoptr[col_start + 8];
	  Anoptr[col_start+1] = (invJ[0]*a1 + invJ[1]*a2)/cellVertH[3*i];
	  Anoptr[col_start+2] = (invJ[2]*a1 + invJ[3]*a2)/cellVertH[3*i];
	  Anoptr[col_start+4] = (invJ[0]*a4 + invJ[1]*a5)/cellVertH[3*i+1];
	  Anoptr[col_start+5] = (invJ[2]*a4 + invJ[3]*a5)/cellVertH[3*i+1];
	  Anoptr[col_start+7] = (invJ[0]*a7 + invJ[1]*a8)/cellVertH[3*i+2];
	  Anoptr[col_start+8] = (invJ[2]*a7 + invJ[3]*a8)/cellVertH[3*i+2];
	}
    }
  
}

void CubicHermite::postApplyTransformationTriangle( const Mesh &mesh,
						    const Array<int> &cellLIDs , 
						    const CellJacobianBatch& JVol,
						    RCP<Array<double> >& A ) const
{
  Array<double> &Anoptr = *A;
  
  Array<double> cellVertH;
  getVertexH( mesh , cellLIDs , cellVertH );
  
  
  const int numRows = Anoptr.size() / 10 / JVol.numCells();
  
  for (int i=0;i<JVol.numCells();i++) 
    {
      const double *invJ = JVol.jVals(i);
      
      const int cell_start = i * numRows * 10;
      // handle columns 1 and 2
      for (int j=0;j<numRows;j++) 
	{
	  const double a = Anoptr[cell_start + numRows + j];
	  const double b = Anoptr[cell_start + 2*numRows + j];
	  Anoptr[cell_start + numRows + j] = (invJ[0] * a + invJ[1] * b)/cellVertH[3*i];
	  Anoptr[cell_start + 2*numRows + j] = (invJ[2] * a + invJ[3] * b)/cellVertH[3*i];
	}
      
      // handle columns 4 and 5
      for (int j=0;j<numRows;j++) 
	{
	  const double a = Anoptr[cell_start + 4*numRows + j];
	  const double b = Anoptr[cell_start + 5*numRows + j];
	  Anoptr[cell_start + 4*numRows + j] = (invJ[0] * a + invJ[1] * b)/cellVertH[3*i+1];
	  Anoptr[cell_start + 5*numRows + j] = (invJ[2] * a + invJ[3] * b)/cellVertH[3*i+1];
	}
      
      // handle columns 7 and 8
      for (int j=0;j<numRows;j++) 
	{
	  const double a = Anoptr[cell_start + 7*numRows + j];
	  const double b = Anoptr[cell_start + 8*numRows + j];
	  Anoptr[cell_start + 7*numRows + j] = (invJ[0] * a + invJ[1] * b)/cellVertH[3*i+2];
	  Anoptr[cell_start + 8*numRows + j] = (invJ[2] * a + invJ[3] * b)/cellVertH[3*i+2];
	}
      
    }
  
  return;
}

void CubicHermite::preApplyTransformationTransposeTriangle( const Mesh &mesh,
							    const Array<int> &cellLIDs,
							    const CellJacobianBatch& JVol,
							    Array<double>& A ) const
{
  // this applies M from the left on each cell
  // A for each cell has 10 rows because it is Hermite
  // so, this gives the number of columns per cell
  const int numCols = A.size() / JVol.numCells() / 10;
  
  
  Array<double> cellVertH;
  getVertexH( mesh , cellLIDs , cellVertH );
  
  for (int i=0;i<JVol.numCells();i++) 
    {
      const int cell_start = i * numCols * 10;
      const double *invJ = JVol.jVals(i);
      
      for (int j=0;j<numCols;j++) 
	{
	  const int col_start = cell_start + 10 * j;
	  const double a1 = A[col_start + 1];
	  const double a2 = A[col_start + 2];
	  const double a4 = A[col_start + 4];
	  const double a5 = A[col_start + 5];
	  const double a7 = A[col_start + 7];
	  const double a8 = A[col_start + 8];
	  A[col_start+1] = (invJ[0]*a1 + invJ[2]*a2)/cellVertH[3*i];
	  A[col_start+2] = (invJ[1]*a1 + invJ[3]*a2)/cellVertH[3*i];
	  A[col_start+4] = (invJ[0]*a4 + invJ[2]*a5)/cellVertH[3*i+1];
	  A[col_start+5] = (invJ[1]*a4 + invJ[3]*a5)/cellVertH[3*i+1];
	  A[col_start+7] = (invJ[0]*a7 + invJ[2]*a8)/cellVertH[3*i+2];
	  A[col_start+8] = (invJ[1]*a7 + invJ[3]*a8)/cellVertH[3*i+2];
	}
    }
  
}

// cellVertexH is 3 * cellLIDs.size and has the appropriate h for each
// vertex of each cell
void CubicHermite::getVertexH( const Mesh &mesh, 
			       const Array<int> &cellLIDs ,
			       Array<double> &cellVertexH ) const
{
  cellVertexH.resize(3*cellLIDs.size());
  const int N = mesh.numCells(2);
  const double h = 1.0 / sqrtl( N );
  for (int i=0;i<cellVertexH.size();i++) 
    {
      cellVertexH[i] = h;
    }
//     // now I get the vertices that are included in the batch of cells


//     map<int,double> vert_h;

//     Array<int> cellVertLIDs(3);
//     Array<int> vertOrient(3);
//     const int cellDim = mesh.spatialDim();
//     // set up map structure
//     for (int i=0;i<cellLIDs.size();i++) 
//       {
// 	const int cellLID = cellLIDs[i];
// 	mesh.getFacetArray( cellDim , cellLID , 0 , cellVertLIDs , vertOrient );
// 	for (int j=0;j<3;j++)
// 	  {
// 	    map<int,double>::iterator vert_h_cur = vert_h.find(cellVertLIDs[j]);
// 	    if (vert_h_cur == vert_h.end())
// 	      {
// 		vert_h[cellVertLIDs[j]] = 0.0;
// 	      }
// 	  }
//       }
    
//     // fill in numbers
//     for (map<int,double>::iterator it=vert_h.begin();it!=vert_h.end();++it)
//       {
// 	const int vert = it->first;
// 	it->second = 0.0;
// 	// loop over maximal cofacets of each vertex
// 	Array<int> incident_cells;
// 	Array<double> cell_diams;
// 	mesh.getCofacets( 0 , vert , cellDim , incident_cells );
// 	mesh.getCellDiameters( 2 , incident_cells , cell_diams );
// 	for (int i=0;i<cell_diams.size();i++) 
// 	  {
// 	    it->second += cell_diams[i];
// 	  }
// 	it->second /= cell_diams.size();

//       }

//     cellVertexH.resize( 3 * cellLIDs.size() );

//     for (int i=0;i<cellLIDs.size();i++) 
//       {
// 	const int cellLID = cellLIDs[i];
// 	mesh.getFacetArray( cellDim , cellLID , 0 , cellVertLIDs , vertOrient );
// 	for (int j=0;j<3;j++) 
// 	  {
// 	    cellVertexH[3*i+j] = 1./vert_h[cellVertLIDs[j]];
// 	  }
//       }




    return;

}
