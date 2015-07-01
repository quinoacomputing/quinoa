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

#include "SundanceLegendre.hpp"
#include "SundanceADReal.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;

// ctor
Legendre::Legendre(int order):order_(order)
{

	// set the nr DOFs

	if (order_ > 1) nrDOF_edge_ = order_ - 1;
	else nrDOF_edge_ = 0;

	if (order_ > 3) nrDOF_face_ = ((order_-2)*(order_-3))/2;
	else nrDOF_face_ = 0;

	nrDOF_brick_ = 0;

	//SUNDANCE_OUT( true , "Legendre ctor nrDOF_edge_:" << nrDOF_edge_ << " , nrDOF_face_:" << nrDOF_face_ << " , nrDOF_brick_:"<<nrDOF_brick_);
}

bool Legendre::supportsCellTypePair(
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
    case QuadCell:
      switch(cellType)
      {
        case QuadCell:
        case LineCell:
        case PointCell:
          return true;
        default:
          return false;
      }
     case BrickCell:
       switch(cellType)
       {
         case BrickCell:
         case QuadCell:
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

void Legendre::print(std::ostream& os) const
{
  os << "Legendre";
}

int Legendre::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
      switch(cellType)
      {
        case PointCell:
          return 1;
        case LineCell:
          return nrDOF_edge_;
        case QuadCell:
            return nrDOF_face_;
        case BrickCell:
            return nrDOF_brick_;
        default:
            TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument , "illegal combination of cell type and maximal cell type" );
            return -1;
      }
}

Array<int> Legendre::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

void Legendre::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{

  typedef Array<int> Aint;

  //switch(order_)

  switch(cellType)
  {
    case PointCell:
    {
        dofs.resize(1);
        dofs[0] = tuple<Aint>(tuple(0));
        return;
    }
    break;
    case LineCell:
    {
        dofs.resize(2);
        dofs[0] = tuple<Aint>(tuple(0), tuple(1));
        if (nrDOF_edge_>0)
        {
        	dofs[1] = tuple<Aint>(makeRange(2, 2+nrDOF_edge_-1));
        }
        else
        {
        	dofs[1] = tuple(Array<int>());
        }
      return;
    }
    break;
    case QuadCell:
    {
    	int offs = 0;
        dofs.resize(3);
        // dofs[0] are DOFs at Points
        dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3));
        offs = 4;
        if (nrDOF_edge_>0)
        {
        	dofs[1].resize(4);
        	dofs[1][0] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        	dofs[1][1] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        	dofs[1][2] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        	dofs[1][3] = makeRange(offs, offs+nrDOF_edge_-1);  offs += nrDOF_edge_;
        }
        else
        {
        	dofs[1] = tuple(Array<int>(), Array<int>(),
	                 Array<int>(), Array<int>());
        }

        if (nrDOF_edge_>0)
        {
        	dofs[2].resize(1);
        	dofs[2][0] = makeRange(offs, offs+nrDOF_face_-1);  offs += nrDOF_face_;
        }
        else
        {
        	dofs[2] = tuple(Array<int>());
        }
        //SUNDANCE_OUT( true , "Legendre::getReferenceDOFs offs:" << offs );
    }
    break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
        << cellType << " not implemented in Legendre basis");
  }
}


void Legendre::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    std::runtime_error,
    "cannot evaluate spatial derivative " << sds << " on Legendre basis");
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
    case QuadCell:
      for (int i=0; i<pts.length(); i++)
      {
        evalOnQuad(pts[i], deriv, result[0][i]);
      }
      return;
    case BrickCell:
      for (int i=0; i<pts.length(); i++)
      {
        evalOnBrick(pts[i], deriv, result[0][i]);
      }
      return;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Legendre::refEval() unimplemented for cell type "
        << cellType);

  }
}

/* ---------- evaluation on different cell types -------------- */

void Legendre::evalOnLine(const Point& pt,
  const MultiIndex& deriv,
  Array<double>& result) const
{
  result.resize(2+nrDOF_edge_);
  ADReal x = ADReal(pt[0],0,1);

  Array<ADReal> refAll(7);

  refAll[0] = 1-x;
  refAll[1] = x;
  refAll[2] = 2.44948974278318 * ( (2*x-1)*(2*x-1) - 1 ) / 4;
  refAll[3] = 3.16227766016838 * ( (2*x-1)*(2*x-1) - 1 ) * (2*x-1) / 4;
  refAll[4] = 3.74165738677394 * ( 5*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 6*(2*x-1)*(2*x-1) + 1) / 16;
  refAll[5] = 4.24264068711929 * (2*x-1) * (7*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 10*(2*x-1)*(2*x-1) + 3) / 16;
  refAll[6] = 4.69041575982343 * (21*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) -
		                          35*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) + 15*(2*x-1)*(2*x-1) - 1) / 32;


  for (int i=0; i<result.length(); i++)
  {
    if (deriv.order()==0) 
    {
      result[i] = refAll[i].value();
    }
    else 
    {
      result[i] = refAll[i].gradient()[0];
    }
  }  
  //SUNDANCE_OUT( true , "Legendre::evalOnLine result.length():" << result.length() );
  return;
}

void Legendre::evalOnQuad(const Point& pt,
  const MultiIndex& deriv,
  Array<double>& result) const

{
  result.resize( 4 + 4*nrDOF_edge_ + nrDOF_face_);
  ADReal x = ADReal(pt[0], 0, 2);
  ADReal y = ADReal(pt[1], 1, 2);
  ADReal one(1.0, 2);
  
  Array<ADReal> refAllx(7);
  Array<ADReal> refAlly(7);

  refAllx[0] = 1-x;
  refAllx[1] = x;
  refAllx[2] = 2.44948974278318 * ( (2*x-1)*(2*x-1) - 1 ) / 4;
  refAllx[3] = 3.16227766016838 * ( (2*x-1)*(2*x-1) - 1 ) * (2*x-1) / 4;
  refAllx[4] = 3.74165738677394 * ( 5*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 6*(2*x-1)*(2*x-1) + 1) / 16;
  refAllx[5] = 4.24264068711929 * (2*x-1) * (7*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) - 10*(2*x-1)*(2*x-1) + 3) / 16;
  refAllx[6] = 4.69041575982343 * (21*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) -
		                          35*(2*x-1)*(2*x-1)*(2*x-1)*(2*x-1) + 15*(2*x-1)*(2*x-1) - 1) / 32;

  refAlly[0] = 1-y;
  refAlly[1] = y;
  refAlly[2] = 2.44948974278318 * ( (2*y-1)*(2*y-1) - 1 ) / 4;
  refAlly[3] = 3.16227766016838 * ( (2*y-1)*(2*y-1) - 1 ) * (2*y-1) / 4;
  refAlly[4] = 3.74165738677394 * ( 5*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) - 6*(2*y-1)*(2*y-1) + 1) / 16;
  refAlly[5] = 4.24264068711929 * (2*y-1) * (7*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) - 10*(2*y-1)*(2*y-1) + 3) / 16;
  refAlly[6] = 4.69041575982343 * (21*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) -
		                          35*(2*y-1)*(2*y-1)*(2*y-1)*(2*y-1) + 15*(2*y-1)*(2*y-1) - 1) / 32;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
    << y.value());

  int sideIndex[4][2] = { {0,0} , {1,0} , {0,1} , {1,1}};
  int edgeIndex[4]    = { 0 , 1 , 1 , 0};
  int edgeMultI[4]    = { 0 , 0 , 1 , 1};
  int offs = 0;
  Array<ADReal> tmp(4 + 4*nrDOF_edge_ + nrDOF_face_);

  // loop over vertexes
  for (int i=0; i < 4 ; i++, offs++){
	  tmp[offs] = refAllx[sideIndex[i][0]] * refAlly[sideIndex[i][1]];
  }

  // loop over edges
  for (int i=0; i < 4 ; i++){
	  // loop over each DOF on the edge
	  if (edgeIndex[i] == 0){
		  for (int j = 0 ; j < nrDOF_edge_ ; j++ , offs++){
			  tmp[offs] = refAllx[2+j] * refAlly[edgeMultI[i]];
		  }
	  }
	  else
	  {
		  for (int j = 0 ; j < nrDOF_edge_ ; j++ , offs++){
			  tmp[offs] = refAllx[edgeMultI[i]] * refAlly[2+j];
		  }
	  }
  }

  // loop over all internal DOFs
  if ( nrDOF_face_ > 0 ){
	  // loop for each hierarchical layer
	  for (int hierarch = 0 ; hierarch < order_ - 3 ; hierarch++)
	  {
		  //SUNDANCE_OUT( true , "Legendre::evalOnQuad hierarch:" << hierarch );
		  // for each layer add the basis function
		  for (int i=0 ; i < hierarch+1 ; i++ , offs++)
		  {
			  //SUNDANCE_OUT( true , "Legendre::evalOnQuad offs:" << offs << " 2+i:" << 2+i << " , 2+(hierarch-1-i):" << 2+(hierarch-i));
			  tmp[offs] = refAllx[2+i] * refAlly[2+(hierarch-i)];
		  }
	  }
  }

  // compute the results
  for (int i=0; i<result.length(); i++)
  {
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
  //SUNDANCE_OUT( true , "Legendre::evalOnQuad result.length():" << result.length() );
}

void Legendre::evalOnBrick(const Point& pt,
  const MultiIndex& deriv,
  Array<double>& result) const
{
  ADReal x = ADReal(pt[0], 0, 3);
  ADReal y = ADReal(pt[1], 1, 3);
  ADReal z = ADReal(pt[2], 2, 3);
  ADReal one(1.0, 3);
  
  TEUCHOS_TEST_FOR_EXCEPT(true);
}
