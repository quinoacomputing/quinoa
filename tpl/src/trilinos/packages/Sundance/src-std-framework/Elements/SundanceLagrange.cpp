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

#include "SundanceLagrange.hpp"
#include "SundanceADReal.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


Lagrange::Lagrange(int order)
  : order_(order) , doFInfromationCalculated_(false)
{
TEUCHOS_TEST_FOR_EXCEPTION(order < 0, std::runtime_error,
                     "invalid polynomial order=" << order
                     << " in Lagrange ctor");
}

bool Lagrange::supportsCellTypePair(
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
    case QuadCell:
        switch(cellType){
          case QuadCell:
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

void Lagrange::print(std::ostream& os) const 
{
  os << "Lagrange(" << order_ << ")";
}

int Lagrange::nReferenceDOFsWithoutFacets(
  const CellType& maximalCellType,
  const CellType& cellType
  ) const
{
  if (order_==0)
  {
    if (maximalCellType != cellType) return 0;
    return 1;
  }

  switch(cellType)
  {
    case PointCell:
      return 1;
    case LineCell:
      return order_-1;
    case TriangleCell:
      if (order_ < 3) return 0;
      if (order_ == 3) return 1;
      TEUCHOS_TEST_FOR_EXCEPTION(order_>3, std::runtime_error, 
        "Lagrange order > 3 not implemented "
        "for triangle cells");
      return 0;
    case QuadCell:
      if (order_==1) return 0;
      if (order_==2) return 1;
      TEUCHOS_TEST_FOR_EXCEPTION(order_>2, std::runtime_error, 
        "Lagrange order > 2 not implemented "
        "for quad cells");
    case TetCell:
      if (order_<=2) return 0;
      TEUCHOS_TEST_FOR_EXCEPTION(order_>2, std::runtime_error, 
        "Lagrange order > 2 not implemented "
        "for tet cells");
    case BrickCell:
      if (order_<=1) return 0;
      if (order_==2) return 1;
      TEUCHOS_TEST_FOR_EXCEPTION(order_>2, std::runtime_error, 
        "Lagrange order > 2 not implemented "
        "for brick cells");
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
        << cellType << " not implemented in Lagrange basis");
      return -1; // -Wall
  }
}

void Lagrange::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  typedef Array<int> Aint;

  if (order_==0)
  {
    int dim = dimension(cellType);
    dofs.resize(dim+1);
    for (int d=0; d<dim; d++)
    {
      Array<Array<int> > dd;
      for (int f=0; f<numFacets(cellType, d); f++) dd.append(Array<int>());
      dofs[d] = dd;
    }
    if (cellType==maximalCellType) dofs[dim] = tuple<Aint>(tuple(0));
    else dofs[dim] = tuple<Aint>(Array<int>());
    return;
  }

  switch(cellType)
    {
    case PointCell:
      dofs.resize(1);
      dofs[0] = tuple<Aint>(tuple(0));
      return;
    case LineCell:
      dofs.resize(2);
      dofs[0] = tuple<Aint>(tuple(0), tuple(1));
      dofs[1] = tuple<Aint>(makeRange(2, order()));
      return;
    case TriangleCell:
      {
        int n = order()-1;
        dofs.resize(3);
        dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2));
        dofs[1] = tuple<Aint>(makeRange(3,2+n), 
                        makeRange(3+n, 2+2*n),
                        makeRange(3+2*n, 2+3*n));
        if (order()<3)
          {
            dofs[2] = tuple(Array<int>());
          }
        else 
          {
            dofs[2] = tuple<Aint>(makeRange(3+3*n, 3+3*n));
          }
        return;
      }
    case QuadCell:{ // TODO: ask Kevin about Q0 (P0)
         dofs.resize(3);
    	 // dofs[0] are DOFs at Points
    	 dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3));
    	 if (order() == 2)
    	 {
    	   // dofs[1] are the DOFs on the line
    	   dofs[1] = tuple<Aint>(tuple(4), tuple(5), tuple(6), tuple(7));
    	   //dofs[2] are DOFs inside the Cell

    	   //dofs[2] = tuple(Array<int>());
           dofs[2] = tuple<Aint>(tuple(8));
    	 } else {
    	   dofs[1] = tuple(Array<int>(), Array<int>(),
    		                 Array<int>(), Array<int>());
           dofs[2] = tuple(Array<int>());
    	 }
         return;
      }
    case TetCell:
      {
        dofs.resize(4);
        dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3));
        if (order() == 2)
          {
            dofs[1] = tuple<Aint>(tuple(4), tuple(5), tuple(6), 
                            tuple(7), tuple(8), tuple(9));
          }
        else
          {
            dofs[1] = tuple(Array<int>(), Array<int>(),
                            Array<int>(), Array<int>(), 
                            Array<int>(), Array<int>());
          }
        dofs[2] = tuple(Array<int>(), Array<int>(), 
                        Array<int>(), Array<int>());
        dofs[3] = tuple(Array<int>());
        return;
      }
    case BrickCell:{
         dofs.resize(4);
    	 // dofs[0] are DOFs at Points
    	 dofs[0] = tuple<Aint>(tuple(0), tuple(1), tuple(2), tuple(3),
    			               tuple(4), tuple(5), tuple(6), tuple(7));
    	 if (order() == 2)
    	 {
    	   // dofs[1] are the DOFs on the line
    	   dofs[1] = tuple<Aint>( tuple(8), tuple(9) , tuple(10), tuple(11),
    			                 tuple(12), tuple(13), tuple(14), tuple(15),
    			                 tuple(16), tuple(17), tuple(18), tuple(19));
    	   //dofs[2] are DOFs on the faces
           dofs[2] = tuple<Aint>(tuple(20), tuple(21), tuple(22), tuple(23),
        		           tuple(24), tuple(25) );
           dofs[3] = tuple<Aint>(tuple(26));
    	 } else {
    	   dofs[1] = tuple(Array<int>(), Array<int>(), Array<int>(), Array<int>(),
    			           Array<int>(), Array<int>(), Array<int>(), Array<int>(),
    			           Array<int>(), Array<int>(), Array<int>(), Array<int>());
           dofs[2] = tuple(Array<int>(), Array<int>(), Array<int>(), Array<int>(),
                           Array<int>(), Array<int>() );
           dofs[3] = tuple(Array<int>());
    	 }
         return;
      }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
                         << cellType << " not implemented in Lagrange basis");
    }
}



Array<int> Lagrange::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

void Lagrange::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    std::runtime_error,
    "cannot evaluate spatial derivative " << sds << " on Lagrange basis");
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
    case QuadCell:
      for (int i=0; i<pts.length(); i++)
       {
          evalOnquad(pts[i], deriv, result[0][i]);
       }
       return;
    case TetCell:
      for (int i=0; i<pts.length(); i++)
        {
          evalOnTet(pts[i], deriv, result[0][i]);
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
                         "Lagrange::refEval() unimplemented for cell type "
                         << cellType);

    }
}

/* ---------- evaluation on different cell types -------------- */

void Lagrange::evalOnLine(const Point& pt, 
													const MultiIndex& deriv,
													Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 1);
	ADReal one(1.0, 1);
	
	result.resize(order()+1);
	Array<ADReal> tmp(result.length());
  Array<double> x0(order()+1);

  if (order_ == 0)
    {
      tmp[0] = one;
    }
  else
    {
      x0[0] = 0.0;
      x0[1] = 1.0;
      for (int i=0; i<order_-1; i++)
        {
          x0[i+2] = (i+1.0)/order_;
        }

      for (int i=0; i<=order_; i++)
        {
          tmp[i] = one;
          for (int j=0; j<=order_; j++)
            {
              if (i==j) continue;
              tmp[i] *= (x - x0[j])/(x0[i]-x0[j]);
            }
        }
    }
  /*
	switch(order_)
		{
		case 0:
			tmp[0] = one;
			break;
		case 1:
			tmp[0] = 1.0 - x;
			tmp[1] = x;
			break;
		case 2:
			tmp[0] = 2.0*(1.0-x)*(0.5-x);
			tmp[1] = 2.0*x*(x-0.5);
			tmp[2] = 4.0*x*(1.0-x);
			break;
    case 3:
      tmp[0] = 9.0/2.0 * (1.0 - x) * (1.0/3.0 - x) * (2.0/3.0 - x);
      tmp[1] = 9.0/2.0 * x * (x - 1.0/3.0) * (x - 2.0/3.0);
      tmp[2] = 27.0/2.0 * x * (1.0 - x) * (2.0/3.0 - x);
      tmp[3] = 27.0/2.0 * x * (1.0 - x) * (x - 1.0/3.0);

      break;
		default:
			SUNDANCE_ERROR("Lagrange::evalOnLine polynomial order > 2 has not been"
                     " implemented");
		}
  */
	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else result[i] = tmp[i].gradient()[0];
		}
}

void Lagrange::evalOnTriangle(const Point& pt, 
															const MultiIndex& deriv,
															Array<double>& result) const



{
	ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
               << y.value());

	switch(order_)
		{
		case 0:
      result.resize(1);
      tmp.resize(1);
			tmp[0] = one;
			break;
		case 1:
      result.resize(3);
      tmp.resize(3);
			tmp[0] = 1.0 - x - y;
			tmp[1] = x;
			tmp[2] = y;
			break;
		case 2:
      result.resize(6);
      tmp.resize(6);
			tmp[0] = (1.0-x-y)*(1.0-2.0*x-2.0*y);
			tmp[1] = 2.0*x*(x-0.5);
			tmp[2] = 2.0*y*(y-0.5);
			tmp[5] = 4.0*x*(1.0-x-y); 
			tmp[3] = 4.0*x*y;
			tmp[4] = 4.0*y*(1.0-x-y);
			break;
    case 3:
      result.resize(10);
      tmp.resize(10);
      {
        ADReal q1 = 1.0 - x - y;
        ADReal q2 = x;
        ADReal q3 = y;
        tmp[0] = 9.0/2.0 * q1 * (q1 - 2.0/3.0) * (q1 - 1.0/3.0);
        tmp[1] = 9.0/2.0 * q2 * (q2 - 2.0/3.0) * (q2 - 1.0/3.0);
        tmp[2] = 9.0/2.0 * q3 * (q3 - 2.0/3.0) * (q3 - 1.0/3.0);
        tmp[7] = 27.0/2.0 * q1 * q2   * (q1 - 1.0/3.0);
        tmp[8] = 27.0/2.0 * q1 * q2   * (q2 - 1.0/3.0);
        tmp[3] = 27.0/2.0 * q2 * q3   * (q2 - 1.0/3.0);
        tmp[4] = 27.0/2.0 * q2 * q3   * (q3 - 1.0/3.0);
        tmp[6] = 27.0/2.0 * q3 * q1   * (q3 - 1.0/3.0);
        tmp[5] = 27.0/2.0 * q3 * q1   * (q1 - 1.0/3.0);
        tmp[9] = 27.0 * q1 * q2 * q3;
      }
      break;
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR("Lagrange::evalOnTriangle poly order > 2");
#else
			SUNDANCE_ERROR7("Lagrange::evalOnTriangle poly order > 2");
#endif
		}

	for (int i=0; i<tmp.length(); i++)
		{
      SUNDANCE_OUT(this->verb() > 3,
                   "tmp[" << i << "]=" << tmp[i].value() 
                   << " grad=" << tmp[i].gradient());
			if (deriv.order()==0) result[i] = tmp[i].value();
			else 
          result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}

void Lagrange::evalOnquad(const Point& pt,
	    const MultiIndex& deriv,
	    Array<double>& result) const {

	ADReal x = ADReal(pt[0], 0, 2);
	ADReal y = ADReal(pt[1], 1, 2);
	ADReal one(1.0, 2);

  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
               << y.value());

	switch(order_)
		{
		case 0:
            result.resize(1);
            tmp.resize(1);
			tmp[0] = one;
			break;
		case 1:
            result.resize(4);
            tmp.resize(4);
			tmp[0] = (1.0 - x)*(1.0 - y);
			tmp[1] = x*(1.0 - y);
			tmp[2] = (1.0 - x)*y;
			tmp[3] = x*y; //These are the 4 basis functions for the first order
			break;
		case 2: //The basis functions calculated with Octave / Matlab

            result.resize(9);
            tmp.resize(9);
			tmp[0] = 9.0 - 9.0*x - 9.0*y + 9.0*x*y - 6.0*(x*x -1)*(y-1)
					- 6.0*(x-1)*(y*y-1)                     + 4.0*(x*x-1)*(y*y-1);
			tmp[1] = 6.0 - 3.0*x - 6.0*y + 3.0*x*y - 6.0*(x*x -1)*(y-1)
                    - 3.0*(x-1)*(y*y-1) + 1.0*(x+1)*(y*y-1) + 4.0*(x*x-1)*(y*y-1);
			tmp[2] = 6.0 - 6.0*x - 3.0*y + 3.0*x*y - 3.0*(x*x -1)*(y-1) + 1.0*(x*x -1)*(y+1)
					-6.0*(x-1)*(y*y-1) +                    + 4.0*(x*x-1)*(y*y-1);
			tmp[3] = 4.0 - 2.0*x - 2.0*y + 1.0*x*y - 3.0*(x*x -1)*(y-1) + 1.0*(x*x -1)*(y+1)
					- 3.0*(x-1)*(y*y-1) + 1.0*(x+1)*(y*y-1) + 4.0*(x*x-1)*(y*y-1);
			tmp[4] = -12.0 + 12.0*x + 12.0*y - 12.0*x*y + 12.0*(x*x -1)*(y-1)
					+ 8.0*(x-1)*(y*y-1)                     - 8.0*(x*x-1)*(y*y-1);
			tmp[5] = -12.0 + 12.0*x + 12.0*y - 12.0*x*y + 8.0*(x*x -1)*(y-1)
					+ 12.0*(x-1)*(y*y-1)                    - 8.0*(x*x-1)*(y*y-1);
			tmp[6] = -8.0 + 4.0*x + 8.0*y - 4.0*x*y + 8.0*(x*x -1)*(y-1)
					+ 6.0*(x-1)*(y*y-1) - 2.0*(x+1)*(y*y-1) - 8.0*(x*x-1)*(y*y-1);
			tmp[7] = -8.0 + 8.0*x + 4.0*y - 4.0*x*y + 6.0*(x*x -1)*(y-1) - 2.0*(x*x -1)*(y+1)
					+ 8.0*(x-1)*(y*y-1)                     - 8.0*(x*x-1)*(y*y-1);
			tmp[8] = 16.0 - 16.0*x - 16.0*y + 16.0*x*y - 16.0*(x*x -1)*(y-1)
					- 16.0*(x-1)*(y*y-1)                    + 16.0*(x*x-1)*(y*y-1);
			break;
		default:{}
		}

	 for (int i=0; i<tmp.length(); i++) {
		 SUNDANCE_OUT(this->verb() > 3 ,
                   "tmp[" << i << "]=" << tmp[i].value()
                   << " grad=" << tmp[i].gradient());
		if (deriv.order()==0)
				result[i] = tmp[i].value();
		else
                result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
	}
}

void Lagrange::evalOnTet(const Point& pt, 
												 const MultiIndex& deriv,
												 Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 3);
	ADReal y = ADReal(pt[1], 1, 3);
	ADReal z = ADReal(pt[2], 2, 3);
	ADReal one(1.0, 3);


	Array<ADReal> tmp(result.length());

	switch(order_)
		{
		case 0:
      tmp.resize(1);
      result.resize(1);
			tmp[0] = one;
			break;
		case 1:
      result.resize(4);
      tmp.resize(4);
			tmp[0] = 1.0 - x - y - z;
			tmp[1] = x;
			tmp[2] = y;
			tmp[3] = z;
			break;
		case 2:
      result.resize(10);
      tmp.resize(10);
			tmp[0] = (1.0-x-y-z)*(1.0-2.0*x-2.0*y-2.0*z);
			tmp[1] = 2.0*x*(x-0.5);
			tmp[2] = 2.0*y*(y-0.5);
			tmp[3] = 2.0*z*(z-0.5);
			tmp[9] = 4.0*x*(1.0-x-y-z);
			tmp[6] = 4.0*x*y;
			tmp[8] = 4.0*y*(1.0-x-y-z);
			tmp[7] = 4.0*z*(1.0-x-y-z);
			tmp[5] = 4.0*x*z;
			tmp[4] = 4.0*y*z;
			break;
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR("Lagrange::evalOnTet poly order > 2");
#else
			SUNDANCE_ERROR7("Lagrange::evalOnTet poly order > 2");
#endif
		}

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else 
				result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}

void Lagrange::evalOnBrick(const Point& pt,
												 const MultiIndex& deriv,
												 Array<double>& result) const
{
	ADReal x = ADReal(pt[0], 0, 3);
	ADReal y = ADReal(pt[1], 1, 3);
	ADReal z = ADReal(pt[2], 2, 3);
	ADReal one(1.0, 3);


	Array<ADReal> tmp(result.length());

	switch(order_)
		{
		case 0:
            tmp.resize(1);
            result.resize(1);
			tmp[0] = one;
			break;
		case 1: // tri linear basis functions for bricks
            result.resize(8);
            tmp.resize(8);
			tmp[0] = (1.0 - x)*(1.0 - y)*(1.0 - z);
			tmp[1] = (x)*(1.0 - y)*(1.0 - z);
			tmp[2] = (1.0 - x)*(y)*(1.0 - z);
			tmp[3] = (x)*(y)*(1.0 - z);
			tmp[4] = (1.0 - x)*(1.0 - y)*(z);
			tmp[5] = (x)*(1.0 - y)*(z);
			tmp[6] = (1.0 - x)*(y)*(z);
			tmp[7] = (x)*(y)*(z);
			break;
		case 2: // second order function in 3D
	        result.resize(27);
	        tmp.resize(27);
                   // 0.0 + 0.0*x + 0.0*y + 0.0*z + 0.0*x*y + 0.0*x*z + 0.0*y*z + 0.0*x*y*z
                   //+ 0.0*(x*x - 1)*(y-1)*(z-1) + 0.0*(x*x - 1)*(y-1)*(z+1) + 0.0*(x*x - 1)*(y+1)*(z-1) + 0.0*(x*x - 1)*(y+1)*(z+1)
                   //+ 0.0*(y*y - 1)*(x-1)*(z-1) + 0.0*(y*y - 1)*(x-1)*(z+1) + 0.0*(y*y - 1)*(x+1)*(z-1) + 0.0*(y*y - 1)*(x+1)*(z+1)
                   //+ 0.0*(z*z - 1)*(x-1)*(y-1) + 0.0*(z*z - 1)*(x-1)*(y+1) + 0.0*(z*z - 1)*(x+1)*(y-1) + 0.0*(z*z - 1)*(x+1)*(y+1)
                   //+ 0.0*(x*x - 1)*(y*y - 1)*(z-1) + 0.0*(x*x - 1)*(y*y - 1)*(z+1)
                   //+ 0.0*(x*x - 1)*(z*z - 1)*(y-1) + 0.0*(x*x - 1)*(z*z - 1)*(y+1)
                   //+ 0.0*(y*y - 1)*(z*z - 1)*(x-1) + 0.0*(y*y - 1)*(z*z - 1)*(x+1)
                   //+ 0.0*(x*x - 1)*(y*y - 1)*(z*z - 1);


            tmp[0] = 27.0 - 27.0*x - 27.0*y - 27.0*z + 27.0*x*y + 27.0*x*z + 27.0*y*z - 27.0*x*y*z
            + 18.0*(x*x - 1)*(y-1)*(z-1)
            + 18.0*(y*y - 1)*(x-1)*(z-1)
            + 18.0*(z*z - 1)*(x-1)*(y-1)
            - 12.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 12.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 12.0*(y*y - 1)*(z*z - 1)*(x-1)
            + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[1] = 18.0 - 9.0*x - 18.0*y - 18.0*z + 9.0*x*y + 9.0*x*z + 18.0*y*z - 9.0*x*y*z
            + 18.0*(x*x - 1)*(y-1)*(z-1)
            + 9.0*(y*y - 1)*(x-1)*(z-1) - 3.0*(y*y - 1)*(x+1)*(z-1)
            + 9.0*(z*z - 1)*(x-1)*(y-1) - 3.0*(z*z - 1)*(x+1)*(y-1)
            - 12.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 12.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 6.0*(y*y - 1)*(z*z - 1)*(x-1) + 2.0*(y*y - 1)*(z*z - 1)*(x+1)
            + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[2] = 18.0 - 18.0*x - 9.0*y - 18.0*z + 9.0*x*y + 18.0*x*z + 9.0*y*z - 9.0*x*y*z
             + 9.0*(x*x - 1)*(y-1)*(z-1) - 3.0*(x*x - 1)*(y+1)*(z-1)
             + 18.0*(y*y - 1)*(x-1)*(z-1)
             + 9.0*(z*z - 1)*(x-1)*(y-1) - 3.0*(z*z - 1)*(x-1)*(y+1)
             - 12.0*(x*x - 1)*(y*y - 1)*(z-1)
             - 6.0*(x*x - 1)*(z*z - 1)*(y-1) + 2.0*(x*x - 1)*(z*z - 1)*(y+1)
             - 12.0*(y*y - 1)*(z*z - 1)*(x-1)
             + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[3] =  12.0 - 6.0*x - 6.0*y - 12.0*z + 3.0*x*y + 6.0*x*z + 6.0*y*z - 3.0*x*y*z
            + 9.0*(x*x - 1)*(y-1)*(z-1) - 3.0*(x*x - 1)*(y+1)*(z-1)
            + 9.0*(y*y - 1)*(x-1)*(z-1) - 3.0*(y*y - 1)*(x+1)*(z-1)
            + 4.5*(z*z - 1)*(x-1)*(y-1) - 1.5*(z*z - 1)*(x-1)*(y+1) - 1.5*(z*z - 1)*(x+1)*(y-1) + 0.5*(z*z - 1)*(x+1)*(y+1)
            - 12.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 6.0*(x*x - 1)*(z*z - 1)*(y-1) + 2.0*(x*x - 1)*(z*z - 1)*(y+1)
            - 6.0*(y*y - 1)*(z*z - 1)*(x-1) + 2.0*(y*y - 1)*(z*z - 1)*(x+1)
            + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[4] = 18.0 - 18.0*x - 18.0*y - 9.0*z + 18.0*x*y + 9.0*x*z + 9.0*y*z - 9.0*x*y*z
             + 9.0*(x*x - 1)*(y-1)*(z-1) - 3.0*(x*x - 1)*(y-1)*(z+1)
             + 9.0*(y*y - 1)*(x-1)*(z-1) - 3.0*(y*y - 1)*(x-1)*(z+1)
             + 18.0*(z*z - 1)*(x-1)*(y-1)
             - 6.0*(x*x - 1)*(y*y - 1)*(z-1) + 2.0*(x*x - 1)*(y*y - 1)*(z+1)
             - 12.0*(x*x - 1)*(z*z - 1)*(y-1)
             - 12.0*(y*y - 1)*(z*z - 1)*(x-1)
             + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[5] =  12.0 - 6.0*x - 12.0*y - 6.0*z + 6.0*x*y + 3.0*x*z + 6.0*y*z - 3.0*x*y*z
            + 9.0*(x*x - 1)*(y-1)*(z-1) - 3.0*(x*x - 1)*(y-1)*(z+1)
            + 4.5*(y*y - 1)*(x-1)*(z-1) - 1.5*(y*y - 1)*(x-1)*(z+1) - 1.5*(y*y - 1)*(x+1)*(z-1) + 0.5*(y*y - 1)*(x+1)*(z+1)
            + 9.0*(z*z - 1)*(x-1)*(y-1) - 3.0*(z*z - 1)*(x+1)*(y-1)
            - 6.0*(x*x - 1)*(y*y - 1)*(z-1) + 2.0*(x*x - 1)*(y*y - 1)*(z+1)
            - 12.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 6.0*(y*y - 1)*(z*z - 1)*(x-1) + 2.0*(y*y - 1)*(z*z - 1)*(x+1)
            + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[6] = 12.0 - 12.0*x - 6.0*y - 6.0*z + 6.0*x*y + 6.0*x*z + 3.0*y*z - 3.0*x*y*z
            + 4.5*(x*x - 1)*(y-1)*(z-1) - 1.5*(x*x - 1)*(y-1)*(z+1) - 1.5*(x*x - 1)*(y+1)*(z-1) + 0.5*(x*x - 1)*(y+1)*(z+1)
            + 9.0*(y*y - 1)*(x-1)*(z-1) - 3.0*(y*y - 1)*(x-1)*(z+1)
            + 9.0*(z*z - 1)*(x-1)*(y-1) - 3.0*(z*z - 1)*(x-1)*(y+1)
            - 6.0*(x*x - 1)*(y*y - 1)*(z-1) + 2.0*(x*x - 1)*(y*y - 1)*(z+1)
            - 6.0*(x*x - 1)*(z*z - 1)*(y-1) + 2.0*(x*x - 1)*(z*z - 1)*(y+1)
            - 12.0*(y*y - 1)*(z*z - 1)*(x-1)
            + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[7] = 8.0 - 4.0*x - 4.0*y - 4.0*z + 2.0*x*y + 2.0*x*z + 2.0*y*z - 1.0*x*y*z
             + 4.5*(x*x - 1)*(y-1)*(z-1) - 1.5*(x*x - 1)*(y-1)*(z+1) - 1.5*(x*x - 1)*(y+1)*(z-1) + 0.5*(x*x - 1)*(y+1)*(z+1)
             + 4.5*(y*y - 1)*(x-1)*(z-1) - 1.5*(y*y - 1)*(x-1)*(z+1) - 1.5*(y*y - 1)*(x+1)*(z-1) + 0.5*(y*y - 1)*(x+1)*(z+1)
             + 4.5*(z*z - 1)*(x-1)*(y-1) - 1.5*(z*z - 1)*(x-1)*(y+1) - 1.5*(z*z - 1)*(x+1)*(y-1) + 0.5*(z*z - 1)*(x+1)*(y+1)
             - 6.0*(x*x - 1)*(y*y - 1)*(z-1) + 2.0*(x*x - 1)*(y*y - 1)*(z+1)
             - 6.0*(x*x - 1)*(z*z - 1)*(y-1) + 2.0*(x*x - 1)*(z*z - 1)*(y+1)
             - 6.0*(y*y - 1)*(z*z - 1)*(x-1) + 2.0*(y*y - 1)*(z*z - 1)*(x+1)
             + 8.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[8] =  -36.0 + 36.0*x + 36.0*y + 36.0*z - 36.0*x*y - 36.0*x*z - 36.0*y*z + 36.0*x*y*z
             - 36.0*(x*x - 1)*(y-1)*(z-1)
             - 24.0*(y*y - 1)*(x-1)*(z-1)
             - 24.0*(z*z - 1)*(x-1)*(y-1)
             + 24.0*(x*x - 1)*(y*y - 1)*(z-1)
             + 24.0*(x*x - 1)*(z*z - 1)*(y-1)
             + 16.0*(y*y - 1)*(z*z - 1)*(x-1)
             - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[9] = -36.0 + 36.0*x + 36.0*y + 36.0*z - 36.0*x*y - 36.0*x*z - 36.0*y*z + 36.0*x*y*z
            - 24.0*(x*x - 1)*(y-1)*(z-1)
            - 36.0*(y*y - 1)*(x-1)*(z-1)
            - 24.0*(z*z - 1)*(x-1)*(y-1)
            + 24.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 16.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 24.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[10] = -36.0 + 36.0*x + 36.0*y + 36.0*z - 36.0*x*y - 36.0*x*z - 36.0*y*z + 36.0*x*y*z
            - 24.0*(x*x - 1)*(y-1)*(z-1)
            - 24.0*(y*y - 1)*(x-1)*(z-1)
            - 36.0*(z*z - 1)*(x-1)*(y-1)
            + 16.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 24.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 24.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[11] = -24.0 + 12.0*x + 24.0*y + 24.0*z - 12.0*x*y - 12.0*x*z - 24.0*y*z + 12.0*x*y*z
            - 24.0*(x*x - 1)*(y-1)*(z-1)
            - 18.0*(y*y - 1)*(x-1)*(z-1) + 6.0*(y*y - 1)*(x+1)*(z-1)
            - 12.0*(z*z - 1)*(x-1)*(y-1) + 4.0*(z*z - 1)*(x+1)*(y-1)
            + 24.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 16.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 12.0*(y*y - 1)*(z*z - 1)*(x-1) - 4.0*(y*y - 1)*(z*z - 1)*(x+1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[12] = -24.0 + 12.0*x + 24.0*y + 24.0*z - 12.0*x*y - 12.0*x*z - 24.0*y*z + 12.0*x*y*z
            - 24.0*(x*x - 1)*(y-1)*(z-1)
            - 12.0*(y*y - 1)*(x-1)*(z-1) + 4.0*(y*y - 1)*(x+1)*(z-1)
            - 18.0*(z*z - 1)*(x-1)*(y-1) + 6.0*(z*z - 1)*(x+1)*(y-1)
            + 16.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 24.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 12.0*(y*y - 1)*(z*z - 1)*(x-1) - 4.0*(y*y - 1)*(z*z - 1)*(x+1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[13] = -24.0 + 24.0*x + 12.0*y + 24.0*z - 12.0*x*y - 24.0*x*z - 12.0*y*z + 12.0*x*y*z
            - 18.0*(x*x - 1)*(y-1)*(z-1) + 6.0*(x*x - 1)*(y+1)*(z-1)
            - 24.0*(y*y - 1)*(x-1)*(z-1)
            - 12.0*(z*z - 1)*(x-1)*(y-1) + 4.0*(z*z - 1)*(x-1)*(y+1)
            + 24.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 12.0*(x*x - 1)*(z*z - 1)*(y-1) - 4.0*(x*x - 1)*(z*z - 1)*(y+1)
            + 16.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[14] = -24.0 + 24.0*x + 12.0*y + 24.0*z - 12.0*x*y - 24.0*x*z - 12.0*y*z + 12.0*x*y*z
            - 12.0*(x*x - 1)*(y-1)*(z-1) + 4.0*(x*x - 1)*(y+1)*(z-1)
            - 24.0*(y*y - 1)*(x-1)*(z-1)
            - 18.0*(z*z - 1)*(x-1)*(y-1) + 6.0*(z*z - 1)*(x-1)*(y+1)
            + 16.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 12.0*(x*x - 1)*(z*z - 1)*(y-1) - 4.0*(x*x - 1)*(z*z - 1)*(y+1)
            + 24.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[15] = -16.0 + 8.0*x + 8.0*y + 16.0*z - 4.0*x*y - 8.0*x*z - 8.0*y*z + 4.0*x*y*z
            - 12.0*(x*x - 1)*(y-1)*(z-1) + 4.0*(x*x - 1)*(y+1)*(z-1)
            - 12.0*(y*y - 1)*(x-1)*(z-1) + 4.0*(y*y - 1)*(x+1)*(z-1)
            - 9.0*(z*z - 1)*(x-1)*(y-1) + 3.0*(z*z - 1)*(x-1)*(y+1) + 3.0*(z*z - 1)*(x+1)*(y-1) - 1.0*(z*z - 1)*(x+1)*(y+1)
            + 16.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 12.0*(x*x - 1)*(z*z - 1)*(y-1) - 4.0*(x*x - 1)*(z*z - 1)*(y+1)
            + 12.0*(y*y - 1)*(z*z - 1)*(x-1) - 4.0*(y*y - 1)*(z*z - 1)*(x+1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[16] = -24.0 + 24.0*x + 24.0*y + 12.0*z - 24.0*x*y - 12.0*x*z - 12.0*y*z + 12.0*x*y*z
            - 18.0*(x*x - 1)*(y-1)*(z-1) + 6.0*(x*x - 1)*(y-1)*(z+1)
            - 12.0*(y*y - 1)*(x-1)*(z-1) + 4.0*(y*y - 1)*(x-1)*(z+1)
            - 24.0*(z*z - 1)*(x-1)*(y-1)
            + 12.0*(x*x - 1)*(y*y - 1)*(z-1) - 4.0*(x*x - 1)*(y*y - 1)*(z+1)
            + 24.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 16.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[17] = -24.0 + 24.0*x + 24.0*y + 12.0*z - 24.0*x*y - 12.0*x*z - 12.0*y*z + 12.0*x*y*z
            - 12.0*(x*x - 1)*(y-1)*(z-1) + 4.0*(x*x - 1)*(y-1)*(z+1)
            - 18.0*(y*y - 1)*(x-1)*(z-1) + 6.0*(y*y - 1)*(x-1)*(z+1)
            - 24.0*(z*z - 1)*(x-1)*(y-1)
            + 12.0*(x*x - 1)*(y*y - 1)*(z-1) - 4.0*(x*x - 1)*(y*y - 1)*(z+1)
            + 16.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 24.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[18] = -16.0 + 8.0*x + 16.0*y + 8.0*z - 8.0*x*y - 4.0*x*z - 8.0*y*z + 4.0*x*y*z
            - 12.0*(x*x - 1)*(y-1)*(z-1) + 4.0*(x*x - 1)*(y-1)*(z+1)
            - 9.0*(y*y - 1)*(x-1)*(z-1) + 3.0*(y*y - 1)*(x-1)*(z+1) + 3.0*(y*y - 1)*(x+1)*(z-1) - 1.0*(y*y - 1)*(x+1)*(z+1)
            - 12.0*(z*z - 1)*(x-1)*(y-1) + 4.0*(z*z - 1)*(x+1)*(y-1)
            + 12.0*(x*x - 1)*(y*y - 1)*(z-1) - 4.0*(x*x - 1)*(y*y - 1)*(z+1)
            + 16.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 12.0*(y*y - 1)*(z*z - 1)*(x-1) - 4.0*(y*y - 1)*(z*z - 1)*(x+1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[19] = -16.0 + 16.0*x + 8.0*y + 8.0*z - 8.0*x*y - 8.0*x*z - 4.0*y*z + 4.0*x*y*z
            - 9.0*(x*x - 1)*(y-1)*(z-1) + 3.0*(x*x - 1)*(y-1)*(z+1) + 3.0*(x*x - 1)*(y+1)*(z-1) - 1.0*(x*x - 1)*(y+1)*(z+1)
            - 12.0*(y*y - 1)*(x-1)*(z-1) + 4.0*(y*y - 1)*(x-1)*(z+1)
            - 12.0*(z*z - 1)*(x-1)*(y-1) + 4.0*(z*z - 1)*(x-1)*(y+1)
            + 12.0*(x*x - 1)*(y*y - 1)*(z-1) - 4.0*(x*x - 1)*(y*y - 1)*(z+1)
            + 12.0*(x*x - 1)*(z*z - 1)*(y-1) - 4.0*(x*x - 1)*(z*z - 1)*(y+1)
            + 16.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 16.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[20] = 48.0 - 48.0*x - 48.0*y - 48.0*z + 48.0*x*y + 48.0*x*z + 48.0*y*z - 48.0*x*y*z
            + 48.0*(x*x - 1)*(y-1)*(z-1)
            + 48.0*(y*y - 1)*(x-1)*(z-1)
            + 32.0*(z*z - 1)*(x-1)*(y-1)
            - 48.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 32.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 32.0*(y*y - 1)*(z*z - 1)*(x-1)
            + 32.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[21] = 48.0 - 48.0*x - 48.0*y - 48.0*z + 48.0*x*y + 48.0*x*z + 48.0*y*z - 48.0*x*y*z
            + 48.0*(x*x - 1)*(y-1)*(z-1)
            + 32.0*(y*y - 1)*(x-1)*(z-1)
            + 48.0*(z*z - 1)*(x-1)*(y-1)
            - 32.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 48.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 32.0*(y*y - 1)*(z*z - 1)*(x-1)
            + 32.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[22] = 48.0 - 48.0*x - 48.0*y - 48.0*z + 48.0*x*y + 48.0*x*z + 48.0*y*z - 48.0*x*y*z
            + 32.0*(x*x - 1)*(y-1)*(z-1)
            + 48.0*(y*y - 1)*(x-1)*(z-1)
            + 48.0*(z*z - 1)*(x-1)*(y-1)
            - 32.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 32.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 48.0*(y*y - 1)*(z*z - 1)*(x-1)
            + 32.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[23] = 32.0 - 16.0*x - 32.0*y - 32.0*z + 16.0*x*y + 16.0*x*z + 32.0*y*z - 16.0*x*y*z
            + 32.0*(x*x - 1)*(y-1)*(z-1)
            + 24.0*(y*y - 1)*(x-1)*(z-1) - 8.0*(y*y - 1)*(x+1)*(z-1)
            + 24.0*(z*z - 1)*(x-1)*(y-1) - 8.0*(z*z - 1)*(x+1)*(y-1)
            - 32.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 32.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 24.0*(y*y - 1)*(z*z - 1)*(x-1) + 8.0*(y*y - 1)*(z*z - 1)*(x+1)
            + 32.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[24] = 32.0 - 32.0*x - 16.0*y - 32.0*z + 16.0*x*y + 32.0*x*z + 16.0*y*z - 16.0*x*y*z
            + 24.0*(x*x - 1)*(y-1)*(z-1)- 8.0*(x*x - 1)*(y+1)*(z-1)
            + 32.0*(y*y - 1)*(x-1)*(z-1)
            + 24.0*(z*z - 1)*(x-1)*(y-1) - 8.0*(z*z - 1)*(x-1)*(y+1)
            - 32.0*(x*x - 1)*(y*y - 1)*(z-1)
            - 24.0*(x*x - 1)*(z*z - 1)*(y-1) + 8.0*(x*x - 1)*(z*z - 1)*(y+1)
            - 32.0*(y*y - 1)*(z*z - 1)*(x-1)
            + 32.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[25] = 32.0 - 32.0*x - 32.0*y - 16.0*z + 32.0*x*y + 16.0*x*z + 16.0*y*z - 16.0*x*y*z
            + 24.0*(x*x - 1)*(y-1)*(z-1) - 8.0*(x*x - 1)*(y-1)*(z+1)
            + 24.0*(y*y - 1)*(x-1)*(z-1) - 8.0*(y*y - 1)*(x-1)*(z+1)
            + 32.0*(z*z - 1)*(x-1)*(y-1)
            - 24.0*(x*x - 1)*(y*y - 1)*(z-1) + 8.0*(x*x - 1)*(y*y - 1)*(z+1)
            - 32.0*(x*x - 1)*(z*z - 1)*(y-1)
            - 32.0*(y*y - 1)*(z*z - 1)*(x-1)
            + 32.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

            tmp[26] = -64.0 + 64.0*x + 64.0*y + 64.0*z - 64.0*x*y - 64.0*x*z - 64.0*y*z + 64.0*x*y*z
            - 64.0*(x*x - 1)*(y-1)*(z-1)
            - 64.0*(y*y - 1)*(x-1)*(z-1)
            - 64.0*(z*z - 1)*(x-1)*(y-1)
            + 64.0*(x*x - 1)*(y*y - 1)*(z-1)
            + 64.0*(x*x - 1)*(z*z - 1)*(y-1)
            + 64.0*(y*y - 1)*(z*z - 1)*(x-1)
            - 64.0*(x*x - 1)*(y*y - 1)*(z*z - 1);

			break;
		default:
#ifndef TRILINOS_7
			SUNDANCE_ERROR("Lagrange::evalOnBrick poly order > 2");
#else
			SUNDANCE_ERROR7("Lagrange::evalOnBrick poly order > 2");
#endif
		}

	for (int i=0; i<tmp.length(); i++)
		{
			if (deriv.order()==0) result[i] = tmp[i].value();
			else
				result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
		}
}

void Lagrange::getConstrainsForHNDoF(
 	const int indexInParent,
 	const int maxCellDim,
 	const int maxNrChild,
 	const int facetDim,
 	const int facetIndex,
 	const int nodeIndex,
 	Array<int>& localDoFs,
  	Array<int>& parentFacetDim,
  	Array<int>& parentFacetIndex,
  	Array<int>& parentFacetNode,
 	Array<double>& coefs
 	) {

	// the index in the Parent is very important
	// - in case of Lagrange first get the position of the evaluation
	// - then evaluate the refElement at that point
	// - see which local DoFs are not zero and return them

	Array<Array<Array<double> > >    result;
	const SpatialDerivSpecifier      deriv;
	int                              index = 0;

	//setVerb(6);

	localDoFs.resize(0);
	coefs.resize(0);
	parentFacetDim.resize(0);
	parentFacetIndex.resize(0);
	parentFacetNode.resize(0);

	CellType maximalCellType;

	switch (maxCellDim){
	case 2:{ // 2D
		maximalCellType = QuadCell;

		switch (maxNrChild){
		case 9:{ //trisection on quads
			Point dofPos( 0.0 , 0.0 );
            // get the position of the DoF on the refCell
			SUNDANCE_OUT(this->verb() > 3 , ",maxCellDim=" << maxCellDim << ",facetDim="
					<< facetDim << ",facetIndex=" << facetIndex << ",nodeIndex=" << nodeIndex
					<< ",dofPos=" << dofPos);
			returnDoFPositionOnRef( maxCellDim,
			  facetDim, facetIndex, nodeIndex, dofPos);
			dofPos = dofPos/3.0;
			SUNDANCE_OUT(this->verb() > 3 ,  "dofPos=" << dofPos );

			// shitf the position of the DOFs according to the position of the child in the parent cell
			switch (indexInParent){
			  case 0: {/* nothing to add */  break;}
			  case 1: {dofPos[0] += 1.0/3.0; break;}
			  case 2: {dofPos[0] += 2.0/3.0; break;}
			  case 3: {dofPos[0] += 2.0/3.0; dofPos[1] += 1.0/3.0; break;}
			  case 5: {dofPos[1] += 1.0/3.0; break;}
			  case 6: {dofPos[1] += 2.0/3.0; break;}
			  case 7: {dofPos[0] += 1.0/3.0; dofPos[1] += 2.0/3.0; break;}
			  case 8: {dofPos[0] += 2.0/3.0; dofPos[1] += 2.0/3.0; break;}
			}

			SUNDANCE_OUT(this->verb() > 3 , "dofPos=" << dofPos );
			Array<Point> tmp_points(1);
			tmp_points[0] = dofPos;
			refEval( maximalCellType , tmp_points ,
			     deriv,  result , 4);

		   break;}
		case 4:{ // bisection on quads
			Point dofPos( 0.0 , 0.0 );
		    // get the position of the DoF on the refCell
			SUNDANCE_OUT(this->verb() > 3 , ",maxCellDim=" << maxCellDim << ",facetDim="
					<< facetDim << ",facetIndex=" << facetIndex << ",nodeIndex=" << nodeIndex
					<< ",dofPos=" << dofPos);
			returnDoFPositionOnRef( maxCellDim,
			  facetDim, facetIndex, nodeIndex, dofPos);
			dofPos = dofPos/2.0;
			SUNDANCE_OUT(this->verb() > 3 ,  "dofPos=" << dofPos );
			// shitf the position of the DOFs according to the position of the child in the parent cell
			switch (indexInParent){
			  case 0: {/* nothing to add */  break;}
			  case 1: {dofPos[0] += 1.0/2.0; break;}
			  case 2: {dofPos[1] += 1.0/2.0; break;}
			  case 3: {dofPos[0] += 1.0/2.0; dofPos[1] += 1.0/2.0; break;}
			}

			SUNDANCE_OUT(this->verb() > 3 , "dofPos=" << dofPos );
			Array<Point> tmp_points(1);
			tmp_points[0] = dofPos;
			refEval( maximalCellType , tmp_points ,
				     deriv,  result , 4);
		   break;
		   }
		}
	break;} // ------ end 2D -------
	case 3:{ //  ----- 3D ------------

		maximalCellType = BrickCell;

		Point dofPos( 0.0 , 0.0 , 0.0);
        // get the position of the DoF on the refCell
		SUNDANCE_OUT(this->verb() > 3 , ",maxCellDim=" << maxCellDim << ",facetDim="
				<< facetDim << ",facetIndex=" << facetIndex << ",nodeIndex=" << nodeIndex
				<< ",dofPos=" << dofPos);
		returnDoFPositionOnRef( maxCellDim,
		  facetDim, facetIndex, nodeIndex, dofPos);

		// shitf the position of the DOFs according to the position of the child in the parent cell
		switch (maxNrChild){
		case 27:{ // tri section on bricks
			dofPos = dofPos/3.0;
			SUNDANCE_OUT(this->verb() > 3 ,  "dofPos=" << dofPos );
			switch (indexInParent){
			  case 0:  {dofPos[0] += 0.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 1:  {dofPos[0] += 1.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 2:  {dofPos[0] += 2.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 3:  {dofPos[0] += 0.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 4:  {dofPos[0] += 1.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 5:  {dofPos[0] += 2.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 6:  {dofPos[0] += 0.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 7:  {dofPos[0] += 1.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 8:  {dofPos[0] += 2.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 0.0/3.0; break;}
			  case 9:  {dofPos[0] += 0.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 10: {dofPos[0] += 1.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 11: {dofPos[0] += 2.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 12: {dofPos[0] += 0.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  // 13 should never occur !!!!
			  case 13: {dofPos[0] += 1.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 14: {dofPos[0] += 2.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 15: {dofPos[0] += 0.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 16: {dofPos[0] += 1.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 17: {dofPos[0] += 2.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 1.0/3.0; break;}
			  case 18: {dofPos[0] += 0.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 19: {dofPos[0] += 1.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 20: {dofPos[0] += 2.0/3.0; dofPos[1] += 0.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 21: {dofPos[0] += 0.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 22: {dofPos[0] += 1.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 23: {dofPos[0] += 2.0/3.0; dofPos[1] += 1.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 24: {dofPos[0] += 0.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 25: {dofPos[0] += 1.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 2.0/3.0; break;}
			  case 26: {dofPos[0] += 2.0/3.0; dofPos[1] += 2.0/3.0; dofPos[2] += 2.0/3.0; break;}
			}
		break;}
		case 8:{ // bisection on bricks
			dofPos = dofPos/2.0;
			SUNDANCE_OUT(this->verb() > 3 ,  "dofPos=" << dofPos );
			switch (indexInParent){
			  case 0:  {dofPos[0] += 0.0/2.0; dofPos[1] += 0.0/2.0; dofPos[2] += 0.0/2.0; break;}
			  case 1:  {dofPos[0] += 1.0/2.0; dofPos[1] += 0.0/2.0; dofPos[2] += 0.0/2.0; break;}
			  case 2:  {dofPos[0] += 0.0/2.0; dofPos[1] += 1.0/2.0; dofPos[2] += 0.0/2.0; break;}
			  case 3:  {dofPos[0] += 1.0/2.0; dofPos[1] += 1.0/2.0; dofPos[2] += 0.0/2.0; break;}
			  case 4:  {dofPos[0] += 0.0/2.0; dofPos[1] += 0.0/2.0; dofPos[2] += 1.0/2.0; break;}
			  case 5:  {dofPos[0] += 1.0/2.0; dofPos[1] += 0.0/2.0; dofPos[2] += 1.0/2.0; break;}
			  case 6:  {dofPos[0] += 0.0/2.0; dofPos[1] += 1.0/2.0; dofPos[2] += 1.0/2.0; break;}
			  case 7:  {dofPos[0] += 1.0/2.0; dofPos[1] += 1.0/2.0; dofPos[2] += 1.0/2.0; break;}
			}
		break;}
		}

		SUNDANCE_OUT(this->verb() > 3 , "dofPos=" << dofPos );
		Array<Point> tmp_points(1);
		tmp_points[0] = dofPos;
		refEval( maximalCellType , tmp_points ,
		     deriv,  result , 4);

	break;} //end 3D
	}

    // evalute the results from the bi or trisection
	Array<double>& res = result[0][0];
	//SUNDANCE_OUT(this->verb() > 3 , "res=" << res );

	// get those data only once (since it is rather costly)
	if (doFInfromationCalculated_ == false)
	{
	  getDoFsOrdered( maximalCellType, res.size() , facetD_, facetI_, facetN_);
	  doFInfromationCalculated_ = true;
	}

	// here take the "result" array and get the nonzero elements
	index = 0;
	for ( int i=0 ; i < res.size() ; i++){
		// take only those DoFs which are non zero
		if ( fabs(res[i]) > 1e-5 ){
			SUNDANCE_OUT(this->verb() > 3 , "found DoF i=" << i );
			localDoFs.resize(index+1);
			coefs.resize(index+1);
			parentFacetDim.resize(index+1);
			parentFacetIndex.resize(index+1);
			parentFacetNode.resize(index+1);
			localDoFs[index] = i;
			coefs[index] = res[i];
			parentFacetDim[index] = facetD_[i];
			parentFacetIndex[index] = facetI_[i];
			parentFacetNode[index] = facetN_[i];
			index++;
		}
	}

	SUNDANCE_OUT(this->verb() > 3 , "localDoFs=" << localDoFs );
	SUNDANCE_OUT(this->verb() > 3 , "coefs=" << coefs );
	SUNDANCE_OUT(this->verb() > 3 , "parentFacetDim=" << parentFacetDim );
	SUNDANCE_OUT(this->verb() > 3 , "parentFacetIndex=" << parentFacetIndex );
	SUNDANCE_OUT(this->verb() > 3 , "parentFacetNode=" << parentFacetNode );

}

/* -------------- Get the position of one DoF ----------------------- */
void Lagrange::returnDoFPositionOnRef(
	const int maxCellDim,
	const int facetDim,
	const int facetIndex,
	const int nodeIndex,
	Point& pos) const {

	switch (facetDim){
	case 0:{  // Point DoFs
		switch (maxCellDim){
		case 2:{
			switch (facetIndex){
			case 0: {pos[0] = 0.0; pos[1] = 0.0; break;}
			case 1: {pos[0] = 1.0; pos[1] = 0.0; break;}
			case 2: {pos[0] = 0.0; pos[1] = 1.0; break;}
			case 3: {pos[0] = 1.0; pos[1] = 1.0; break;}
			}
			break;}
		case 3:{
			switch (facetIndex){
			case 0: {pos[0] = 0.0; pos[1] = 0.0; pos[2] = 0.0; break;}
			case 1: {pos[0] = 1.0; pos[1] = 0.0; pos[2] = 0.0; break;}
			case 2: {pos[0] = 0.0; pos[1] = 1.0; pos[2] = 0.0; break;}
			case 3: {pos[0] = 1.0; pos[1] = 1.0; pos[2] = 0.0; break;}
			case 4: {pos[0] = 0.0; pos[1] = 0.0; pos[2] = 1.0; break;}
			case 5: {pos[0] = 1.0; pos[1] = 0.0; pos[2] = 1.0; break;}
			case 6: {pos[0] = 0.0; pos[1] = 1.0; pos[2] = 1.0; break;}
			case 7: {pos[0] = 1.0; pos[1] = 1.0; pos[2] = 1.0; break;}
			}
			break;}
		}
		break;}
	case 1:{  // Edge DoFs
		switch (maxCellDim){
		case 2:{
			// here we calculate the position of the DoF on the edge
			double posOnEdge = (1.0/ (double)(order())) +(double)(nodeIndex) / (double)(order()-1);
			switch (facetIndex){
			case 0: {pos[0] = 0.0+posOnEdge; pos[1] = 0.0; break;}
			case 1: {pos[0] = 0.0; pos[1] = 0.0+posOnEdge; break;}
			case 2: {pos[0] = 1.0; pos[1] = 0.0+posOnEdge; break;}
			case 3: {pos[0] = 0.0+posOnEdge; pos[1] = 1.0; break;}
			}
			break;}
		case 3:{ // Edge in the brick test this
			double posOnEdge = (1.0/ (double)(order())) +(double)(nodeIndex) / (double)(order()-1);
			switch (facetIndex){
			case 0: {pos[0] = 0.0+posOnEdge; pos[1] = 0.0; pos[2] = 0.0; break;}
			case 1: {pos[0] = 0.0; pos[1] = 0.0+posOnEdge; pos[2] = 0.0; break;}
			case 2: {pos[0] = 0.0; pos[1] = 0.0; pos[2] = 0.0+posOnEdge; break;}
			case 3: {pos[0] = 1.0; pos[1] = 0.0+posOnEdge; pos[2] = 0.0; break;}
			case 4: {pos[0] = 1.0; pos[1] = 0.0; pos[2] = 0.0+posOnEdge; break;}
			case 5: {pos[0] = 0.0+posOnEdge; pos[1] = 1.0; pos[2] = 0.0; break;}
			case 6: {pos[0] = 0.0; pos[1] = 1.0; pos[2] = 0.0+posOnEdge; break;}
			case 7: {pos[0] = 1.0; pos[1] = 1.0; pos[2] = 0.0+posOnEdge; break;}
			case 8: {pos[0] = 0.0+posOnEdge; pos[1] = 0.0; pos[2] = 1.0; break;}
			case 9: {pos[0] = 0.0; pos[1] = 0.0+posOnEdge; pos[2] = 1.0; break;}
			case 10: {pos[0] = 1.0; pos[1] = 0.0+posOnEdge; pos[2] = 1.0; break;}
			case 11: {pos[0] = 0.0+posOnEdge; pos[1] = 1.0; pos[2] = 1.0; break;}
			}
			break;}
		}
		break;}
	case 2:{ //face DoFs only in 3D
		// DoFs position on the face
		int nrDoFoFace = (order()-1)*(order()-1);
		double posOnFace1 = (1.0/ (double)(order())) +(double)(nodeIndex % nrDoFoFace) / (double)(order()-1);
		double posOnFace2 = (1.0/ (double)(order())) +(double)(nodeIndex / nrDoFoFace) / (double)(order()-1);
		switch (facetIndex){
		case 0: {pos[0] = 0.0+posOnFace1; pos[1] = 0.0+posOnFace2; pos[2] = 0.0; break;}
		case 1: {pos[0] = 0.0+posOnFace1; pos[1] = 0.0; pos[2] = 0.0+posOnFace2; break;}
		case 2: {pos[0] = 0.0; pos[1] = 0.0+posOnFace1; pos[2] = 0.0+posOnFace2; break;}
		case 3: {pos[0] = 1.0; pos[1] = 0.0+posOnFace1; pos[2] = 0.0+posOnFace2; break;}
		case 4: {pos[0] = 0.0+posOnFace1; pos[1] = 1.0; pos[2] = 0.0+posOnFace2; break;}
		case 5: {pos[0] = 0.0+posOnFace1; pos[1] = 0.0+posOnFace2; pos[2] = 1.0; break;}
		}
		break;}
	}

}

void Lagrange::getDoFsOrdered(
		const CellType maxCellDim,
		int nrDoF,
		Array<int>& facetD,
		Array<int>& facetI,
		Array<int>& facetN)
{
	// this could be done once at initialization
	// ----------------------------------------------------------------
	Array<Array<Array<int> > > dofs_struct;
	getReferenceDOFs(maxCellDim , maxCellDim , dofs_struct);
	bool getNext;
	facetD.resize(nrDoF); facetI.resize(nrDoF); facetN.resize(nrDoF);
	int dimo=0 , indexo=0 , nodo=0;
	for ( int i=0 ; i < nrDoF ; i++){
		facetD[i] = dimo;
		facetI[i] = indexo;
		facetN[i] = nodo;
		// found the next DoF and change the variables
		// nodo,indexo,dimo
		getNext = true;
		while ( (getNext) && (dimo < dofs_struct.size() )){
			//SUNDANCE_OUT(this->verb() > 3, "dimo=" << dimo << " ,indexo:" << indexo << ",nodo:" << nodo <<
			//		",maxCellDim:" << maxCellDim << " ,nrDoF:" << nrDoF);
			if (dofs_struct[dimo].size() > (indexo+1))
			{
			    if (dofs_struct[dimo][indexo].size() > (nodo+1)){
				  nodo++; getNext = false;
			   }else{
				  nodo = 0;
				  if (dofs_struct[dimo].size() > (indexo+1)){
					  indexo++; getNext = false;
				  }else{
					  indexo = 0;
					  if (dofs_struct.size() > (dimo+1)){
						  dimo++; getNext = false;
					  }else{
						  // actually we should not arrived here
						  getNext = false;
					  }
				  }
			   }
		    }
		    else
		    {
               dimo++; indexo = 0; nodo = 0; getNext = false;
			}
			//SUNDANCE_OUT(this->verb() > 3, "dofs_struct.size()=" << dofs_struct.size() << " ,dofs_struct[dimo].size():"
			//		<< dofs_struct[dimo].size() << ",dofs_struct[dimo][indexo].size():" << dofs_struct[dimo][indexo].size());
			if (dofs_struct.size() <= dimo ) {dimo++; getNext = true; continue;}
			if (dofs_struct[dimo].size() <= indexo ) {dimo++; getNext = true; continue;}
			if (dofs_struct[dimo][indexo].size() <= nodo ) {dimo++; getNext = true; continue;}
		}
	}
	SUNDANCE_OUT(this->verb() > 3, "facetD=" << facetD );
	SUNDANCE_OUT(this->verb() > 3, "facetI=" << facetI );
	SUNDANCE_OUT(this->verb() > 3, "facetN=" << facetN );
    // ------------------------------------------------
}
