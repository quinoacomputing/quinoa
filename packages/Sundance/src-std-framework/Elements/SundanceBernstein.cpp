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

#include "SundanceBernstein.hpp"
#include "SundanceADReal.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundancePoint.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Teuchos;


Bernstein::Bernstein(int order)
  : order_(order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(order < 0, std::runtime_error,
    "invalid polynomial order=" << order
    << " in Bernstein ctor");
}

bool Bernstein::supportsCellTypePair(
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

void Bernstein::print(std::ostream& os) const 
{
  os << "Bernstein(" << order_ << ")";
}

int Bernstein::nReferenceDOFsWithoutFacets(
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
      if (order_ >= 3) return (order_-1)*(order_-2)/2;
      break;
    case QuadCell:
      if (order_==1) return 0;
      if (order_==2) return 1;
      TEUCHOS_TEST_FOR_EXCEPTION(order_>2, std::runtime_error, 
        "Bernstein order > 2 not implemented "
        "for quad cells");
    case TetCell:
      if (order_<=2) return 0;
      TEUCHOS_TEST_FOR_EXCEPTION(order_>2, std::runtime_error, 
        "Bernstein order > 2 not implemented "
        "for tet cells");
    case BrickCell:
      if (order_<=1) return 0;
      if (order_==2) return 1;
      TEUCHOS_TEST_FOR_EXCEPTION(order_>2, std::runtime_error, 
        "Bernstein order > 2 not implemented "
        "for brick cells");
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
        << cellType << " not implemented in Bernstein basis");
      return -1; // -Wall
  }
  return -1; // -Wall
}

void Bernstein::getReferenceDOFs(
  const CellType& maximalCellType,
  const CellType& cellType,
  Array<Array<Array<int> > >& dofs) const 
{
  const int N=(order()+1)*(order()+2)/2;
  typedef Array<int> Aint;
  switch(cellType)
  {
    case PointCell:
      dofs.resize(1);
      dofs[0] = tuple<Aint>(tuple(0));
      return;
    case LineCell:
      dofs.resize(2);
      dofs[0] = tuple<Aint>(tuple(0), tuple(order()));
      dofs[1] = tuple<Aint>(makeRange(1, order()-1));
      return;
    case TriangleCell:
    {
      dofs.resize(3);
      dofs[0] = tuple<Aint>(tuple(0),tuple(N-order()-1),tuple(N-1));
      if (order()>1)
      {
        // order() - 1 dof per edge
        Aint dof0(order()-1);
        Aint dof1; // will fill in with range
        Aint dof2(order()-1);

        // first edge (from vertex 0 to 1
        int inc = 2;
        dof0[0] = 1;
        for (int i=1;i<order()-1;i++) 
	      {
          dof0[i] = dof0[i-1]+inc;
          inc++;
	      }

        // second edge runs from vertex 1 to 2
        dof1 = makeRange(N-order(),N-2);

        // third edge runs from vertex 2 to vertex 0
        inc = 3;
        dof2[order()-2] = 2;
        for (int i=order()-3;i>=0;i--)
	      {
          dof2[i] = dof2[i+1]+inc;
          inc++;
	      }
        /** permuted from (0,1,2) to (1,2,0) for UFC ordering */
        dofs[1] = tuple<Aint>(dof1,dof2,dof0);
      }
      else
      {
        dofs[1] = tuple(Array<int>());
      }
      if (order()>2)
      {
        Array<int> internaldofs;
        int bfcur = 0;
        //int internalbfcur=0; // Commented out unused variable -- KL
        for (int alpha1=order();alpha1>=0;alpha1--)
	      {
          for (int alpha2=order()-alpha1;alpha2>=0;alpha2--)
          {
            int alpha3 = order()-alpha1-alpha2;
            if (alpha1*alpha2*alpha3>0) {
              internaldofs.append(bfcur);
            }
            bfcur++;
          }
	      }
        dofs[2] = tuple(internaldofs);
      }
      else
      {
        dofs[2] = tuple(Array<int>());
      }
      return;
    }
    case TetCell:
    {
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cell type "
        << cellType << " not implemented in Bernstein basis");
  }
}



Array<int> Bernstein::makeRange(int low, int high)
{
  if (high < low) return Array<int>();

  Array<int> rtn(high-low+1);
  for (int i=0; i<rtn.length(); i++) rtn[i] = low+i;
  return rtn;
}

void Bernstein::refEval(
  const CellType& cellType,
  const Array<Point>& pts,
  const SpatialDerivSpecifier& sds,
  Array<Array<Array<double> > >& result,
  int verbosity) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(sds.isPartial() || sds.isIdentity()), 
    std::runtime_error,
    "cannot evaluate spatial derivative " << sds << " on Bernstein basis");
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
    case TetCell:
      for (int i=0; i<pts.length(); i++)
      {
        evalOnTet(pts[i], deriv, result[0][i]);
      }
      return;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Bernstein::refEval() unimplemented for cell type "
        << cellType);

  }
}

/* ---------- evaluation on different cell types -------------- */

void Bernstein::evalOnLine(const Point& pt, 
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
    double binom_cur = 1.0;
    for (int i=0; i<=order_; i++)
    {
      tmp[i] = one;

      for (int j=0;j<order_-i;j++) 
	    {
	      tmp[i] *= (1-x);
	    }
      for (int j=0;j<i;j++) 
	    {
	      tmp[i] *= x;
	    }
      tmp[i] *= binom_cur;
      binom_cur *= double(order()-i) / double(i+1);
    }
  }
  
  for (int i=0; i<tmp.length(); i++)
  {
    if (deriv.order()==0) result[i] = tmp[i].value();
    else result[i] = tmp[i].gradient()[0];
  }
}

void Bernstein::evalOnTriangle(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const

{
  ADReal x = ADReal(pt[0], 0, 2);
  ADReal y = ADReal(pt[1], 1, 2);
  ADReal one(1.0, 2);
  
  Array<ADReal> tmp;

  SUNDANCE_OUT(this->verb() > 3, "x=" << x.value() << " y="
    << y.value());
 
  if(order_==0) {
    result.resize(1);
    tmp.resize(1);
    tmp[0] = one;
  }
  else {
    int N = (order()+1)*(order()+2)/2;
    result.resize(N);
    tmp.resize(N);
    // these are the barycentric coordinates
    ADReal b1 = 1.0 - x - y;
    ADReal b2 = x;
    ADReal b3 = y;

    // will hold \binom{n}{\alpha_1}
    int bfcur = 0;

    for (int alpha1=order();alpha1>=0;alpha1--) 
    {
      for (int alpha2 = order()-alpha1;alpha2>=0;alpha2--)
      {
        int alpha3 = order() - alpha1 - alpha2;
        tmp[bfcur] = one;
        for (int i=0;i<alpha1;i++)
	      {
          tmp[bfcur] *= b1;
	      }
        for (int i=0;i<alpha2;i++) 
	      {
          tmp[bfcur] *= b2;
	      }
        for (int i=0;i<alpha3;i++)
	      {
          tmp[bfcur] *= b3;
	      }
        for (int i=2;i<=order();i++)
	      {
          tmp[bfcur] *= (double) i;
	      }
        for (int i=2;i<=alpha1;i++) 
	      {
          tmp[bfcur] /= (double) i;
	      }
        for (int i=2;i<=alpha2;i++) 
	      {
          tmp[bfcur] /= (double) i;
	      }
        for (int i=2;i<=alpha3;i++) 
	      {
          tmp[bfcur] /= (double) i;
	      }
        bfcur++;
      }
    }
  }

  for (int i=0; i<tmp.length(); i++)
  {
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
}

void Bernstein::evalOnTet(const Point& pt, 
  const MultiIndex& deriv,
  Array<double>& result) const
{
  ADReal x = ADReal(pt[0], 0, 3);
  ADReal y = ADReal(pt[1], 1, 3);
  ADReal z = ADReal(pt[2], 2, 3);
  ADReal one(1.0, 3);
  
  Array<ADReal> tmp(result.length());

  if(order_==0)
  {
    tmp.resize(1);
    result.resize(1);
    tmp[0] = one;
  }
  else
  {
  }

  for (int i=0; i<tmp.length(); i++)
  {
    if (deriv.order()==0) result[i] = tmp[i].value();
    else 
      result[i] = tmp[i].gradient()[deriv.firstOrderDirection()];
  }
}




