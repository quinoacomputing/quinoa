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

#include "SundanceBinaryCellFilter.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOrderedTuple.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

BinaryCellFilter::BinaryCellFilter(const CellFilter& left,
                                   const CellFilter& right,
                                   const CellFilterOpType& op)
  : CellFilterBase(), op_(op), left_(left), right_(right)
{
  std::string str;
  switch(op)
  {
    case Union:
      str = "Union(";
      break;
    case Intersection:
      str = "Intersection(";
      break;
    default:
      str = "SetDifference(";
    }

  setName(str + left.toString() + ", " + right.toString() + ")");
}

int BinaryCellFilter::dimension(const Mesh& mesh) const
{
  int d1 = left_.dimension(mesh);
  int d2 = right_.dimension(mesh);

  TEUCHOS_TEST_FOR_EXCEPTION(d1 != d2, std::runtime_error,
                     "BinaryCellFilter::dimension() mismatched dimensions. "
                     "Left filter has dimension d1=" << d1 << " but "
                     "right filter has dimension d2=" << d2);

  return d1;
}

CellSet BinaryCellFilter::internalGetCells(const Mesh& mesh) const
{
  Tabs tab;

  SUNDANCE_OUT(this->verb() > 2,
               "cell filter " << toXML().toString() << " is getting its cells for mesh " 
               << mesh.id());

  CellSet L = left_.getCells(mesh);
  CellSet R = right_.getCells(mesh);


  SUNDANCE_OUT(this->verb() > 2,
               "cell filter " << toXML().toString() << " is performing its operation");
  
  switch(op_)
    {
    case Union:
      return L.setUnion(R);
    case Intersection:
      return L.setIntersection(R);
    case Difference:
      return L.setDifference(R);
    }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "unknown cell filter op type" << op_ 
                     << " in BinaryCellFilter::internalGetCells()");
  return L; // -Wall
}

string BinaryCellFilter::opName() const 
{
  switch(op_)
    {
    case Union:
      return "UnionCellFilter";
    case Intersection:
      return "IntersectionCellFilter";
    case Difference:
      return "DifferenceCellFilter";
    }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "unknown cell filter op type" << op_ 
                     << " in BinaryCellFilter::opName()");
  return "UnionCellFilter";//-Wall
}

XMLObject BinaryCellFilter::toXML() const 
{
  XMLObject rtn(opName());
  rtn.addChild(left_.toXML());
  rtn.addChild(right_.toXML());
  return rtn;
}


#ifdef OLD_CELL_FILTER

bool BinaryCellFilter::lessThan(const CellFilterStub* other) const
{
  const BinaryCellFilter* B 
    = dynamic_cast<const BinaryCellFilter*>(other);

  TEUCHOS_TEST_FOR_EXCEPTION(B==0,
                     std::logic_error,
                     "argument " << other->toXML() 
                     << " to BinaryCellFilter::lessThan() should be "
                     "a BinaryCellFilter pointer.");

  return OrderedTriple<CellFilterOpType, CellFilter, CellFilter>(op_, left_, right_) 
    < OrderedTriple<CellFilterOpType, CellFilter, CellFilter>(B->op_, B->left_, B->right_) ;
}

#endif
