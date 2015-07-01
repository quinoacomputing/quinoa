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

#include "SundanceCellSet.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"
#include <algorithm>
#include <iterator>

using namespace Sundance;
using namespace Teuchos;
using Playa::Handle;
using Playa::Handleable;


CellSet::CellSet(const Mesh& mesh, int cellDim,
  const CellType& cellType,
  const Set<int>& cellLIDs)
  : Handle<CellSetBase>(rcp(new ExplicitCellSet(mesh, cellDim, cellType, cellLIDs)))
{}

CellSet CellSet::setUnion(const CellSet& other) const
{
  if (isNull()) return other;
  if (other.isNull()) return *this;

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("union", other);
  
  Set<int>& cells = rtn->cells();

  std::set_union(this->begin(), this->end(), other.begin(), other.end(), 
    std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setIntersection(const CellSet& other) const
{
  if (isNull()) return *this;
  if (other.isNull()) return other;

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("intersection", other);
  
  Set<int>& cells = rtn->cells();

  std::set_intersection(this->begin(), this->end(), other.begin(), other.end(), 
    std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}

CellSet CellSet::setDifference(const CellSet& other) const
{
  if (isNull()) return *this;
  if (other.isNull()) return *this;

  ExplicitCellSet* rtn = new ExplicitCellSet(mesh(), dimension(), cellType());

  checkCompatibility("difference", other);
  
  Set<int>& cells = rtn->cells();

  std::set_difference(this->begin(), this->end(), other.begin(), other.end(), 
    std::insert_iterator<Set<int> >(cells, cells.begin()));
  
  return rtn;
}


void CellSet::checkCompatibility(const std::string& op, const CellSet& other) const 
{

  TEUCHOS_TEST_FOR_EXCEPTION(meshID() != other.meshID(), std::runtime_error,
    "CellSet::checkCompatibility(): "
    "incompatible mesh ID numbers in " << op
    << ". LHS=" << meshID() << " RHS=" << other.meshID());

  TEUCHOS_TEST_FOR_EXCEPTION(dimension() != other.dimension(), std::runtime_error,
    "CellSet::checkCompatibility() incompatible dimensions in " << op
    << "LHS has "
    "dimension=" << dimension() << " but RHS has dimension="
    << other.dimension());
  
  TEUCHOS_TEST_FOR_EXCEPTION(cellType() != other.cellType(), std::runtime_error,
    "CellSet::checkCompatibility() incompatible cell types. "
    " in " << op << " LHS has "
    "cellType=" << cellType() << " but RHS has cellType="
    << other.cellType());


  SUNDANCE_OUT(this->verb() > 2,
    "Set operation: " << op << std::endl
    << "LHS cells: " << *this << std::endl
    << "RHS cells: " << other);
               
}


bool CellSet::areFacetsOf(const CellSet& other) const
{
  Array<int> cofacetLIDs;
  int myDim = dimension();
  int cofacetDim = other.dimension();
  CellType cofacetType = other.cellType();
  if (myDim >= cofacetDim) return false;

  for (CellIterator i=begin(); i!=end(); i++)
  {
    int cellLID = *i;
    mesh().getCofacets(myDim, cellLID, cofacetDim, cofacetLIDs);
    Set<int> cofacetSet;
    for (int c=0; c<cofacetLIDs.size(); c++)
    {
      int cf = cofacetLIDs[c];
      cofacetSet.put(cf);
    }
    CellSet myCofacetSet(mesh(), cofacetDim, cofacetType, cofacetSet); 
    CellSet intersection = myCofacetSet.setIntersection(other);      
    /* if the intersection is empty, then we have found a cell
     * that is not a facet of one of the other cells */
    if (intersection.begin()==intersection.end()) return false;
  }
  return true;
}

bool CellSet::operator<(const CellSet& other) const
{
  Tabs tab;
  bool rtn = ptr()->lessThan(other.ptr().get());
  return rtn;
}


int CellSet::numCells() const 
{
  int count = 0;
  for (CellIterator i=begin(); i!=end(); i++)
  {
    count ++;
  }
  return count;
}

