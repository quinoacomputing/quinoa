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

#include "SundanceCellFilter.hpp"
#include "SundanceCellFilterBase.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceNullCellFilter.hpp"
#include "PlayaTabs.hpp"
#include "SundanceSubsetManager.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

bool CellFilter::isNullCellFilter() const 
{
  return dynamic_cast<const NullCellFilterStub*>(ptr().get()) != 0;
}

bool CellFilter::isNull() const
{
  return ptr().get() == 0 || isNullCellFilter();
}

void CellFilter::setName(const std::string& name)
{
  nonConstCfbPtr()->setName(name);
}

CellSet CellFilter::getCells(const Mesh& mesh) const
{
  if (isNull() || isNullCellFilter())
  {
    return new ExplicitCellSet(mesh, -1, 
      NullCell);
  }
  return cfbPtr()->getCells(mesh);
}



int CellFilter::dimension(const Mesh& mesh) const
{
  if (isNullCellFilter())
  {
    return -1;
  }
  return cfbPtr()->dimension(mesh);
}



CellFilter CellFilter::operator+(const CellFilter& other) const 
{
  if (isNull())
  {
    return other;
  }
  else if (other.isNull())
  {
    return *this;
  }
  else
  {
    CellFilter rtn 
      = new BinaryCellFilter(*this, other, BinaryCellFilter::Union);
    rtn.registerSubset(*this);
    rtn.registerSubset(other);
    return rtn;
  }
}



CellFilter CellFilter::operator-(const CellFilter& other) const 
{
  if (other.isNull())
  {
    return *this;
  }
  else if (isKnownDisjointWith(other) || other.isKnownDisjointWith(*this))
  {
    return *this;
  }
  else if (isKnownSubsetOf(other))
  {
    CellFilter rtn;
    return rtn;
  }
  else if (*this == other)
  {
    CellFilter rtn;
    return rtn;
  }
  else
  {
    CellFilter rtn 
      = new BinaryCellFilter(*this, other, BinaryCellFilter::Difference);
    rtn.registerDisjoint(other);
    this->registerSubset(rtn);
    return rtn;
  }
}



CellFilter CellFilter::intersection(const CellFilter& other) const 
{
  if (isNull() || other.isNull())
  {
    CellFilter rtn;
    return rtn;
  }
  else if (isKnownDisjointWith(other) || other.isKnownDisjointWith(*this))
  {
    CellFilter rtn;
    return rtn;
  }
  else if (isKnownSubsetOf(other))
  {
    return *this;
  }
  else if (other.isKnownSubsetOf(*this))
  {
    return other;
  }
  else if (*this==other)
  {
    return *this;
  }
  else
  {
    CellFilter rtn 
      = new BinaryCellFilter(*this, other, BinaryCellFilter::Intersection);
    other.registerSubset(rtn);
    this->registerSubset(rtn);
    
    return rtn;
  }
}



CellFilter CellFilter::labeledSubset(int label) const
{
  return labeledSubset(tuple(label));
}

CellFilter CellFilter::labeledSubset(const Array<int>& labels) const
{
  Set<int> labelSet = makeSet(labels);
  CellPredicate pred = new LabelCellPredicate(labelSet);
  CellFilter rtn = new SubsetCellFilter(*this, pred);
  this->registerLabeledSubset(labelSet, rtn);
  this->registerSubset(rtn);

  return rtn;
}

CellFilter CellFilter::coordSubset(int dir, const double& coordVal) const
{
  CellPredicate pred = new CoordinateValueCellPredicate(dir, coordVal);
  CellFilter rtn = new SubsetCellFilter(*this, pred);
  this->registerSubset(rtn);

  return rtn;
}

CellFilter CellFilter::subset(const CellPredicate& pred) const
{
  CellFilter rtn = new SubsetCellFilter(*this, pred);
  this->registerSubset(rtn);
  return rtn;
}


CellFilter CellFilter::subset(const RCP<CellPredicateFunctorBase>& test) const
{
  CellFilter rtn = new SubsetCellFilter(*this, CellPredicate(test));
  this->registerSubset(rtn);
  return rtn;
}

bool CellFilter::isKnownSubsetOf(const CellFilter& other) const
{
  if (other.knownSubsets().contains(*this)) return true;
  return false;
}

bool CellFilter::isKnownDisjointWith(const CellFilter& other) const
{
  if (other.knownDisjoints().contains(*this)) return true;
  if (this->knownDisjoints().contains(other)) return true;

  return false;
}

bool CellFilter::isSubsetOf(const CellFilter& other,
  const Mesh& mesh) const
{
  if (isKnownSubsetOf(other)) 
  {
    return true;
  }
  else
  {
    CellSet myCells = getCells(mesh);
    CellSet yourCells = other.getCells(mesh);
    CellSet inter = myCells.setIntersection(yourCells);
    if (inter.begin() == inter.end()) return false;
    CellSet diff = myCells.setDifference(inter);
    return (diff.begin() == diff.end());
  }
}



bool CellFilter::operator==(const CellFilter& other) const
{
  if (*this < other) return false;
  if (other < *this) return false;
  return true;
}

bool CellFilter::operator!=(const CellFilter& other) const
{
  return !( *this == other );
}


const Set<CellFilter>& CellFilter::knownSubsets() const
{
  return SubsetManager::getSubsets(*this);
}
const Set<CellFilter>& CellFilter::knownDisjoints() const
{
  return SubsetManager::getDisjoints(*this);
}

void CellFilter::registerSubset(const CellFilter& sub) const
{
  SubsetManager::registerSubset(*this, sub);
  
  for (Set<CellFilter>::const_iterator 
         i=sub.knownSubsets().begin(); i!=sub.knownSubsets().end(); i++)
  {
    SubsetManager::registerSubset(*this, *i);
  }
}


void CellFilter::registerDisjoint(const CellFilter& sub) const
{
  SubsetManager::registerDisjoint(*this, sub);
  
  for (Set<CellFilter>::const_iterator 
         i=sub.knownDisjoints().begin(); i!=sub.knownDisjoints().end(); i++)
  {
    SubsetManager::registerDisjoint(*this, *i);
  }
}

void CellFilter::registerLabeledSubset(const Set<int>& label, 
  const CellFilter& sub) const
{
  SubsetManager::registerLabeledSubset(*this, label, sub);
  
  const Map<Set<int>, CellFilter>& subsub = SubsetManager::getLabeledSubsets(sub);

  for (Map<Set<int>, CellFilter>::const_iterator 
         iter=subsub.begin(); iter!=subsub.end(); iter++)
  {
    if (iter->first == label) continue;
    SubsetManager::registerDisjoint(sub, iter->second);
    SubsetManager::registerDisjoint(iter->second, sub);
  }
}


XMLObject CellFilter::toXML() const 
{
  return ptr()->toXML();
}

string CellFilter::toString() const 
{
  return cfbPtr()->toString();
}

const CellFilterBase* CellFilter::cfbPtr() const
{
  const CellFilterBase* rtn = dynamic_cast<const CellFilterBase*>(ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(rtn==0, std::logic_error, "CellFilter::cfbPtr() cast failed");
  return rtn;
}

CellFilterBase* CellFilter::nonConstCfbPtr()
{
  CellFilterBase* rtn = dynamic_cast<CellFilterBase*>(ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(rtn==0, std::logic_error, "CellFilter::nonConstCfbPtr() cast failed");
  return rtn;
}
