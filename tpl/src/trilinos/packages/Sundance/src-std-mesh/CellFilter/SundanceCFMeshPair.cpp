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

#include "SundanceCFMeshPair.hpp"
#include "PlayaTabs.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceBinaryCellFilter.hpp"
#include "SundanceSubsetCellFilter.hpp"
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceNullCellFilter.hpp"

#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


static Time& csPartitionTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("cell set partitioning"); 
  return *rtn;
}


CFMeshPair::CFMeshPair()
  : filter_(),
    mesh_(),
    cellSet_(),
    funcs_()
{}

CFMeshPair::CFMeshPair(const CellFilter& cf,
                       const Mesh& mesh,
                       const Set<int>& funcs)
  : filter_(cf),
    mesh_(mesh),
    cellSet_(),
    funcs_(funcs)
{
  if (filter_.ptr().get() != 0) cellSet_ = filter_.getCells(mesh);
}

bool CFMeshPair::operator<(const CFMeshPair& other) const
{
  if (isEmpty()) return true;
  if (other.isEmpty()) return false;

  TEUCHOS_TEST_FOR_EXCEPTION(mesh_.id() != other.mesh_.id(),
                     std::runtime_error,
                     "mismatched meshes!");

  return cellSet_ < other.cellSet_;
}

bool CFMeshPair::isEmpty() const
{
  return filter_.ptr().get() == 0 
    || cellSet_.ptr().get() == 0
    || cellSet_.begin()==cellSet_.end();
}

CFMeshPair CFMeshPair::setMinus(const CFMeshPair& other) const
{

  if (isEmpty()) return CFMeshPair();
  if (other.isEmpty()) return *this;

  TEUCHOS_TEST_FOR_EXCEPTION(mesh().id() != other.mesh().id(),
                     std::runtime_error,
                     "mismatched meshes!");

  CellFilter diff = filter() - other.filter();
  return CFMeshPair(diff, mesh(), funcs_);
}

CFMeshPair CFMeshPair::intersection(const CFMeshPair& other) const
{
  if (isEmpty() || other.isEmpty()) return CFMeshPair();

  TEUCHOS_TEST_FOR_EXCEPTION(mesh().id() != other.mesh().id(),
                     std::runtime_error,
                     "mismatched meshes!");

  CellFilter inter = filter().intersection(other.filter());

  return CFMeshPair(inter, mesh(), funcs_.setUnion(other.funcs_));
}


namespace Sundance
{
  Array<CFMeshPair> resolvePair(const CFMeshPair& a, 
                                const CFMeshPair& b)
  {
    Tabs tab0;

    CFMeshPair inter = a.intersection(b);
    CFMeshPair aMinusB = a.setMinus(b);
    CFMeshPair bMinusA = b.setMinus(a);

    return tuple(bMinusA, inter, aMinusB);
  }

  
  Array<CFMeshPair> resolveSets(const Array<CFMeshPair>& s)
  {

    if (s.size() == 1) return s;

    Array<CFMeshPair> T = tuple(s[0]);
    
    for (int i=0; i<s.size(); i++)
      {
        CFMeshPair A = s[i];
        Array<CFMeshPair> TN;
        for (int j=0; j<T.size(); j++)
          {
            Array<CFMeshPair> p = resolvePair(A, T[j]);
            if (!p[0].isEmpty())
              {
                TN.append(p[0]);
              }
            if (!p[1].isEmpty())
              {
                TN.append(p[1]);
              }
            A = p[2];
          }
        if (!A.isEmpty())
          {
            TN.append(A);
          }
        T = TN;
      }
    return T;
  }


  Array<CFMeshPair>
  findDisjointFilters(const Array<CellFilter>& filters,
                      const Array<Set<int> >& funcs,
                      const Mesh& mesh)
  {
    TimeMonitor timer(csPartitionTimer() );
    Array<CFMeshPair> cf(filters.size());

    TEUCHOS_TEST_FOR_EXCEPT(filters.size() != funcs.size());
    for (int i=0; i<filters.size(); i++)
      {
        cf[i] = CFMeshPair(filters[i], mesh, funcs[i]);
      }
    return resolveSets(cf);
  }

}
