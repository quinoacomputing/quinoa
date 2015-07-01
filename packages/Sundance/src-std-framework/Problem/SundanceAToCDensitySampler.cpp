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

#include "SundanceAToCDensitySampler.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceGeomUtils.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include <queue>

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#endif

using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


static Time& densitySamplingTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("density sampling"); 
  return *rtn;
}

AToCDensitySampler::AToCDensitySampler(const AToCPointLocator& locator,
                                       const VectorType<double>& vecType)
  : discSpace_(locator.mesh(), new Lagrange(0), locator.subdomain(), vecType),
    dim_(locator.mesh().spatialDim()),
    mesh_(locator.mesh()),
    elemToVecIndexMap_(),
    elemWeights_(new DiscreteFunction(discSpace_, 0.0)),
    elemWeightVec_(),
    locator_(locator),
    isAxisymmetric_(false),
    origin_(),
    axis_()
{
  init();
}


AToCDensitySampler::AToCDensitySampler(const AToCPointLocator& locator,
                                       const std::vector<double>& origin,
                                       const std::vector<double>& rotationalAxis,
                                       const VectorType<double>& vecType)
  : discSpace_(locator.mesh(), new Lagrange(0), locator.subdomain(), vecType),
    dim_(locator.mesh().spatialDim()),
    mesh_(locator.mesh()),
    elemToVecIndexMap_(),
    elemWeights_(new DiscreteFunction(discSpace_, 0.0)),
    elemWeightVec_(),
    locator_(locator),
    isAxisymmetric_(true),
    origin_(vec2point(origin)),
    axis_(normPoint(vec2point(rotationalAxis)))
{
  init();
}



void AToCDensitySampler::init()
{
  const CellFilter& domain = discSpace_.cellFilters(0);

  elemWeightVec_ = DiscreteFunction::discFunc(elemWeights_)->getVector();

  elemToVecIndexMap_ = rcp(new Array<int>(mesh_.numCells(dim_), -1));

  Array<int>& a = *elemToVecIndexMap_;

  CellSet cells = domain.getCells(mesh_);

  Array<int> cellLID;
  cellLID.reserve(mesh_.numCells(dim_));

  for (CellIterator i=cells.begin(); i!=cells.end(); i++)
    {
      cellLID.append(*i);
    }

  const RCP<DOFMapBase>& dofMap = discSpace_.map();

  Set<int> funcs = makeSet(0);
  Array<Array<int> > dofs;
  Array<int> nNodes;
  dofMap->getDOFsForCellBatch(dim_, cellLID, funcs, dofs, nNodes,0);
  
  const Array<int>& dofs0 = dofs[0];
  for (int c=0; c<cellLID.size(); c++)
    {
      int vecIndex = dofs0[c];
      int lid = cellLID[c];
      a[lid] = vecIndex;
      double vol = volume(mesh_, dim_, lid);
      if (isAxisymmetric_)
        {
          Point xCell = mesh_.centroid(dim_, lid) - origin_;
          double dPerp = ::sqrt(xCell*xCell - (xCell*axis_)*(xCell*axis_));
          vol = vol * dPerp;
        }
      elemWeightVec_[vecIndex] = vol;
    }
}

Point AToCDensitySampler::vec2point(const std::vector<double>& x) const
{
  if (x.size()==1) return Point(x[0]);
  else if (x.size()==2U) return Point(x[0], x[1]);
  else if (x.size()==3U) return Point(x[0], x[1], x[2]);
  TEUCHOS_TEST_FOR_EXCEPT(x.size() < 1 || x.size() > 3U);
  return Point();
}

Point AToCDensitySampler::normPoint(const Point& x) const
{
  return x/sqrt(x*x);
}


Expr AToCDensitySampler::sample(const std::vector<double>& positions,
                                const double& particleWeight) const 
{
  TimeMonitor timer(densitySamplingTimer());

  TEUCHOS_TEST_FOR_EXCEPTION(positions.size() % dim_ != 0, std::runtime_error,
                     "vector of coordinates should by an integer multiple "
                     "of the spatial dimension");

  Expr rtn = new DiscreteFunction(discSpace_, 0.0);
  Vector<double> density = DiscreteFunction::discFunc(rtn)->getVector();

  int nPts = positions.size() / dim_;

  for (int i=0; i<nPts; i++)
    {
      const double* x = &(positions[dim_*i]);

      int guess = locator_.guessCell(x);

      TEUCHOS_TEST_FOR_EXCEPTION(guess < 0, std::runtime_error, "particle #" << i << " position="
                         << AToCPointLocator::makePoint(dim_, x) 
                         << " is out of search grid");

      int cellLID = locator_.findEnclosingCell(guess, x);

      int vecIndex = (*elemToVecIndexMap_)[cellLID];
      double vol = elemWeightVec_[vecIndex];
      density[vecIndex] += particleWeight/vol;
    }

  return rtn;
}


Expr AToCDensitySampler::resetCounts() const 
{
  Expr rtn = new DiscreteFunction(discSpace_, 0.0);

  return rtn;
}

void AToCDensitySampler::addToCounts(const std::vector<double>& positions,
                                     const double& particleWeight,
                                     Expr density) const 
{
  TimeMonitor timer(densitySamplingTimer());

  TEUCHOS_TEST_FOR_EXCEPTION(positions.size() % dim_ != 0, std::runtime_error,
                     "vector of coordinates should by an integer multiple "
                     "of the spatial dimension");

  Vector<double> vec = DiscreteFunction::discFunc(density)->getVector();

  int nPts = positions.size() / dim_;

  for (int i=0; i<nPts; i++)
    {
      const double* x = &(positions[dim_*i]);

      int guess = locator_.guessCell(x);

      TEUCHOS_TEST_FOR_EXCEPTION(guess < 0, std::runtime_error, "particle #" << i << " position="
                         << AToCPointLocator::makePoint(dim_, x) 
                         << " is out of search grid");

      int cellLID = locator_.findEnclosingCell(guess, x);

      int vecIndex = (*elemToVecIndexMap_)[cellLID];
      double vol = elemWeightVec_[vecIndex];
      vec[vecIndex] += particleWeight/vol;
    }
}


