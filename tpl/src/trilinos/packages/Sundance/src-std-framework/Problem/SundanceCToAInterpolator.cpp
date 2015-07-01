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

#include "SundanceCToAInterpolator.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


static Time& particleInterpolationTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("particle interpolation"); 
  return *rtn;
}

CToAInterpolator::CToAInterpolator(const AToCPointLocator& locator,
  const Expr& field)
  : dim_(locator.mesh().spatialDim()),
    nFacets_(dim_+1),
    rangeDim_(-1),
    elemToVecValuesMap_(),
    locator_(locator)
{
  updateField(field);
}

void CToAInterpolator::updateField(const Expr& field)
{
  int newRangeDim = field.size();
  if (newRangeDim != rangeDim_)
  {
    rangeDim_ = newRangeDim;
    elemToVecValuesMap_ 
      = rcp(new Array<double>(rangeDim_ * locator_.mesh().numCells(dim_) * nFacets_));
  }

  int nCells = locator_.mesh().numCells(dim_);

  const DiscreteFunction* df = DiscreteFunction::discFunc(field);
  TEUCHOS_TEST_FOR_EXCEPTION(df == 0, std::runtime_error,
    "discrete function expected in "
    "CToAInterpolator::updateField()");
  
  Vector<double> vec = df->getVector();      

  Array<int> cellLID(nCells);
  for (int i=0; i<nCells; i++)
  {
    cellLID[i] = i;
  }
  Array<Array<int> > dofs;
  Array<int> nNodes;
  Set<int> funcs;
  for (int i=0; i<rangeDim_; i++) funcs.put(i);

  
  df->map()->getDOFsForCellBatch(dim_, cellLID, funcs, dofs, nNodes,0);

  const Array<int>& dofs0 = dofs[0];

  for (int c=0; c<cellLID.size(); c++)
  {
    for (int n=0; n<nFacets_; n++)
    {
      for (int f=0; f<rangeDim_; f++)
      {
        int dof = dofs0[(c*rangeDim_ + f)*nFacets_ + n];
        (*elemToVecValuesMap_)[(cellLID[c]*nFacets_ + n)*rangeDim_ + f]
          = vec[dof];
      }
    }
  }
}


void CToAInterpolator::interpolate(const Teuchos::Array<double>& positions,
  Teuchos::Array<double>& results) const
{
  TimeMonitor timer(particleInterpolationTimer());

  TEUCHOS_TEST_FOR_EXCEPTION(positions.size() % dim_ != 0, std::runtime_error,
    "vector of coordinates should by an integer multiple "
    "of the spatial dimension");

  int nPts = positions.size() / dim_;

  results.resize(rangeDim_ * nPts);

  Array<double> xLocal(dim_);

  for (int i=0; i<nPts; i++)
  {
    const double* x = &(positions[dim_*i]);

    int guess = locator_.guessCell(x);

    TEUCHOS_TEST_FOR_EXCEPTION(guess < 0, std::runtime_error, "particle position "
      << x << " out of search grid");

      
    int cellLID = locator_.findEnclosingCell(guess, x, &(xLocal[0]));

    if (dim_==2)
    {
      double s = xLocal[0];
      double t = xLocal[1];
      Array<double> phi(nFacets_);
      phi[0] = 1.0-s-t;
      phi[1] = s;
      phi[2] = t;
      for (int f=0; f<rangeDim_; f++) results[rangeDim_*i+f] = 0.0;
      for (int n=0; n<nFacets_; n++)
      {
        for (int f=0; f<rangeDim_; f++)
        {
          results[rangeDim_*i+f] += phi[n]*(*elemToVecValuesMap_)[(cellLID*nFacets_ + n)*rangeDim_ + f];
        }
      }
    }
  }

}

