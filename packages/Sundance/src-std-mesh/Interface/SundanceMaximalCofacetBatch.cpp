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

#include "SundanceMaximalCofacetBatch.hpp"
#include "PlayaExceptions.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;

MaximalCofacetBatch::MaximalCofacetBatch()
  : cofacetLIDs_(2),
    facetIndices_(2),
    numCells_(-1),
    numCofacets_(-1)
{
  for (int i=0; i<2; i++)
  {
    cofacetLIDs_[i] = rcp(new Array<int>());
    facetIndices_[i] = rcp(new Array<int>());
  }
}

void MaximalCofacetBatch::reset(int numCells, int numCofacets)
{
  TEUCHOS_TEST_FOR_EXCEPTION(numCofacets < 0 || numCofacets > 2, std::runtime_error,
    "invalid number of maximal cofacets = " << numCofacets);
  numCofacets_ = numCofacets;
  reset(numCells);
}

void MaximalCofacetBatch::reset(int numCells)
{
  numCells_ = numCells;
  TEUCHOS_TEST_FOR_EXCEPTION(numCells_ <= 0, std::runtime_error,
    "invalid number of cells = " << numCells_);
  for (int i=0; i<numCofacets_; i++)
  {
    cofacetLIDs_[i]->resize(numCells_);
    facetIndices_[i]->resize(numCells_);
  }
}

void MaximalCofacetBatch::addSingleCofacet(int c, 
  int cofacetLID, int facetIndex)
{
  TEUCHOS_TEST_FOR_EXCEPTION(numCofacets_ != 1, std::runtime_error,
    "addSingleCofacet called for a batch configured with " << numCofacets_ 
    << " cofacets");

  cofacetLIDs_[0]->operator[](c) = cofacetLID;
  facetIndices_[0]->operator[](c) = facetIndex;
}
void MaximalCofacetBatch::addTwoCofacets(int c, 
  int cofacet1, int facetIndex1,
  int cofacet2, int facetIndex2
)
{
  TEUCHOS_TEST_FOR_EXCEPTION(numCofacets_ != 2, std::runtime_error,
    "addTwoCofacets called for a batch configured with " << numCofacets_ 
    << " cofacets");

  cofacetLIDs_[0]->operator[](c) = cofacet1;
  facetIndices_[0]->operator[](c) = facetIndex1;

  cofacetLIDs_[1]->operator[](c) = cofacet2;
  facetIndices_[1]->operator[](c) = facetIndex2;
}

int MaximalCofacetBatch::cofacetLID(int c, int n, int& facetIndex) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(n >= numCofacets_, std::runtime_error,
    "invalid cofacet number n=" << n);
  TEUCHOS_TEST_FOR_EXCEPTION(c >= numCells_, std::runtime_error,
    "invalid cell number c=" << c);

  facetIndex = facetIndices_[n]->operator[](c);
  return cofacetLIDs_[n]->operator[](c);
}

void MaximalCofacetBatch::getSpecifiedCofacets(
  const Array<int>& cofacetNumbers,
  RCP<Array<int> >& cofacets,
  RCP<Array<int> >& facetIndices) const
{
  TEUCHOS_TEST_FOR_EXCEPTION((int) cofacetNumbers.size() != numCells(),
    std::runtime_error,
    "mismatch between cofacet batch size (" << numCells() << ") and "
    "requested number of cofacets (" << cofacetNumbers.size() << ")");

  cofacets->resize(numCells());
  facetIndices->resize(numCells());

  for (int c=0; c<numCells(); c++)
  {
    (*cofacets)[c] = cofacetLID(c, cofacetNumbers[c], (*facetIndices)[c]);
  }
}

void MaximalCofacetBatch::getSpecifiedCofacets(
  int cofacetNumber,
  RCP<Array<int> >& cofacets,
  RCP<Array<int> >& facetIndices) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(cofacetNumber < 0 || cofacetNumber>1,    
    std::runtime_error,
    "invalid cofacet number=" << cofacetNumber);

  cofacets = cofacetLIDs_[cofacetNumber];
  facetIndices = facetIndices_[cofacetNumber];
}
