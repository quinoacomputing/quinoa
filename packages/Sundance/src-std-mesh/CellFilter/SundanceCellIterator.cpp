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

#include "SundanceCellIterator.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

CellIterator::CellIterator()
  :  isImplicit_(true),
    currentLID_(-1),
    reorderer_(0),
     iter_(dummy().begin())
{;}

CellIterator::CellIterator(const CellIterator& other)
  :  isImplicit_(other.isImplicit_),
     currentLID_(other.currentLID_),
     reorderer_(other.reorderer_),
     iter_()
{
  if (!isImplicit_) iter_ = other.iter_;
}

CellIterator::CellIterator(const Mesh& mesh, 
                           int cellDim, 
                           CellIteratorPos pos)
  : isImplicit_(true),
    currentLID_(-1),
    reorderer_(mesh.reorderer()),
    iter_(dummy().begin())
{
  if (cellDim == mesh.spatialDim() && reorderer_ != 0)
    {
      switch(pos)
        {
        case Begin:
          currentLID_ = reorderer_->begin(); 
          break;
        case End:
          currentLID_ = mesh.numCells(cellDim);
        }
      SUNDANCE_OUT(mesh.verb() > 2, 
                   "created implicit cell iterator with LID=" << currentLID_);
    }
  else 
    {
      switch(pos)
        {
        case Begin:
          currentLID_ = 0;
          break;
        case End:
          currentLID_ = mesh.numCells(cellDim);
        }
      SUNDANCE_OUT(mesh.verb() > 2, 
                   "created implicit cell iterator with LID=" << currentLID_);
    }


}



CellIterator::CellIterator(const Set<int>* cells, CellIteratorPos pos)
  : isImplicit_(false),
    currentLID_(-1),
    reorderer_(0),
    iter_(dummy().begin())
{
  switch(pos)
  {
    case Begin:
      iter_ = cells->begin();
      break;
    case End:
      iter_ = cells->end();
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(1);
  }
}



CellIterator& CellIterator::operator=(const CellIterator& other)
{
  if (*this!=other) 
  {
    isImplicit_ = other.isImplicit_;
    currentLID_ = other.currentLID_;
    reorderer_=other.reorderer_;
    if (!isImplicit_) iter_ = other.iter_;
  }
  return *this;
}

    
