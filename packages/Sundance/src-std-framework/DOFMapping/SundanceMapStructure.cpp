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

#include "SundanceMapStructure.hpp"
#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceBasisDOFTopologyBase.hpp"
#include "SundanceOut.hpp"


namespace Sundance
{
using namespace Teuchos;
using std::endl;
using std::setw;

MapStructure::MapStructure(int nTotalFuncs,
                           const Array<RCP<BasisDOFTopologyBase> >& bases,
                           const Array<Array<int> >& funcs)
{
  init(nTotalFuncs, bases, funcs);
}


MapStructure::MapStructure(int nTotalFuncs,
                           const RCP<BasisDOFTopologyBase>& bases,
                           const Array<Array<int> >& funcs)
{
  init(nTotalFuncs, replicate(bases, funcs.size()), funcs);
}

MapStructure::MapStructure(int nTotalFuncs,
                           const RCP<BasisDOFTopologyBase>& bases)
{
  Array<int> f(nTotalFuncs);
  for (int i=0; i<nTotalFuncs; i++) f[i] = i;
  init(nTotalFuncs, tuple(bases), tuple(f));
}



void MapStructure::init(int nTotalFuncs,
                        const Array<RCP<BasisDOFTopologyBase> >& bases,
                        const Array<Array<int> >& funcs)
{
  bases_ = bases;
  funcs_ = funcs;
  chunkForFuncID_.resize(nTotalFuncs);
  indexForFuncID_.resize(nTotalFuncs);

  TEUCHOS_TEST_FOR_EXCEPTION(bases.size() != funcs.size(), std::logic_error,
                     "mismatched number of basis chunks=" << bases.size()
                     << " and number of function chunks=" << funcs.size());

  for (int f=0; f<indexForFuncID_.size(); f++) 
    {
      indexForFuncID_[f] = -1;
      chunkForFuncID_[f] = -1;
    }

  for (int b=0; b<funcs_.size(); b++)
    {
      for (int f=0; f<funcs_[b].size(); f++)
        {
          int fid = funcs_[b][f];
          TEUCHOS_TEST_FOR_EXCEPTION(fid >= nTotalFuncs, std::logic_error,
                             "bad funcID=" << fid 
                             << ". nTotalFuncs=" << nTotalFuncs);
          indexForFuncID_[fid] = f;
          chunkForFuncID_[fid] = b;
        }
    }
}

int MapStructure::chunkForFuncID(int funcID) const 
{
  int rtn = chunkForFuncID_[funcID];
  TEUCHOS_TEST_FOR_EXCEPTION(rtn < 0, std::logic_error,
                     "funcID=" << funcID << " not defined in map chunk."
                     " The functions defined there are " 
                     << funcs_ << ". The most likely cause of this error is "
                     "that you are trying to access a discrete function on "
                     "subdomain for which it is not defined.");
  return rtn;
}

int MapStructure::indexForFuncID(int funcID) const 
{
  int rtn = indexForFuncID_[funcID];
  TEUCHOS_TEST_FOR_EXCEPTION(rtn < 0, std::logic_error,
                     "funcID=" << funcID << " not defined in map chunk."
                     " The functions defined there are " 
                     << funcs_ << ". The most likely cause of this error is "
                     "that you are trying to access a discrete function on "
                     "subdomain for which it is not defined.");

  return rtn;
}


std::ostream& MapStructure::print(std::ostream& os) const
{
  Tabs tab;
  os << tab << "Map structure: nChunks=" << numBasisChunks() << std::endl;
  for (int c=0; c<numBasisChunks(); c++) 
  {
    Tabs tab1;
    os << tab1 << "chunk " << c << " funcIDs=" << funcs(c) << std::endl;
  }
  os << tab << "## end map structure" << std::endl;
  return os;
}


Array<RCP<BasisDOFTopologyBase> > replicate(
  const RCP<BasisDOFTopologyBase>& model,
  int n)
{
  Array<RCP<BasisDOFTopologyBase> > rtn(n);
  for (int i=0; i<n; i++) rtn[i] = model;
  return rtn;
}

} // end namespace Sundance

