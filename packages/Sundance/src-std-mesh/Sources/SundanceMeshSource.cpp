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


#include "SundanceMeshSource.hpp"
#include "SundanceOut.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMeshType.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Teuchos;
using namespace Sundance;
using Playa::Handle;
using Playa::Handleable;


static Time& getMeshTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("get mesh"); 
  return *rtn;
}

MeshSource::MeshSource()
  : Handle<MeshSourceBase>()
{}

MeshSource::MeshSource(Handleable<MeshSourceBase>* rawPtr)
  : Handle<MeshSourceBase>(rawPtr)
{}


MeshSource::MeshSource(const RCP<MeshSourceBase>& smartPtr)
  : Handle<MeshSourceBase>(smartPtr)
{}

Mesh MeshSource::getMesh() const
{
  TimeMonitor timer(getMeshTimer());

  Mesh rtn;
  try
    {
      Tabs tabs;
      int nProc = ptr()->comm().getNProc();
      SUNDANCE_OUT(ptr()->verb() > 0, 
                   "MeshSource::getMesh()");
      if (staggerOutput() && nProc > 1)
        {
          int myRank = ptr()->comm().getRank();
          for (int p=0; p<nProc; p++)
            {
              ptr()->comm().synchronize();
              if (p != myRank) continue;
              SUNDANCE_OUT(ptr()->verb() > 0, 
                           "========= Building local mesh on processor " 
                           << p << " ============ ");
              rtn = ptr()->getMesh();
            }
        }
      else 
        {
          rtn = ptr()->getMesh();
        }

      if (rtn.spatialDim() > 1) rtn.assignIntermediateCellGIDs(1);
      if (rtn.spatialDim() > 2) rtn.assignIntermediateCellGIDs(2);
    }
  catch(std::exception& e)
    {
      SUNDANCE_TRACE(e);
    }
  return rtn;
}

void MeshSource::getAttributes(RCP<Array<Array<double> > >& nodeAttributes,
                               RCP<Array<Array<double> > >& elemAttributes) const
{
  getMesh();
  ptr()->getAttributes(nodeAttributes, elemAttributes);
}


MeshType& MeshSource::defaultMeshType() 
{
  static MeshType rtn = new BasicSimplicialMeshType(); 
  return rtn;
}

const MPIComm& MeshSource::comm() const
{
  return ptr()->comm();
}
