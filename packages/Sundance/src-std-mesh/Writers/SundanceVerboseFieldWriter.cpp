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

#include "SundanceVerboseFieldWriter.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"


using namespace Sundance;
using namespace Teuchos;


void VerboseFieldWriter::write() const 
{
  int nProc = mesh().comm().getNProc();
  int myRank = mesh().comm().getRank();

  RCP<std::ostream> osp;
  if (filename().length()==0)
    {
      osp = rcp(&std::cout, false);
    }
  else 
    {
      std::string f = filename() + ".txt";
      if (nProc > 1) f = f + "." + Teuchos::toString(myRank);
      osp = rcp(new std::ofstream(f.c_str()));
    }
  std::ostream& os = *osp;

  if (myRank==0) os << "VerboseFieldWriter output" << std::endl;
  for (int p=0; p<nProc; p++)
    {
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      if (p != myRank) continue;
      os << "======== processor " << p << " ============================ "
         << std::endl;
      Tabs tab0;
      int dim = mesh().spatialDim();
      int nPts = mesh().numCells(0);
      int nElems = mesh().numCells(dim);
      os << tab0 << "spatial dimension = " << dim << std::endl;
      os << tab0 << "num points = " << nPts << std::endl;
      os << tab0 << "num elements = " << nElems << std::endl;
      os << tab0 << "Point list: " << std::endl;

      int dummy;
      for (int i=0; i<nPts; i++)
        {
          Tabs tab1;
          os << tab1 << "L=" << i 
             << " G=" << mesh().mapLIDToGID(0, i) 
             << " x=" << mesh().nodePosition(i) 
             << " owner=" << mesh().ownerProcID(0,i) 
             << " label=" << mesh().label(0,i) << std::endl;
          int nc = mesh().numMaxCofacets(0,i);
          Tabs tab2;
          os << tab2 << "num cofacets=" << nc << " cofs = {";
          for (int c=0; c<nc; c++)
            {
              if (c==0) os << " " ;
              else os << ", ";
              os << mesh().mapLIDToGID(dim, mesh().maxCofacetLID(0,i,c,dummy));
            }
          os << "}" << std::endl;
        }

      
      os << tab0 << "Element list: " << std::endl;

      for (int i=0; i<nElems; i++)
        {
          int facetSign;
          Tabs tab1;
          os << tab1 << "L=" << i 
             << " G=" << mesh().mapLIDToGID(dim, i) 
             << ", nodes L={";
          int numNodes = mesh().numFacets(dim, i, 0);
          for (int j=0; j<numNodes; j++)
            {
              if (j != 0) os << ", ";
              os << mesh().facetLID(dim, i, 0, j, facetSign);
            }
          os << "}, G={";
          for (int j=0; j<numNodes; j++)
            {
              if (j != 0) os << ", ";
              os << mesh().mapLIDToGID(0, mesh().facetLID(dim, i, 0, j, facetSign));
            }
          os << "}, owner=" << mesh().ownerProcID(dim,i)
             << ", label=" << mesh().label(dim,i) << std::endl;
          for (int fd=1; fd<dim; fd++)
            {
              Tabs tab2;
              os << tab2 << "facets of dimension " << fd << std::endl;
              int nf = mesh().numFacets(dim, i, fd);
              for (int f=0; f<nf; f++)
                {

                  Tabs tab3;
                  int flid = mesh().facetLID(dim, i, fd, f, facetSign);
                  int fgid = -1;
                  int fowner = -1;
                  if (mesh().hasIntermediateGIDs(fd))
                    { 
                      fgid = mesh().mapLIDToGID(fd, flid);
                      fowner = mesh().ownerProcID(fd, flid);
                    }
                  os << tab3 << "f#=" << f << " L=" << flid
                     << " G=" << fgid << " owner=" << fowner
                     << " nodes={";
                  int nfn = mesh().numFacets(fd, flid, 0);
                  for (int fn=0; fn<nfn; fn++)
                    {
                      if (fn != 0) os << ", ";
                      os << mesh().facetLID(fd, flid, 0, fn, facetSign);
                    }
                  os << "} sign=" << facetSign;
                  os << " label=" << mesh().label(fd,flid) << std::endl;
                }
              
            }
        }
      
    }
}



