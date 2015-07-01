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


#include "SundanceExpr.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceHermiteSpectralBasis.hpp"
#include "SundanceStokhosBasisWrapper.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifdef HAVE_SUNDANCE_STOKHOS
#include "Stokhos_HermiteBasis.hpp"
#endif

using std::cout;
using std::exception;
using std::setw;
using namespace Sundance;
using namespace Teuchos;
#ifdef HAVE_SUNDANCE_STOKHOS
using namespace Stokhos;
#endif


static Time& totalTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}



int main(int argc, char** argv)
{
#ifdef HAVE_SUNDANCE_STOKHOS
  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());

    SpectralBasis h1 = new HermiteSpectralBasis(1, 4);
    SpectralBasis h2 = new Stokhos::HermiteBasis<int, double>(4);

    bool fail = false;
    for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++)
      {
        for (int k=0; k<4; k++)
        {
          double c1 = h1.expectation(i,j,k);
          double c2 = h2.expectation(i,j,k);
          double err = fabs(c1 - c2);
          cout << setw(4) << i << setw(4) << j << setw(4) << k
               << setw(16) << c1 << setw(16) << c2 << setw(16) << err;
          if (err > 1.0e-12) 
          {
            cout << " ***** FAILED!" ;
            fail = true;
          }
          cout << std::endl;
        }
      }
    }

    
    TimeMonitor::summarize();
    if (fail) 
    {
      cout << "PCE test FAILED" << std::endl;
      return -1;
    }
    else
    {
      cout << "PCE test PASSED" << std::endl;
    }
  }
	catch(std::exception& e)
  {
    Out::println(e.what());
    return -1;
  }
#else
  std::cout << "test disabled because Stokhos has not been enabled" << std::endl;

#endif
  return 0;  
}

