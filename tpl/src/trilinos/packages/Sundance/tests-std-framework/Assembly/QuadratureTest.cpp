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

#include "Sundance.hpp"
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceFeketeTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace Playa;


int main(int argc, char** argv)
{
  int stat = 0 ;
	try
		{
			GlobalMPISession session(&argc, &argv);

      Array<int> validTetOrders = tuple(1, 2, 4, 6);
      int maxorder = 15;
      Array<int> validTriOrders(maxorder);
      for (int i=0; i<maxorder; i++)
    	  validTriOrders[i] = i+1;

      Array<int> validFeketeTriOrders = tuple(1, 2, 3, 4, 5, 6, 9);


      Array<int> triFailures;
      Array<int> tetFailures;
      Array<int> FeketeTriFailures;

      std::cerr << "------------- testing triangle rules -------------------"  << std::endl;
      for (int i=0; i<validTriOrders.size(); i++)
				{
          int p = validTriOrders[i];
					bool pass = TriangleQuadrature::test(p);
					if (pass) std::cerr << "order " << p << " PASSED" << std::endl;
					else 
            {
              std::cerr << "order " << p << " FAILED" << std::endl;
              triFailures.append(p);
            }
				}
      std::cerr << "------------- testing tet rules -------------------"  << std::endl;
          
      for (int i=0; i<validTetOrders.size(); i++)
				{
          int p = validTetOrders[i];
					bool pass = TetQuadrature::test(p);
					if (pass) std::cerr << "order " << p << " PASSED" << std::endl;
					else 
            {
              std::cerr << "order " << p << " FAILED" << std::endl;
              tetFailures.append(p);
            }
				}

      std::cerr << "--------- testing Fekete triangle rules ----------------" << std::endl;
      for (int i = 0; i < validFeketeTriOrders.size(); i++)
		{
			int p = validFeketeTriOrders[i];
			bool pass = FeketeTriangleQuadrature::test(p);
			if (pass)
				cerr << "order " << p << " PASSED" << std::endl;
			else
			{
				cerr << "order " << p << " FAILED" << std::endl;
				FeketeTriFailures.append(p);
			}
		}


      if (tetFailures.size()>0) 
        {
          cout << "failures detected for tets: orders " << tetFailures << std::endl;
          cout << "tet tests FAILED" << std::endl;
          stat = -1;
        }
      else
        {
          cout << "tet tests PASSED" << std::endl;
        }

      if (triFailures.size()>0) 
        {
          cout << "failures detected for tris: orders " << triFailures << std::endl;
          cout << "tri tests FAILED" << std::endl;
          stat = -1;
        }
      else
        {
          cout << "tri tests PASSED" << std::endl;
        }

		if (FeketeTriFailures.size() > 0)
		{
			cout << "failures detected for Fekete tris: orders "
					<< FeketeTriFailures << std::endl;
			cout << "Fekete tri tests FAILED" << std::endl;
			stat = -1;
		}
		else
		{
			cout << "Fekete tri tests PASSED" << std::endl;
		}

		}
	catch(std::exception& e)
		{
      std::cerr << "Detected exception: " << e.what() << std::endl;
      std::cerr << "Quadrature test FAILED" << std::endl;
      stat = -1;
		}

  return stat;
  
}
