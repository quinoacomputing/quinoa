/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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


// $Id$ 
// $Source$ 


//   


#include "NOX_Playa_StatusTestBuilder.hpp"         
#include "NOX_StatusTest_NormF.H"         
#include "NOX_StatusTest_NormUpdate.H"         
#include "NOX_StatusTest_SafeCombo.hpp"         
#include "NOX_StatusTest_MaxIters.H"         
#include "Teuchos_Assert.hpp"   

using namespace NOX;
using namespace NOX::NOXPlaya;
using namespace Teuchos;
using std::runtime_error;

RCP<StatusTest::Generic> 
StatusTestBuilder::makeStatusTest(const ParameterList& params)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!params.isSublist("Status Test"), runtime_error,
                     "did not find Status Test sublist in " << params);

  ParameterList testSublist = params.sublist("Status Test");

  double fTol = 1.0e-15;
  double dxTol = 1.0e-15;
  int maxiters = 20;
  if (testSublist.isParameter("Tolerance"))
    {
      fTol = getParameter<double>(testSublist, "Tolerance");
    }
  if (testSublist.isParameter("Residual Tolerance"))
    {
      fTol = getParameter<double>(testSublist, "Residual Tolerance");
    }
  if (testSublist.isParameter("Step Tolerance"))
    {
      dxTol = getParameter<double>(testSublist, "Step Tolerance");
    }
  if (testSublist.isParameter("Max Iterations"))
    {
      maxiters = getParameter<int>(testSublist, "Max Iterations");
    }

  RCP<StatusTest::Generic> A = rcp(new StatusTest::NormF(fTol));
  RCP<StatusTest::Generic> B = rcp(new StatusTest::MaxIters(maxiters));
  RCP<StatusTest::Generic> C = rcp(new StatusTest::NormUpdate(dxTol));
  RCP<StatusTest::Generic> AB 
    = rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, A, B));
  RCP<StatusTest::Generic> ABC 
    = rcp(new StatusTest::SafeCombo(StatusTest::SafeCombo::OR, AB, C));
  
  return ABC;
}



