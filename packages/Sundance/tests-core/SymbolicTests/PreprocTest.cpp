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


#include "SundanceSymbPreprocessor.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "Teuchos_TestingHelpers.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using std::cout;
using std::exception;

using Sundance::List;


#define TEST_THROW(code, passFail) \
  TEUCHOS_TEST_THROW( code, std::exception, Out::os(), passFail)

#define TEST_NOTHROW(code, passFail) \
  TEUCHOS_TEST_NOTHROW( code, Out::os(), passFail)

bool validateFuncTypeChecking()
{
  Expr ux = new UnknownFunctionStub("ux");
  Expr vx = new TestFunctionStub("vx");
  Expr uy = new UnknownFunctionStub("uy");
  Expr vy = new TestFunctionStub("vy");
  Expr uz = new UnknownFunctionStub("uz");
  Expr vz = new TestFunctionStub("vz");

  Expr v = List(vx,vy,vz);
  Expr u = List(ux,uy,uz);

  Expr mixup = List(vx, uy, vz); // mix of test & unknown
  Expr dup = List(vx, vx, vz); // list with duplicates

  bool passFail = true;
  /* */
  Out::os() << "Testing detection of mixed-up function types" << std::endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<UnknownFuncElement>(mixup, makeZeros(v)), passFail);
  
  /* */
  Out::os() << "Testing detection of duplicated functions" << std::endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(dup, makeZeros(v)), passFail);
  
  /* */
  Out::os() << "Testing detection of invalid evaluation points" << std::endl;
  TEST_THROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(v, u), passFail);
  

  /* */
  Out::os() << "Testing processing of good input" << std::endl;
  TEST_NOTHROW(
    SymbPreprocessor::processInputFuncs<SymbolicFuncElement>(v, makeZeros(v)), passFail);
  
  return passFail;
} 

int main(int argc, char** argv)
{
  bool pass = true;
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      pass = pass && validateFuncTypeChecking();
    }
	catch(std::exception& e)
		{
      pass = false;
			Out::println(e.what());
		}

  if (pass)
  {
    Out::os() << "test PASSED" << std::endl;
    return 0;
  }
  else 
  {
    Out::os() << "test FAILED" << std::endl;
    return -1;
  }

  
}
