// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos {


TEUCHOS_UNIT_TEST( UnitTestHarness, nonRootFails ) {
  out << "Pass on even procs but fail on other procs!\n";
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  const int procRank = comm->getRank();
  TEST_EQUALITY_CONST(procRank%2, 0);
}


TEUCHOS_UNIT_TEST( UnitTestHarness, nonRootThrowsTeuchosExcept ) {
  out << "Pass on even procs but throws Teuchos exception on other processes!\n";
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  const int procRank = comm->getRank();
  TEUCHOS_ASSERT_EQUALITY(procRank%2, 0); // Throws on non-root processes
  int myval = 1;
  TEST_EQUALITY_CONST(myval, 1);
}


TEUCHOS_UNIT_TEST( UnitTestHarness, nonRootThrowsIntExcept ) {
  out << "Pass on even procs but throws int exception on other processes!\n";
  const RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  const int procRank = comm->getRank();
  if (procRank%2 != 0) {
    throw procRank;
  }
  int myval = 1;
  TEST_EQUALITY_CONST(myval, 1);
}


} // namespace Teuchos



