/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _test_FEI_Impl_hpp_
#define _test_FEI_Impl_hpp_

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

/** Unit-tester for fei::FEI_Impl.
*/
class test_FEI_Impl : public tester {
 public:
  test_FEI_Impl(MPI_Comm comm);
  virtual ~test_FEI_Impl();

  const char* getName()
    {
      static const char name[] = "fei::FEI_Impl";
      return((const char*)name);
    }

  int runtests();

  int serialtest1();

  int test1();

  int test2();

  int test3();

  int test4();

  int compareCoefs(int n,
		   const double*const* coefs1,
		   const double*const* coefs2);
};


#endif // _test_FEI_Impl_hpp_

