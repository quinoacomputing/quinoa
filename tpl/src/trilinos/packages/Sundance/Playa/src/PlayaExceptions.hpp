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

#ifndef PLAYA_EXCEPTIONS_H
#define PLAYA_EXCEPTIONS_H

#include "PlayaDefs.hpp"
#include "PlayaDebug.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"
#include <stdexcept>


//bvbw for backard compatibility reasons
//     I could not get this to work with ifdefs hence the hack

#define PLAYA_ERROR7(msg) \
{ \
  Teuchos::TestForException_break(); \
  std::ostringstream omsg; \
	omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define PLAYA_ERROR(msg) \
{ \
  std::ostringstream omsg; \
	omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  const std::string &omsgstr = omsg.str(); \
  Teuchos::TestForException_break(omsgstr); \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}


#define PLAYA_TRACE(e) \
{ \
  std::ostringstream omsg; \
	omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define PLAYA_TRACE_MSG(e, msg)                      \
{ \
  std::ostringstream omsg; \
	omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  omsg << msg << std::endl; \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define PLAYA_BOUNDSCHECK(i, low, high, location) \
{ \
  TEUCHOS_TEST_FOR_EXCEPTION( i < low || i > high, Playa::RuntimeError, \
    "Bounds violation in " << location << ": " \
    << #i << "is out of range [" \
    << #low << ", " << #high << "]") \
}

#define PLAYA_CHECK_ARRAY_SIZE_MATCH(a1, a2) \
  {\
    TEUCHOS_TEST_FOR_EXCEPTION(a1.size() != a2.size(), Playa::RuntimeError, \
y      "Mismatched array sizes: size(" << #a1 << ")=" << a1.size() \
      << " size(" << #a2 << ")=" << a2.size() << ". Expected equal sizes");\
  }



namespace Playa
{
  /**
   * InternalError is thrown when an "impossible" case is detected
   * in Playa's internals. Occurance of an InternalError indicates
   * either a bug in Playa or runtime memory corruption that is
   * scrambling an object's virtual tables.
   */
  class InternalError : public std::logic_error
    {
    public:
      /** */
      InternalError(const std::string& msg);
    };

  /**
   * RuntimeError is an exception that occurs as a result of invalid
   * user-level code.
   */
  class RuntimeError : public std::runtime_error
    {
    public:
      /** */
      RuntimeError(const std::string& msg);
    };
  
}

                  

#endif
