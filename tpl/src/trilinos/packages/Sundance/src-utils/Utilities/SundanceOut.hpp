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

#ifndef SUNDANCE_OUT_H
#define SUNDANCE_OUT_H

#include "SundanceDefs.hpp"
#include "Teuchos_Assert.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "PlayaMPIComm.hpp"

using Playa::Out;
using Playa::Tabs;

#define SUNDANCE_OUT(test, msg) PLAYA_OUT(test, msg)


#define SUNDANCE_VERB_EXTREME(msg) PLAYA_MSG4(this->verb(), msg)
#define SUNDANCE_VERB_HIGH(msg) PLAYA_MSG3(this->verb(), msg)
#define SUNDANCE_VERB_MEDIUM(msg) PLAYA_MSG2(this->verb(), msg)
#define SUNDANCE_VERB_LOW(msg) PLAYA_MSG1(this->verb(), msg)

#define SUNDANCE_HEADER_LINE "\n------------------------------------------------------------------\n"

#define SUNDANCE_MSG(context, level, msg) PLAYA_OUT(this->verbLevel(context) >= level, msg)

#define SUNDANCE_LEVEL1(context, msg) PLAYA_MSG(context, 1, msg)

#define SUNDANCE_LEVEL2(context, msg) PLAYA_MSG(context, 2, msg)

#define SUNDANCE_LEVEL3(context, msg) PLAYA_MSG(context, 3, msg)

#define SUNDANCE_LEVEL4(context, msg) PLAYA_MSG(context, 4, msg)

#define SUNDANCE_LEVEL5(context, msg) PLAYA_MSG(context, 5, msg)


#define SUNDANCE_MSG1(level, msg) PLAYA_OUT(level >= 1, msg)

#define SUNDANCE_MSG2(level, msg) PLAYA_OUT(level >= 2, msg)

#define SUNDANCE_MSG3(level, msg) PLAYA_OUT(level >= 3, msg)

#define SUNDANCE_MSG4(level, msg) PLAYA_OUT(level >= 4, msg)

#define SUNDANCE_MSG5(level, msg) PLAYA_OUT(level >= 5, msg)

#define SUNDANCE_BANNER1(level, tab, msg) \
  PLAYA_MSG1(level, tab \
    << "===================================================================");\
  PLAYA_MSG1(level, tab << std::endl << tab \
    << "  " << msg); \
  PLAYA_MSG1(level, tab << std::endl << tab\
    << "===================================================================");


#define SUNDANCE_BANNER2(level, tab, msg) \
  SUNDANCE_MSG2(level, tab \
    << "-------------------------------------------------------------------");\
  SUNDANCE_MSG2(level, tab << msg); \
  SUNDANCE_MSG2(level, tab\
    << "-------------------------------------------------------------------");



#define SUNDANCE_BANNER3(level, tab, msg) \
  SUNDANCE_MSG3(level, tab \
    << "-------------------------------------------------------------------");\
  SUNDANCE_MSG3(level, tab << std::endl << tab \
    << msg); \
  SUNDANCE_MSG3(level, tab << std::endl << tab\
    << "-------------------------------------------------------------------");

#define SUNDANCE_TRACE(e) \
{ \
  TeuchosOStringStream omsg; \
        omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
        throw std::runtime_error(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define SUNDANCE_TRACE_MSG(e, msg)                      \
{ \
  TeuchosOStringStream omsg; \
        omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  omsg << msg << std::endl; \
  throw std::runtime_error(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}


#define SUNDANCE_ERROR(msg) \
{ \
  TeuchosOStringStream omsg; \
        omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  const std::string &omsgstr = omsg.str(); \
  Teuchos::TestForException_break(omsgstr); \
  throw std::runtime_error(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}


#endif
