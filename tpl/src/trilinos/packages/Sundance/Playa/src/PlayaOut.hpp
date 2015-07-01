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

#ifndef PLAYA_OUT_HPP
#define PLAYA_OUT_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "PlayaMPIComm.hpp"

namespace Teuchos
{
template <class T> class Array;
}

namespace Playa
{
class Tabs;
using namespace Teuchos;
  using std::cout;
  using std::endl;
  using std::cerr;

/**
 * Class Out provides standardized access to fancy streams for writing
 * diagnostic output. 
 */
class Out
{
public:
      
  /** Print a line followed by termination */
  static void println(const std::string& str) 
    {
      if (!suppressStdout()) os() << str << std::endl;
    }

  /** */
  static std::ostream*& basicStream()
    {
      static std::ostream* rtn = &std::cout;
      return rtn;
    }

  /** Provide a fancy ostream wrapper to cout */
  static FancyOStream& os();

  /** Provide a fancy ostream wrapper to cout, ignoring output from
   * non-root processors*/
  static FancyOStream& root();

  /** */
  static bool& suppressStdout() {static bool rtn=false; return rtn;}
};


/** Write an array formatted to show a specified number of columns */
void writeTable(std::ostream& os, const Tabs& tab,
  const Array<double>& a, int cols);
/** Write an array formatted to show a specified number of columns */
void writeTable(std::ostream& os, const Tabs& tab, 
  const Array<int>& a, int cols);

}

#define PLAYA_OUT(test, msg) \
  { \
    if (test) \
    { \
      TeuchosOStringStream omsg; \
      omsg << msg; \
      Out::println(omsg.str());                 \
    } \
  }

#define PLAYA_ROOT_OUT(test, msg) \
  { \
    if (test) \
    { \
      TeuchosOStringStream omsg; \
      omsg << msg; \
      Out::root() << omsg.str() << std::endl;                 \
    } \
  }


#define PLAYA_VERB_EXTREME(msg) PLAYA_MSG4(this->verb(), msg)
#define PLAYA_VERB_HIGH(msg) PLAYA_MSG3(this->verb(), msg)
#define PLAYA_VERB_MEDIUM(msg) PLAYA_MSG2(this->verb(), msg)
#define PLAYA_VERB_LOW(msg) PLAYA_MSG1(this->verb(), msg)

#define PLAYA_HEADER_LINE "\n------------------------------------------------------------------\n"

#define PLAYA_MSG(context, level, msg) PLAYA_OUT(this->verbLevel(context) >= level, msg)

#define PLAYA_LEVEL1(context, msg) PLAYA_MSG(context, 1, msg)

#define PLAYA_LEVEL2(context, msg) PLAYA_MSG(context, 2, msg)

#define PLAYA_LEVEL3(context, msg) PLAYA_MSG(context, 3, msg)

#define PLAYA_LEVEL4(context, msg) PLAYA_MSG(context, 4, msg)

#define PLAYA_LEVEL5(context, msg) PLAYA_MSG(context, 5, msg)


#define PLAYA_MSG1(level, msg) PLAYA_OUT(level >= 1, msg)

#define PLAYA_MSG2(level, msg) PLAYA_OUT(level >= 2, msg)

#define PLAYA_MSG3(level, msg) PLAYA_OUT(level >= 3, msg)

#define PLAYA_MSG4(level, msg) PLAYA_OUT(level >= 4, msg)

#define PLAYA_MSG5(level, msg) PLAYA_OUT(level >= 5, msg)



#define PLAYA_ROOT_MSG1(level, msg) PLAYA_ROOT_OUT(level >= 1, msg)

#define PLAYA_ROOT_MSG2(level, msg) PLAYA_ROOT_OUT(level >= 2, msg)

#define PLAYA_ROOT_MSG3(level, msg) PLAYA_ROOT_OUT(level >= 3, msg)

#define PLAYA_ROOT_MSG4(level, msg) PLAYA_ROOT_OUT(level >= 4, msg)

#define PLAYA_ROOT_MSG5(level, msg) PLAYA_ROOT_OUT(level >= 5, msg)




#define PLAYA_BANNER1(level, tab, msg) \
  PLAYA_ROOT_MSG1(level, tab \
    << "===================================================================");\
  PLAYA_ROOT_MSG1(level, tab << std::endl << tab \
    << "  " << msg); \
  PLAYA_ROOT_MSG1(level, tab << std::endl << tab\
    << "===================================================================");


#define PLAYA_BANNER2(level, tab, msg) \
  PLAYA_ROOT_MSG2(level, tab \
    << "-------------------------------------------------------------------");\
  PLAYA_ROOT_MSG2(level, tab << msg); \
  PLAYA_MSG2(level, tab\
    << "-------------------------------------------------------------------");



#define PLAYA_BANNER3(level, tab, msg) \
  PLAYA_ROOT_MSG3(level, tab \
    << "-------------------------------------------------------------------");\
  PLAYA_ROOT_MSG3(level, tab << std::endl << tab \
    << msg); \
  PLAYA_ROOT_MSG3(level, tab << std::endl << tab\
    << "-------------------------------------------------------------------");

#endif
