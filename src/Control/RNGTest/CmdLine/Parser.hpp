// *****************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Parser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     RNGTest's command line parser
  \details   This file declares the command-line argument parser for the random
     number generator test suite, RNGTest.
*/
// *****************************************************************************
#ifndef RNGTestCmdLineParser_h
#define RNGTestCmdLineParser_h

#include "StringParser.hpp"
#include "RNGTest/CmdLine/CmdLine.hpp"

namespace tk { class Print; }

namespace rngtest {

//! \brief Command-line parser for RNGTest.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the random number generator test suite,
//!   RNGTest.
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // rngtest::

#endif // RNGTestCmdLineParser_h
