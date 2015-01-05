//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 06:17:27 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     RNGTest's command line parser
  \details   This file declares the command-line argument parser for the random
     number generator test suite, RNGTest.
*/
//******************************************************************************
#ifndef RNGTestCmdLineParser_h
#define RNGTestCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <RNGTest/CmdLine/CmdLine.h>

namespace rngtest {

//! \brief Command-line parser for RNGTest.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the random number generator test suite,
//!   RNGTest.
//! \author J. Bakosi
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
