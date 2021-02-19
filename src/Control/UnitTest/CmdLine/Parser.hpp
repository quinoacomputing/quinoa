// *****************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Parser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's command line parser
  \details   This file declares the command-line argument parser for the unit
     test suite, UnitTest.
*/
// *****************************************************************************
#ifndef UnitTestCmdLineParser_h
#define UnitTestCmdLineParser_h

#include "StringParser.hpp"
#include "UnitTest/CmdLine/CmdLine.hpp"

namespace tk { class Print; }

namespace unittest {

//! \brief Command-line parser for UnitTest.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the unit test suite, UnitTest.
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline,
                            bool& helped );
};

} // unittest::

#endif // UnitTestCmdLineParser_h
