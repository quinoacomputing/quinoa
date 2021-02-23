// *****************************************************************************
/*!
  \file      src/Control/Inciter/CmdLine/Parser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's command line parser
  \details   This file declares the command-line argument parser for the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#ifndef InciterCmdLineParser_h
#define InciterCmdLineParser_h

#include "StringParser.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"

namespace tk { class Print; }

namespace inciter {

//! \brief Command-line parser for Inciter.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the computational shock hydrodynamics
//!   tool, Inciter.
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc, char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // inciter::

#endif // InciterCmdLineParser_h
