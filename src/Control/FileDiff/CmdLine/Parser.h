// *****************************************************************************
/*!
  \file      src/Control/FileDiff/CmdLine/Parser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     FileDiff's command line parser
  \details   This file declares the command-line argument parser for the mesh
     file converter, FileDiff.
*/
// *****************************************************************************
#ifndef FileDiffCmdLineParser_h
#define FileDiffCmdLineParser_h

#include "StringParser.h"
#include "FileDiff/CmdLine/CmdLine.h"

namespace tk { class Print; }

namespace filediff {

//! \brief Command-line parser for FileDiff.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the file converter, FileDiff.
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // filediff::

#endif // FileDiffCmdLineParser_h
