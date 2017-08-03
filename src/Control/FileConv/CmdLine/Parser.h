// *****************************************************************************
/*!
  \file      src/Control/FileConv/CmdLine/Parser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     FileConv's command line parser
  \details   This file declares the command-line argument parser for the mesh
     file converter, FileConv.
*/
// *****************************************************************************
#ifndef FileConvCmdLineParser_h
#define FileConvCmdLineParser_h

#include "StringParser.h"
#include "FileConv/CmdLine/CmdLine.h"

namespace tk { class Print; }

namespace fileconv {

//! \brief Command-line parser for FileConv.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the file converter, FileConv.
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // fileconv::

#endif // FileConvCmdLineParser_h
