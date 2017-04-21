// *****************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Parser.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     MeshConv's command line parser
  \details   This file declares the command-line argument parser for the mesh
     file converter, MeshConv.
*/
// *****************************************************************************
#ifndef MeshConvCmdLineParser_h
#define MeshConvCmdLineParser_h

#include "StringParser.h"
#include "MeshConv/CmdLine/CmdLine.h"

namespace tk { class Print; }

namespace meshconv {

//! \brief Command-line parser for MeshConv.
//! \details This class is used to interface with PEGTL, for the purpose of
//!   parsing command-line arguments for the mesh converter, MeshConv.
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc,
                            char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // meshconv::

#endif // MeshConvCmdLineParser_h
