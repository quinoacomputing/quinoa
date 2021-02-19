// *****************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Parser.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     MeshConv's command line parser
  \details   This file declares the command-line argument parser for the mesh
     file converter, MeshConv.
*/
// *****************************************************************************
#ifndef MeshConvCmdLineParser_h
#define MeshConvCmdLineParser_h

#include "StringParser.hpp"
#include "MeshConv/CmdLine/CmdLine.hpp"

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
