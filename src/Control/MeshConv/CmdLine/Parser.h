//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 04:05:34 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     MeshConv's command line parser
  \details   MeshConv's command line parser
*/
//******************************************************************************
#ifndef MeshConvCmdLineParser_h
#define MeshConvCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <MeshConv/CmdLine/CmdLine.h>

namespace meshconv {

//! CmdLineParser : StringParser
class CmdLineParser : public tk::StringParser {

  public:
    //! Constructor
    explicit CmdLineParser( int argc, char** argv,
                            const tk::Print& print,
                            ctr::CmdLine& cmdline );
};

} // meshconv::

#endif // MeshConvCmdLineParser_h
