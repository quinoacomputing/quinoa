//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Tue 08 Apr 2014 09:27:41 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
class CmdLineParser : public tk::StringParser{

  public:
    //! Constructor
    explicit CmdLineParser(int argc, char** argv,
                           const tk::Print& print,
                           std::unique_ptr< ctr::CmdLine >& cmdline);

    //! Destructor
    ~CmdLineParser() noexcept override = default;

  private:
    //! Don't permit copy constructor
    CmdLineParser(const CmdLineParser&) = delete;
    //! Don't permit copy assigment
    CmdLineParser& operator=(const CmdLineParser&) = delete;
    //! Don't permit move constructor
    CmdLineParser(CmdLineParser&&) = delete;
    //! Don't permit move assigment
    CmdLineParser& operator=(CmdLineParser&&) = delete;
};

} // meshconv::

#endif // MeshConvCmdLineParser_h
