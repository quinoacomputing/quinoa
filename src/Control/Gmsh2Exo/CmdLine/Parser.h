//******************************************************************************
/*!
  \file      src/Control/Gmsh2Exo/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Wed Mar 19 10:16:36 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh2Exo's command line parser
  \details   Gmsh2Exo's command line parser
*/
//******************************************************************************
#ifndef Gmsh2ExoCmdLineParser_h
#define Gmsh2ExoCmdLineParser_h

#include <Print.h>
#include <StringParser.h>
#include <Gmsh2Exo/CmdLine/CmdLine.h>

namespace gmsh2exo {

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

} // gmsh2exo::

#endif // Gmsh2ExoCmdLineParser_h
