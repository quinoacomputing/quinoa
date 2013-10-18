//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Fri Oct 18 12:11:01 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's command line parser
  \details   RNGTest's command line parser
*/
//******************************************************************************
#ifndef RNGTestCmdLineParser_h
#define RNGTestCmdLineParser_h

#include <StringParser.h>
#include <RNGTest/CmdLine/CmdLine.h>

namespace rngtest {

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

} // rngtest::

#endif // RNGTestCmdLineParser_h
