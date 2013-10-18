//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Fri Oct 18 12:08:31 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's command line parser
  \details   Quinoa's command line parser
*/
//******************************************************************************
#ifndef QuinoaCmdLineParser_h
#define QuinoaCmdLineParser_h

#include <StringParser.h>
#include <Quinoa/CmdLine/CmdLine.h>

namespace quinoa {

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

} // quinoa::

#endif // QuinoaCmdLineParser_h
