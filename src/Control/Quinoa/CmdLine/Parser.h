//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Thu Oct  3 07:39:19 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's command line parser
  \details   Quinoa's command line parser
*/
//******************************************************************************
#ifndef QuinoaCmdLineParser_h
#define QuinoaCmdLineParser_h

#include <StringParser.h>

namespace quinoa {

//! CmdLineParser : StringParser
class CmdLineParser : public StringParser{

  public:
    //! Constructor from std::string
    explicit CmdLineParser(const std::string& cmdline, Base& base) :
      StringParser(cmdline, base) {}

    //! Constructor from argc, argv
    explicit CmdLineParser(int argc, char** argv, Base& base) :
      StringParser(argc, argv, base) {}

    //! Destructor
    ~CmdLineParser() noexcept override = default;

    //! Parse quinoa control file
    void parse() override;

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

} // namespace quinoa

#endif // QuinoaCmdLineParser_h
