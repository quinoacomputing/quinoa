//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Parser.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:17:48 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's command line parser
  \details   RNGTest's command line parser
*/
//******************************************************************************
#ifndef RNGTestCmdLineParser_h
#define RNGTestCmdLineParser_h

#include <StringParser.h>

namespace rngtest {

//! CmdLineParser : StringParser
class CmdLineParser : public tk::StringParser{

  public:
    //! Constructor from std::string
    explicit CmdLineParser(const std::string& cmdline, Base& base) :
      tk::StringParser(cmdline),
      m_base(base) {}

    //! Constructor from argc, argv
    explicit CmdLineParser(int argc, char** argv, Base& base) :
      tk::StringParser(argc, argv),
      m_base(base) {}

    //! Destructor
    ~CmdLineParser() noexcept override = default;

    //! Parse rngtest control file
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

    const Base& m_base;                  //!< Essentials
};

} // rngtest::

#endif // RNGTestCmdLineParser_h
