//******************************************************************************
/*!
  \file      src/Control/FileParser.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:11:30 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     File parser
  \details   File parser
*/
//******************************************************************************
#ifndef FileParser_h
#define FileParser_h

#include <string>

#include <Parser.h>

namespace tk {

//! FileParser : Parser
class FileParser : public Parser {

  protected:
    //! Constructor
    explicit FileParser(const std::string& filename);

    //! Destructor
    ~FileParser() noexcept override = default;

    //! Parser interface
    void parse() override = 0;

    const std::string m_filename;             //!< Name of file to parse

  private:
    //! Don't permit copy constructor
    FileParser(const FileParser&) = delete;
    //! Don't permit copy assigment
    FileParser& operator=(const FileParser&) = delete;
    //! Don't permit move constructor
    FileParser(FileParser&&) = delete;
    //! Don't permit move assigment
    FileParser& operator=(FileParser&&) = delete;
};

} // tk::

#endif // FileParser_h
