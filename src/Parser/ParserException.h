//******************************************************************************
/*!
  \file      src/Parser/ParserException.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:00:32 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ParserException
  \details   ParserException
*/
//******************************************************************************
#ifndef ParserException_h
#define ParserException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! Parser exception types
enum ParserExceptType { CMDLINE_EXCEPT=0,           //!< Command line exception
                        UNKNOWN_KEYWORD,            //!< Unknown keyword
                        UNKNOWN_HYDRO,              //!< Unknown hydro model
                        NUM_PARSER_EXCEPT
};

//! Parser exception error messages
const string ParserMsg[NUM_PARSER_EXCEPT] = {
  "Exactly one command line argument required: filename.q",
  "Unknown keyword",
  "Unknown hydrodynamics model: ",
};

//! ParserException : Exception
class ParserException : public Exception {

  public:
    //! Constructor without message
    explicit ParserException(const ExceptType except,
                             const ParserExceptType parserExcept,
                             const string& file,
                             const string& func,
                             const unsigned int& line) :
      Exception(except, file, func, line), m_except(parserExcept) {}

    //! Constructor with message from thrower
    explicit ParserException(const ExceptType except,
                             const ParserExceptType parserExcept,
                             const string throwerMsg,
                             const string& file,
                             const string& func,
                             const unsigned int& line) :
      Exception(except, file, func, line), m_throwerMsg(throwerMsg),
      m_except(parserExcept) {}

    //! Move constructor for throws, default compiler generated
    ParserException(ParserException&&) = default;

    //! Destructor
    virtual ~ParserException() noexcept = default;

    //! Handle ParserException
    virtual ErrCode handleException(Driver* const driver);

  protected:
    //! Message from thrower
    string m_throwerMsg;

  private:
    //! Don't permit copy constructor
    ParserException(const ParserException&) = delete;
    //! Don't permit copy assignment
    ParserException& operator=(const ParserException&) = delete;
    //! Don't permit move assignment
    ParserException& operator=(ParserException&&) = delete;

    //! Parser exception type (MIXMODEL, etc.)
    const ParserExceptType m_except;
};

} // namespace Quinoa

#endif // ParserException_h
