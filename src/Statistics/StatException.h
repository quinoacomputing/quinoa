//******************************************************************************
/*!
  \file      src/Statistics/StatException.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:01:31 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics exception
  \details   Statistics Exception
*/
//******************************************************************************
#ifndef StatException_h
#define StatException_h

#include <string>

using namespace std;

#include <QuinoaTypes.h>
#include <Exception.h>

namespace Quinoa {

//! StatException types
enum StatExceptType {
  STATEXCEPT_UNIMPLEMENTED=0,  //!< function unimplemented
  STATEXCEPT_BAD_SAMPLE,       //!< Wrong sample space dimension
  STATEXCEPT_NO_SUCH_MOMENT,   //!< Wrong moment index
  STATEXCEPT_NO_SUCH_MEAN,     //!< Cannot find mean for fluctuation
  NUM_STAT_EXCEPT
};

//! StatException error messages
const string StatMsg[NUM_STAT_EXCEPT] = {
  "Method unimplemented",
  "Sample incompatible with sample space",
  "No such moment",
  "Cannot find mean: "
};

//! StatException : Exception
class StatException : public Exception {

  public:
    //! Constructor without message
    explicit StatException(const ExceptType except,
                           const StatExceptType statExcept,
                           const string& file,
                           const string& func,
                           const unsigned int& line) :
      Exception(except, file, func, line), m_except(statExcept) {}

    //! Constructor with message from thrower
    explicit StatException(const ExceptType except,
                           const StatExceptType statExcept,
                           const string throwerMsg,
                           const string& file,
                           const string& func,
                           const unsigned int& line) :
      Exception(except, file, func, line),
      m_throwerMsg(throwerMsg),
      m_except(statExcept) {}

    //! Move constructor for throws, default compiler generated
    StatException(StatException&&) = default;

    //! Destructor
    virtual ~StatException() noexcept = default;

    //! Handle StatException
    virtual ErrCode handleException(Driver* const driver);

  protected:
    //! Message from thrower
    string m_throwerMsg;

  private:
    //! Don't permit copy constructor
    StatException(const StatException&) = delete;
    //! Don't permit copy assignment
    StatException& operator=(const StatException&) = delete;
    //! Don't permit move assignment
    StatException& operator=(StatException&&) = delete;

    //! Statistrics exception type (UNIMPLEMENTED, etc.)
    const StatExceptType m_except;
};

} // namespace Quinoa

#endif // StatException_h
