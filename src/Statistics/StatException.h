//******************************************************************************
/*!
  \file      src/Statistics/StatException.h
  \author    J. Bakosi
  \date      Sat 27 Oct 2012 11:39:08 AM MDT
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
enum StatExceptType { STATEXCEPT_UNIMPLEMENTED=0,  //!< function unimplemented
                      NUM_STAT_EXCEPT
};

//! StatException error messages
const string StatMsg[NUM_STAT_EXCEPT] = {
  "Method unimplemented"
};

//! StatException : Exception
class StatException : Exception {

  public:
    //! Constructor with filename
    StatException(ExceptType except, StatExceptType statExcept) :
      Exception(except), m_except(statExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    StatException(StatException&&) = default;

    //! Destructor
    ~StatException() = default;

    //! Handle StatException
    ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy constructor
    StatException(const StatException&) = delete;
    //! Don't permit copy assignment
    StatException& operator=(const StatException&) = delete;
    //! Don't permit move assignment
    StatException& operator=(StatException&&) = delete;

    //! Statistrics exception type (UNIMPLEMENTED, etc.)
    StatExceptType m_except;
};

} // namespace Quinoa

#endif // StatException_h
