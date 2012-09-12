//******************************************************************************
/*!
  \file      src/Base/IOException.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 04:32:04 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     IOException class declaration
  \details   IOException class declaration
*/
//******************************************************************************
#ifndef IOException_h
#define IOException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! IOException types
enum IOExceptionType { FAILED_OPEN=0,  //!< failed to open file
                       FAILED_CLOSE,   //!< failed to close file
                       NUM_IO_EXCEPTIONS
};

//! IOException error messages
const string IOMessage[NUM_IO_EXCEPTIONS] = {
  "Failed to open file: ",
  "Failed to close file: "
};

//! IOException : Exception
class IOException : Exception {

  public:
    //! Constructor with filename
    IOException(ExceptionType exception,
                IOExceptionType ioException,
                string filename) :
      Exception(exception), m_filename(filename), m_exception(ioException) {}

    //! Constructor without filename
    IOException(ExceptionType exception,
                  IOExceptionType ioException) :
      Exception(exception), m_exception(ioException) {}

    //! Copy constructor
    IOException(const IOException&);

    //! Destructor
    ~IOException() {}

    //! Handle IOException
    ErrorCode handleException(Driver* driver);

  protected:
    //! File name
    string m_filename;

  private:
    //! Dont' permit assigment operator
    IOException& operator=(const IOException&);

    //! IOException (FAILED_OPEN, FAILED_CLOSE, etc.)
    IOExceptionType m_exception;
};

} // namespace Quinoa

#endif // IOException_h
