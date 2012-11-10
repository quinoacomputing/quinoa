//******************************************************************************
/*!
  \file      src/IO/IOException.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 06:18:58 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     IOException class declaration
  \details   IOException class declaration
*/
//******************************************************************************
#ifndef IOException_h
#define IOException_h

#include <string>

using namespace std;

#include <QuinoaTypes.h>
#include <Exception.h>

namespace Quinoa {

//! IOException types
enum IOExceptType { IO_FAILED_OPEN=0,  //!< failed to open file
                    IO_FAILED_CLOSE,   //!< failed to close file
                    IO_FAILED_WRITE,   //!< failed to write to file
                    NUM_IO_EXCEPT
};

//! IOException error messages
const string IOMsg[NUM_IO_EXCEPT] = {
  "Failed to open file: ",
  "Failed to close file: ",
  "Failed to write to file: "
};

//! IOException : Exception
class IOException : public Exception {

  public:
    //! Constructor with filename
    IOException(ExceptType except, IOExceptType ioExcept, string filename) :
      Exception(except), m_filename(filename), m_except(ioExcept) {}

    //! Constructor without filename
    // ICC: no delegate constructors
    //IOException(ExceptType except, IOExceptType ioExcept) :
    //  IOException(except, ioExcept, 0) {}

    //! Move constructor, necessary for throws, default compiler generated
    IOException(IOException&&) = default;

    //! Destructor
    virtual ~IOException() {}

    //! Handle IOException
    virtual ErrCode handleException(Driver* driver);

  protected:
    //! File name
    string m_filename;

  private:
    //! Don't permit copy constructor
    IOException(const IOException&) = delete;
    //! Don't permit copy assignment
    IOException& operator=(const IOException&) = delete;
    //! Don't permit move assignment
    IOException& operator=(IOException&&) = delete;

    //! IO exception type (IO_FAILED_OPEN, IO_FAILED_CLOSE, etc.)
    IOExceptType m_except;
};

} // namespace Quinoa

#endif // IOException_h
