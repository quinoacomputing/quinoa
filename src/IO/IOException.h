//******************************************************************************
/*!
  \file      src/IO/IOException.h
  \author    J. Bakosi
  \date      Fri Apr 26 14:48:34 2013
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
    //! Constructor
    explicit IOException(const ExceptType except,
                         const IOExceptType ioExcept,
                         const string filename,
                         const string& file,
                         const string& func,
                         const unsigned int& line) :
      Exception(except, file, func, line), m_filename(filename),
      m_except(ioExcept) {}

    //! Move constructor for throws, default compiler generated
    IOException(IOException&&) = default;

    //! Destructor
    virtual ~IOException() noexcept = default;

    //! Handle IOException
    virtual ErrCode handleException(Driver* const driver);

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
    const IOExceptType m_except;
};

} // namespace Quinoa

#endif // IOException_h
