//******************************************************************************
/*!
  \file      src/Base/IOException.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 07:11:46 PM MDT
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
enum class IOExceptType { FAILED_OPEN=0,  //!< failed to open file
                          FAILED_CLOSE,   //!< failed to close file
                          NUM_IO_EXCEPT
};
//! Number of IO exception types
const Int NUM_IO_EXCEPT = static_cast<Int>(IOExceptType::NUM_IO_EXCEPT);

//! IOException error messages
const string IOMsg[NUM_IO_EXCEPT] = {
  "Failed to open file: ",
  "Failed to close file: "
};

//! IOException : Exception
class IOException : Exception {

  public:
    //! Constructor with filename
    IOException(ExceptType except, IOExceptType ioExcept, string filename) :
      Exception(except), m_filename(filename), m_except(ioExcept) {}

    //! Constructor without filename
    IOException(ExceptType except, IOExceptType ioExcept) :
      IOException(except, ioExcept, 0) {}

    //! Move constructor, necessary for throws, default compiler generated
    IOException(IOException&&) = default;

    //! Destructor
    ~IOException() = default;

    //! Handle IOException
    ErrCode handleException(Driver* driver);

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

    //! IO exception type (FAILED_OPEN, FAILED_CLOSE, etc.)
    IOExceptType m_except;
};

} // namespace Quinoa

#endif // IOException_h
