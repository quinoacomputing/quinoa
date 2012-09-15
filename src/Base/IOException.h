//******************************************************************************
/*!
  \file      src/Base/IOException.h
  \author    J. Bakosi
  \date      Fri Sep 14 17:29:23 2012
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
      Exception(except), m_except(ioExcept) {}

    //! Copy constructor
    IOException(const IOException&);

    //! Destructor
    ~IOException() {}

    //! Handle IOException
    ErrCode handleException(Driver* driver);

  protected:
    //! File name
    string m_filename;

  private:
    //! Dont' permit assigment operator
    IOException& operator=(const IOException&);

    //! IO exception type (FAILED_OPEN, FAILED_CLOSE, etc.)
    IOExceptType m_except;
};

} // namespace Quinoa

#endif // IOException_h
