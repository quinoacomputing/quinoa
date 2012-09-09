//******************************************************************************
/*!
  \file      src/Base/GmshException.h
  \author    J. Bakosi
  \date      Mon 10 Sep 2012 04:10:12 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshException class declaration
  \details   GmshException class declaration
*/
//******************************************************************************
#ifndef GmshException_h
#define GmshException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! GmshException types
enum GmshExceptionType { FAILED_OPEN=0,  //!< failed to open file
                         FAILED_CLOSE,   //!< failed to close file
                         BAD_FORMAT,     //!< unsupported Gmsh mesh format
                         BAD_NUMBERING,  //!< point numbering starts with 0
                         BAD_ELEMENT,    //!< unknown element type
                         NUM_GMSH_EXCEPTIONS
};

//! GmshException error messages
const string GmshMessage[NUM_GMSH_EXCEPTIONS] = {
  "Failed to open mesh file: ",
  "Failed to close mesh file: ",
  "Unsupported mesh format: ",
  "Point numbering starts with 0 in ",
  "Unknown element type in mesh file "
};

//! GmshException : Exception class declaration
class GmshException : Exception {

  public:
    //! Constructor
    GmshException(ExceptionType exception,
                  GmshExceptionType gmshException,
                  string filename) :
      Exception(exception), m_exception(gmshException), m_filename(filename) {}

    //! Copy constructor
    GmshException(const GmshException&);

    //! Destructor
    ~GmshException() {}

    //! Handle GmshException
    ErrorCode handleException(Driver* driver);

  private:
    //! Dont' permit assigment operator
    GmshException& operator=(const GmshException&);

    //! GmshException (BAD_INPUTFILE, BAD_FORMAT, etc.)
    GmshExceptionType m_exception;

    //! Mesh file name
    string m_filename;
};

} // namespace Quinoa

#endif // GmshException_h
