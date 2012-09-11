//******************************************************************************
/*!
  \file      src/Base/MeshException.h
  \author    J. Bakosi
  \date      Tue 11 Sep 2012 06:53:24 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshException class declaration
  \details   MeshException class declaration
*/
//******************************************************************************
#ifndef MeshException_h
#define MeshException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! MeshException types
enum MeshExceptionType { FAILED_OPEN=0,  //!< failed to open file
                         FAILED_CLOSE,   //!< failed to close file
                         BAD_FORMAT,     //!< unsupported Gmsh mesh format
                         BAD_ELEMENT,    //!< unknown element type
                         NUM_MESH_EXCEPTIONS
};

//! MeshException error messages
const string MeshMessage[NUM_MESH_EXCEPTIONS] = {
  "Failed to open mesh file: ",
  "Failed to close mesh file: ",
  "Unsupported mesh format: ",
  "Unknown element type in mesh file "
};

//! MeshException : Exception class declaration
class MeshException : Exception {

  public:
    //! Constructor
    MeshException(ExceptionType exception,
                  MeshExceptionType meshException,
                  string filename) :
      Exception(exception), m_filename(filename), m_exception(meshException) {}

    //! Copy constructor
    MeshException(const MeshException&);

    //! Destructor
    ~MeshException() {}

    //! Handle MeshException
    ErrorCode handleException(Driver* driver);

  protected:
    //! Mesh file name
    string m_filename;

  private:
    //! Dont' permit assigment operator
    MeshException& operator=(const MeshException&);

    //! MeshException (BAD_INPUTFILE, BAD_FORMAT, etc.)
    MeshExceptionType m_exception;
};

} // namespace Quinoa

#endif // MeshException_h
