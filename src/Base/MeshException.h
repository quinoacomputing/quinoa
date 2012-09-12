//******************************************************************************
/*!
  \file      src/Base/MeshException.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 04:35:04 AM KST
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
enum MeshExceptionType { BAD_FORMAT=0,   //!< unsupported Gmsh mesh format
                         BAD_ELEMENT,    //!< unknown element type
                         EMPTY_SET,      //!< no element/node sets
                         NUM_MESH_EXCEPTIONS
};

//! MeshException error messages
const string MeshMessage[NUM_MESH_EXCEPTIONS] = {
  "Unsupported mesh format: ",
  "Unknown element type in mesh file ",
  "No element/node sets in mesh",
};

//! MeshException : Exception
class MeshException : Exception {

  public:
    //! Constructor with filename
    MeshException(ExceptionType exception,
                  MeshExceptionType meshException,
                  string filename) :
      Exception(exception), m_filename(filename), m_exception(meshException) {}

    //! Constructor without filename
    MeshException(ExceptionType exception,
                  MeshExceptionType meshException) :
      Exception(exception), m_exception(meshException) {}

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

    //! MeshException (BAD_FORMAT, BAD_ELEMENT, etc.)
    MeshExceptionType m_exception;
};

} // namespace Quinoa

#endif // MeshException_h
