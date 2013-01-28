//******************************************************************************
/*!
  \file      src/Mesh/MeshException.h
  \author    J. Bakosi
  \date      Sun 27 Jan 2013 07:55:46 PM MST
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
enum MeshExceptType { BAD_FORMAT=0,   //!< unsupported Gmsh mesh format
                      BAD_ELEMENT,    //!< unknown element type
                      EMPTY_SET,      //!< no element/node sets
                                      //! mesh file section unimplemented
                      MESHEXCEPT_UNIMPLEMENTED,
                      NUM_MESH_EXCEPT
};

//! Mesh exception error messages
const string MeshMsg[NUM_MESH_EXCEPT] = {
  "Unsupported mesh format: ",
  "Unknown element type in mesh file ",
  "No element/node sets in mesh",
  "Section not yet implemented: "
};

//! MeshException : Exception
class MeshException : public Exception {

  public:
    //! Constructor with message from thrower
    MeshException(ExceptType except,
                  MeshExceptType mshExcept,
                  const string& file,
                  const string& func,
                  const unsigned int& line) :
      Exception(except, file, func, line), m_except(mshExcept) {}

    //! Constructor without message
    MeshException(ExceptType except,
                  MeshExceptType mshExcept,
                  const string throwerMsg,
                  const string& file,
                  const string& func,
                  const unsigned int& line) :
      Exception(except, file, func, line), m_throwerMsg(throwerMsg),
      m_except(mshExcept) {}

    // ICC: the above can be done with delegate constructors, like
    //MeshException(ExceptType except, MeshExceptType meshExcept) :
    //  MeshException(except, meshExcept, 0) {}

    //! Move constructor, necessary for throws, default compiler generated
    MeshException(MeshException&&) = default;

    //! Destructor
    virtual ~MeshException() {}

    //! Handle MeshException
    virtual ErrCode handleException(Driver* driver);

  protected:
    //! Message from thrower
    string m_throwerMsg;

  private:
    //! Don't permit copy constructor
    MeshException(const MeshException&) = delete;
    //! Don't permit copy assignment
    MeshException& operator=(const MeshException&) = delete;
    //! Don't permit move assignment
    MeshException& operator=(MeshException&&) = delete;

    //! Mesh exception type (BAD_FORMAT, BAD_ELEMENT, etc.)
    MeshExceptType m_except;
};

} // namespace Quinoa

#endif // MeshException_h
