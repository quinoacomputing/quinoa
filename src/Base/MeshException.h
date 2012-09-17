//******************************************************************************
/*!
  \file      src/Base/MeshException.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 08:24:06 PM MDT
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
enum class MeshExceptType { BAD_FORMAT=0,   //!< unsupported Gmsh mesh format
                            BAD_ELEMENT,    //!< unknown element type
                            EMPTY_SET,      //!< no element/node sets
                            UNIMPLEMENTED,  //!< mesh file section unimplemented
                            NUM_MESH_EXCEPT
};
//! Number of mesh exception types
const Int NUM_MESH_EXCEPT = static_cast<Int>(MeshExceptType::NUM_MESH_EXCEPT);

//! Mesh exception error messages
const string MeshMsg[NUM_MESH_EXCEPT] = {
  "Unsupported mesh format: ",
  "Unknown element type in mesh file ",
  "No element/node sets in mesh",
  "Section not yet implemented: ",
};

//! MeshException : Exception
class MeshException : Exception {

  public:
    //! Constructor with message from thrower
    MeshException(ExceptType except,
                  MeshExceptType meshExcept,
                  string throwerMsg) : Exception(except),
      m_throwerMsg(throwerMsg), m_except(meshExcept) {}

    //! Constructor without message from thrower
    MeshException(ExceptType except, MeshExceptType meshExcept) :
      MeshException(except, meshExcept, 0) {}

    //! Move constructor, necessary for throws, default compiler generated
    MeshException(MeshException&&) = default;

    //! Destructor
    ~MeshException() = default;

    //! Handle MeshException
    ErrCode handleException(Driver* driver);

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
