//******************************************************************************
/*!
  \file      src/Mesh/MeshException.h
  \author    J. Bakosi
  \date      Mon Apr 29 15:51:11 2013
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
    explicit MeshException(const ExceptType except,
                           const MeshExceptType mshExcept,
                           const string throwerMsg,
                           const string& file,
                           const string& func,
                           unsigned int line) noexcept :
      Exception(except,
                MeshMsg[static_cast<int>(mshExcept)] + throwerMsg,
                file,
                func,
                line) {}

    //! Move constructor for throws, default compiler generated
    MeshException(MeshException&&) = default;

    //! Destructor
    virtual ~MeshException() noexcept = default;

  private:
    //! Don't permit copy constructor
    MeshException(const MeshException&) = delete;
    //! Don't permit copy assignment
    MeshException& operator=(const MeshException&) = delete;
    //! Don't permit move assignment
    MeshException& operator=(MeshException&&) = delete;
};

} // namespace Quinoa

#endif // MeshException_h
