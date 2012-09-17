//******************************************************************************
/*!
  \file      src/Base/StrMesh.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 06:25:25 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Structured mesh class declaration
  \details   Structured mesh class declaration
*/
//******************************************************************************
#ifndef StrMesh_h
#define StrMesh_h

#include <Mesh.h>

namespace Quinoa {

//! StrMesh : Mesh
class StrMesh : Mesh {

  public:
    //! Constructor, default compiler generated
    StrMesh() = default;

    //! Destructor, default compiler generated
    ~StrMesh() = default;

  private:
    //! Don't permit copy constructor
    StrMesh(const StrMesh&) = delete;
    //! Don't permit assigment constructor
    StrMesh& operator=(const StrMesh&) = delete;
    //! Don't permit move constructor
    StrMesh(StrMesh&&) = delete;
    //! Don't permit move assignment
    StrMesh& operator=(StrMesh&&) = delete;
};

} // namespace Quinoa

#endif // StrMesh_h
