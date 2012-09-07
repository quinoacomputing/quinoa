//******************************************************************************
/*!
  \file      src/Base/StrMesh.h
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 12:32:50 PM MDT
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
    //! Constructor
    StrMesh();

    //! Destructor
    ~StrMesh();

  private:
    //! Don't permit copy operator
    StrMesh(const StrMesh&);

    //! Dont' permit assigment operator
    StrMesh& operator=(const StrMesh&);
};

} // namespace Quinoa

#endif // StrMesh_h
