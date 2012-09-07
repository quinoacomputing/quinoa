//******************************************************************************
/*!
  \file      src/Base/TriMesh.h
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 12:35:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Structured triangle mesh class declaration
  \details   Structured triangle mesh class declaration
*/
//******************************************************************************
#ifndef TriMesh_h
#define TriMesh_h

#include <UnsMesh.h>

namespace Quinoa {

//! TriMesh : UnsMesh
class TriMesh : UnsMesh {

  public:
    //! Constructor
    TriMesh();

    //! Destructor
    ~TriMesh();

  private:
    //! Don't permit copy operator
    TriMesh(const TriMesh&);

    //! Dont' permit assigment operator
    TriMesh& operator=(const TriMesh&);
};

} // namespace Quinoa

#endif // TriMesh_h
