//******************************************************************************
/*!
  \file      src/Base/UnsMesh.h
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 04:20:14 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class declaration
  \details   Unstructured mesh class declaration
*/
//******************************************************************************
#ifndef UnsMesh_h
#define UnsMesh_h

#include <Mesh.h>

namespace Quinoa {

//! UnsMesh : Mesh
class UnsMesh : Mesh {

  public:
    //! Constructor
    UnsMesh();

    //! Destructor
    virtual ~UnsMesh();

  private:
    //! Don't permit copy operator
    UnsMesh(const UnsMesh&);

    //! Dont' permit assigment operator
    UnsMesh& operator=(const UnsMesh&);
};

} // namespace Quinoa

#endif // UnsMesh_h
