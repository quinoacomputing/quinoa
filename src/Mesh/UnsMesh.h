//******************************************************************************
/*!
  \file      src/Base/UnsMesh.h
  \author    J. Bakosi
  \date      Tue 11 Sep 2012 12:20:41 PM KST
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
class UnsMesh : public Mesh {

  public:
    //! Constructor
    UnsMesh(Memory* memory) : Mesh(memory) {};

    //! Destructor
    ~UnsMesh() {};

  private:
    //! Don't permit copy operator
    UnsMesh(const UnsMesh&);

    //! Dont' permit assigment operator
    UnsMesh& operator=(const UnsMesh&);
};

} // namespace Quinoa

#endif // UnsMesh_h
