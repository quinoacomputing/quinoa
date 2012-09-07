//******************************************************************************
/*!
  \file      src/Mesh/Mesh.h
  \author    J. Bakosi
  \date      Thu 06 Sep 2012 09:28:56 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh base class declaration
  \details   Mesh base class declaration
*/
//******************************************************************************
#ifndef Mesh_h
#define Mesh_h

namespace Quinoa {

//! Mesh base class
class Mesh {

  protected:
    //! Constructor
    Mesh();

    //! Destructor
    ~Mesh();

  private:
    //! Don't permit copy operator
    Mesh(const Mesh&);
    //! Don't permit assigment operator
    Mesh& operator=(const Mesh&);
};

} // namespace Quinoa

#endif // Mesh_h
