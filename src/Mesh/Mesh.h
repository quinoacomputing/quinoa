//******************************************************************************
/*!
  \file      src/Mesh/Mesh.h
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:12:50 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh base class declaration
  \details   Mesh base class declaration
*/
//******************************************************************************
#ifndef Mesh_h
#define Mesh_h

#include <QuinoaTypes.h>

namespace Quinoa {

//! Mesh dimension
enum class MeshDim { TWOD=0,
                     THREED,
                     NUM_MESH_DIM
};
//! Number of mesh dimensions
const int NUM_MESH_DIM = static_cast<int>(MeshDim::NUM_MESH_DIM);

//! Mesh base class
class Mesh {

  public:
    //! Constructor, default compiler generated
    // ICC wants it public
    Mesh() = default;

  protected:
    //! Destructor, default compiler generated
    ~Mesh() = default;

    //! Set mesh dimension
    void setDim(MeshDim dim);

  private:
    //! Don't permit copy constructor
    //Mesh(const Mesh&) = delete;
    //! Don't permit assigment constructor
    Mesh& operator=(const Mesh&) = delete;
    //! Don't permit move constructor
    Mesh(Mesh&&) = delete;
    //! Don't permit move assignment
    Mesh& operator=(Mesh&&) = delete;

    //! Mesh dimension
    MeshDim m_dim;
};

} // namespace Quinoa

#endif // Mesh_h
