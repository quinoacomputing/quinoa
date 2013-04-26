//******************************************************************************
/*!
  \file      src/Mesh/Mesh.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:48:26 2013
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
enum class MeshDim { TWO=0,
                     THREE,
                     NUM_MESH_DIM
};
//! Number of mesh dimensions
const int NUM_MESH_DIM = static_cast<int>(MeshDim::NUM_MESH_DIM);

//! Mesh base class
class Mesh {

  protected:
    //! Constructor, default compiler generated
    explicit Mesh() : m_dim(MeshDim::THREE) {}

    //! Destructor, default compiler generated
    virtual ~Mesh() noexcept = default;

    //! Set mesh dimension
    void setDim(const MeshDim dim);

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
