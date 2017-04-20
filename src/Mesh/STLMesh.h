// *****************************************************************************
/*!
  \file      src/Mesh/STLMesh.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     STL (STereoLithography) mesh class declaration
  \details   STL (STereoLithography) mesh class declaration.
*/
// *****************************************************************************
#ifndef STLMesh_h
#define STLMesh_h

#include <memory>
#include <string>
#include <cstddef>
#include <iosfwd>

#include "Types.h"

namespace tk {

//! STLMesh
class STLMesh {

  public:
    //! Constructor
    explicit STLMesh() : m_name(), m_x(), m_y(), m_z(), m_nnode( 0 ) {}

    //! Allocate memory for mesh
    void alloc( std::size_t num );

    //! Set mesh name
    void setName( const std::string& n ) { m_name = n; }
    //! Get mesh name
    const std::string& name() const noexcept { return m_name; }

    //! Coordinate array accessors
    tk::real* getx() const noexcept { return m_x.get(); }
    tk::real* gety() const noexcept { return m_y.get(); }
    tk::real* getz() const noexcept { return m_z.get(); }

    //! Number of nodes accessor
    std::size_t nnode() const noexcept { return m_nnode; }

  private:
    std::string m_name;                         //!< Name of the mesh
    std::unique_ptr< tk::real[] > m_x;          //!< Vertex x coordinates
    std::unique_ptr< tk::real[] > m_y;          //!< Vertex y coordinates
    std::unique_ptr< tk::real[] > m_z;          //!< Vertex z coordinates

    std::size_t m_nnode;                        //!< Number of nodes
};

} // tk::

#endif // STLMesh_h
