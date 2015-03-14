//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.h
  \author    J. Bakosi
  \date      Sat 14 Mar 2015 07:00:35 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     STL (STereoLithography) mesh class declaration
  \details   STL (STereoLithography) mesh class declaration.
*/
//******************************************************************************
#ifndef STLMesh_h
#define STLMesh_h

#include <memory>
#include <string>

#include <Types.h>

namespace tk {

//! STLMesh
class STLMesh {

  public:
    //! Constructor
    explicit STLMesh() : m_nnode( 0 ) {}

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

    //! Node list accessor
    int* nodelist() const noexcept { return m_nodelist.get(); }

    //! Number of nodes accessor
    std::size_t nnode() const noexcept { return m_nnode; }

  private:
    std::string m_name;                         //!< Name of the mesh
    std::unique_ptr< tk::real[] > m_x;          //!< Vertex x coordinates
    std::unique_ptr< tk::real[] > m_y;          //!< Vertex y coordinates
    std::unique_ptr< tk::real[] > m_z;          //!< Vertex z coordinates
    std::unique_ptr< int[] > m_nodelist;        //!< Node indices for facets

    std::size_t m_nnode;                        //!< Number of nodes
};

} // tk::

#endif // STLMesh_h
