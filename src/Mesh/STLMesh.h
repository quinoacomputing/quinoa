//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:01:15 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     STL (STereoLithography) mesh class declaration
  \details   STL (STereoLithography) mesh class declaration
*/
//******************************************************************************
#ifndef STLMesh_h
#define STLMesh_h

#include <memory>
#include <string>

#include <Types.h>

namespace quinoa {

//! STLMesh
class STLMesh {

  public:
    //! Constructor
    explicit STLMesh() : m_nnodes(0) {}

    //! Destructor
    virtual ~STLMesh() = default;

    //! Allocate memory for mesh
    void alloc(const size_t num);

    //! Set/get mesh name
    void setName(const std::string& name) { m_name = name; }
    const std::string& name() const noexcept { return m_name; }

    //! Coordinate array accessors
    tk::real* getx() const noexcept { return m_x.get(); }
    tk::real* gety() const noexcept { return m_y.get(); }
    tk::real* getz() const noexcept { return m_z.get(); }

    //! Node list accessor
    int* nodelist() const noexcept { return m_nodelist.get(); }

    //! Number of nodes accessor
    size_t nnodes() const noexcept { return m_nnodes; }

  private:
    //! Don't permit copy constructor
    STLMesh(const STLMesh&) = delete;
    //! Don't permit assigment constructor
    STLMesh& operator=(const STLMesh&) = delete;
    //! Don't permit move constructor
    STLMesh(STLMesh&&) = delete;
    //! Don't permit move assignment
    STLMesh& operator=(STLMesh&&) = delete;

    std::string m_name;                      //!< Name of the mesh
    std::unique_ptr<tk::real[]> m_x;         //!< Vertex x coordinates
    std::unique_ptr<tk::real[]> m_y;         //!< Vertex y coordinates
    std::unique_ptr<tk::real[]> m_z;         //!< Vertex z coordinates
    std::unique_ptr<int[]> m_nodelist;       //!< Node indices describing facets

    size_t m_nnodes;                         //!< Number of nodes
};

} // quinoa::

#endif // STLMesh_h
