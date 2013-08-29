//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:22:15 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     STL (STereoLithography) mesh class declaration
  \details   STL (STereoLithography) mesh class declaration
*/
//******************************************************************************
#ifndef STLMesh_h
#define STLMesh_h

#include <QuinoaTypes.h>
#include <Memory.h>

namespace quinoa {

//! STLMesh
class STLMesh {

  public:
    //! Constructor
    explicit STLMesh(Memory* const memory);

    //! Destructor
    virtual ~STLMesh() noexcept;

    //! Allocate memory for mesh
    void alloc(const size_t num);

    //! Set/get mesh name
    void setName(const std::string& name) { m_name = name; }
    const std::string& name() const noexcept { return m_name; }

    //! Coordinate array accessors
    real* getx() const noexcept { return m_x.ptr; }
    real* gety() const noexcept { return m_y.ptr; }
    real* getz() const noexcept { return m_z.ptr; }

    //! Node list accessor
    int* nodelist() const noexcept { return m_nodelist.ptr; }

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

    Memory* const m_memory;                  //!< Memory object pointer

    std::string m_name;                      //!< Name of the mesh
    Data<real> m_x;                          //!< Vertex x coordinates
    Data<real> m_y;                          //!< Vertex y coordinates
    Data<real> m_z;                          //!< Vertex z coordinates
    Data<int> m_nodelist;                    //!< Node indices describing facets

    size_t m_nnodes;                         //!< Number of nodes
};

} // namespace quinoa

#endif // STLMesh_h
