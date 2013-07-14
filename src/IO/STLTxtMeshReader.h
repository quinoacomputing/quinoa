//******************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.h
  \author    J. Bakosi
  \date      Sat 13 Jul 2013 10:32:44 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) reader class declaration
  \details   ASCII STL (STereoLithographu) reader class declaration
*/
//******************************************************************************
#ifndef STLTxtMeshReader_h
#define STLTxtMeshReader_h

#include <vector>

#include <Reader.h>

class STLMesh;

namespace Quinoa {

//! STLTxtMeshReader : Reader
class STLTxtMeshReader : public Reader {

  public:
    //! Constructor
    explicit STLTxtMeshReader(const std::string filename,
                              STLMesh* const mesh) :
      Reader(filename),
      m_mesh(mesh),
      m_name(),
      m_Triangles() {
      Assert(m_mesh != nullptr, ExceptType::FATAL,
            "Uninitialized mesh object passed to STLTxtMeshReader constructor");
    }

    //! Destructor, default compiler generated
    virtual ~STLTxtMeshReader() noexcept = default;

    //! Read ASCII STL mesh
    void read();

  private:
    //! Don't permit copy constructor
    STLTxtMeshReader(const STLTxtMeshReader&) = delete;
    //! Don't permit copy assigment
    STLTxtMeshReader& operator=(const STLTxtMeshReader&) = delete;
    //! Don't permit move constructor
    STLTxtMeshReader(STLTxtMeshReader&&) = delete;
    //! Don't permit move assigment
    STLTxtMeshReader& operator=(STLTxtMeshReader&&) = delete;

    //! Vertex
    struct Vertex {
      real x;
      real y;
      real z;
    };

    // Triangle: 3 vertices: A, B, C
    struct Triangle {
      Vertex A;
      Vertex B;
      Vertex C;
    };

    STLMesh* const m_mesh;                      //!< Mesh object pointer

    std::string m_name;                         //!< Model name
    std::vector<Triangle> m_Triangles;          //!< Vector of triangles
};

} // namespace Quinoa

#endif // STLTxtMeshReader_h
