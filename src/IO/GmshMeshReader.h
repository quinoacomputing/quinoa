//******************************************************************************
/*!
  \file      src/IO/GmshMeshReader.h
  \author    J. Bakosi
  \date      Thu 10 Apr 2014 09:35:49 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshMeshReader_h
#define GmshMeshReader_h

#include <map>

#include <Reader.h>
#include <GmshMesh.h>
#include <GmshMeshIO.h>
#include <Exception.h>

namespace quinoa {

//! GmshMeshReader : Reader
class GmshMeshReader : public tk::Reader {

  public:
    //! Constructor
    explicit GmshMeshReader( const std::string filename, GmshMesh& mesh ) :
      Reader(filename),
      m_version( 0.0 ),                        // 0.0: uninitialized
      m_datasize( 0 ),                         //   0: uninitialized
      m_type( GmshFileType::UNDEFINED ),       //  -1: uninitialized
      m_mesh( mesh ) {}

    //! Destructor, default compiler generated
    ~GmshMeshReader() noexcept override = default;

    //! Read Gmsh mesh
    void read() override;

  private:
    //! Don't permit copy constructor
    GmshMeshReader(const GmshMeshReader&) = delete;
    //! Don't permit copy assigment
    GmshMeshReader& operator=(const GmshMeshReader&) = delete;
    //! Don't permit move constructor
    GmshMeshReader(GmshMeshReader&&) = delete;
    //! Don't permit move assigment
    GmshMeshReader& operator=(GmshMeshReader&&) = delete;

    //! Read mandatory "$MeshFormat--$EndMeshFormat" section
    void readMeshFormat();

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();

    // Get mesh type queries
    bool isASCII() const {
      Assert( m_type != GmshFileType::UNDEFINED, tk::ExceptType::FATAL,
              "Mesh type is undefined");
      return m_type == GmshFileType::ASCII ? true : false;
    }
    bool isBinary() const {
      Assert( m_type != GmshFileType::UNDEFINED, tk::ExceptType::FATAL,
              "Mesh type is undefined");
      return m_type == GmshFileType::BINARY ? true : false;
    }

    tk::real m_version;                 //!< Mesh version in mesh file
    int m_datasize;                     //!< Data size in mesh file
    GmshFileType m_type;                //!< Mesh file type: 0:ASCII, 1:binary
    GmshMesh& m_mesh;                   //!< Mesh object

    //! Push back p for different element types
    template< class Pushed, class ElmType1, class ElmType2, class ElmType3 >
    void push_back( int etype, const Pushed& p, ElmType1& e1, ElmType2& e2,
                    ElmType3& e3 )
    {
      using tk::operator+;
      switch ( etype ) {
        case GmshElemType::LIN: e1.push_back( p ); break;
        case GmshElemType::TRI: e2.push_back( p ); break;
        case GmshElemType::TET: e3.push_back( p ); break;
        case GmshElemType::PNT: break;     // ignore 1-node 'point element' type
        default: Throw( tk::ExceptType::FATAL,
                        std::string("Unsupported element type ") + etype +
                        " in mesh file: " + m_filename );

      }
    }

    //! Gmsh element types and their number of nodes,
    //! all Gmsh-supported listed, Quinoa-supported uncommented,
    //! See Gmsh documentation for element ids as keys
    const std::map< int, int > m_elemNodes {
      { GmshElemType::LIN,   2 },  // 2-node line
      { GmshElemType::TRI,   3 },  // 3-node triangle
    //{             3,   4 },      // 4-node quadrangle
      { GmshElemType::TET,   4 },  // 4-node tetrahedron
    //{             5,   8 },      // 8-node hexahedron
    //{             6,   6 },      // 6-node prism
    //{             7,   5 },      // 5-node pyramid
    //{             8,   3 },      // 3-node second order line
    //{             9,   6 },      // 6-node second order triangle
    //{            10,   9 },      // 9-node second order quadrangle
    //{            11,  10 },      // 10-node second order tetrahedron
    //{            12,  27 },      // 27-node second order hexahedron
    //{            13,  18 },      // 18-node second order prism
    //{            14,  14 },      // 14-node second order pyramid
      { GmshElemType::PNT,   1 }   // 1-node point
    //{            16,   8 },      // 8-node second order quadrangle
    //{            17,  20 },      // 20-node second order hexahedron
    //{            18,  15 },      // 15-node second order prism
    //{            19,  13 },      // 13-node second order pyramid
    //{            20,   9 },      // 9-node third order incomplete triangle
    //{            21,  10 },      // 10-node third order triangle
    //{            22,  12 },      // 12-node fourth order incomplete triangle
    //{            23,  15 },      // 15-node fourth order triangle
    //{            24,  15 },      // 15-node fifth order incomplete triangle
    //{            25,  21 },      // 21-node fifth order complete triangle
    //{            26,   4 },      // 4-node third order edge
    //{            27,   5 },      // 5-node fourth order edge
    //{            28,   6 },      // 6-node fifth order edge
    //{            29,  20 },      // 20-node third order tetrahedron
    //{            30,  35 },      // 35-node fourth order tetrahedron
    //{            31,  56 },      // 56-node fifth order tetrahedron
    //{            92,  64 },      // 64-node third order hexahedron
    //{            93, 125 }       // 125-node fourth order hexahedron
    };
};

} // quinoa::

#endif // GmshMeshReader_h
