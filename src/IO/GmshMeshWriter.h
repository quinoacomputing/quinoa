//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.h
  \author    J. Bakosi
  \date      Wed 09 Apr 2014 07:22:09 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshMeshWriter class declaration
  \details   GmshMeshWriter class declaration
*/
//******************************************************************************
#ifndef GmshMeshWriter_h
#define GmshMeshWriter_h

#include <string>
#include <iostream>     // NOT NEEDED

#include <Writer.h>
#include <GmshMesh.h>

namespace quinoa {

//! GmshMeshWriter : Writer
class GmshMeshWriter : public tk::Writer {

  public:
    //! Constructor
    explicit GmshMeshWriter( const std::string& filename,
                             GmshMesh& mesh,
                             tk::real version = 2.2,
                             GmshMesh::FileType type =
                               GmshMesh::FileType::BINARY, // default: binary
                             int datasize = sizeof(double) );

    //! Destructor, default compiler generated
    ~GmshMeshWriter() noexcept override = default;

    //! Write Gmsh mesh to file
    void write();

  private:
    //! Don't permit copy constructor
    GmshMeshWriter(const GmshMeshWriter&) = delete;
    //! Don't permit copy assigment
    GmshMeshWriter& operator=(const GmshMeshWriter&) = delete;
    //! Don't permit move constructor
    GmshMeshWriter(GmshMeshWriter&&) = delete;
    //! Don't permit move assigment
    GmshMeshWriter& operator=(GmshMeshWriter&&) = delete;

    //! Write "$Nodes--$EndNodes" section
    void writeNodes();

    //! Write "$Elements--$EndElements" section
    void writeElements();

    //! Write "$PhysicalNames--$EndPhysicalNames" section
    void writePhysicalNames();

    // Get mesh type
    bool isASCII() const {
      return m_type == GmshMesh::FileType::ASCII ? true : false;
    }
    bool isBinary() const {
      return m_type == GmshMesh::FileType::BINARY ? true : false;
    }

    // Write out element ids, tags, and connectivity (node list)
    template< class ElmId, class ElmTag, class ElmInpoel >
    void writeElemBlock( GmshMesh::ElmType type, ElmId& id, ElmTag& tag,
                         ElmInpoel& inpoel )
    {
      if (tag.size() == 0 || id.size() == 0 || inpoel.size() == 0) return;

      auto n = inpoel.size();
      if (isASCII()) {
        for (std::size_t i=0; i<n; i++) {
          // elm-number elm-type number-of-tags < tag > ... node-number-list
          m_outFile << id[i] << " " << type << " " << tag[i].size() << " ";

          copy( tag[i].begin(), tag[i].end()-1,
                std::ostream_iterator< int >( m_outFile, " " ) );
          m_outFile << tag[i].back() << " ";

          copy( inpoel[i].begin(), inpoel[i].end()-1,
                std::ostream_iterator< int >( m_outFile, " ") );
          m_outFile << inpoel[i].back() << std::endl;
        }
      } else {
        int ntags = tag[0].size();
        // elm-type num-of-elm-follow number-of-tags
        m_outFile.write( reinterpret_cast<char*>(&type), sizeof(int) );
        m_outFile.write( reinterpret_cast<char*>(&n), sizeof(int) );
        m_outFile.write( reinterpret_cast<char*>(&ntags), sizeof(int) );
        for (std::size_t i=0; i<n; i++) {
          // element id
          m_outFile.write(
            reinterpret_cast<char*>(&id[i]), sizeof(int) );
          // element tags
          m_outFile.write( reinterpret_cast<char*>(tag[i].data()),
                           tag[i].size()*sizeof(int) );
          // element node list (i.e. connectivity)
          m_outFile.write( reinterpret_cast<char*>(inpoel[i].data()),
                           inpoel[i].size()*sizeof(int) );
        }
      }
    }

    GmshMesh& m_mesh;                   //!< Mesh object
    GmshMesh::FileType m_type;          //!< Mesh file type: 0:ASCII, 1:binary
};

} // quinoa::

#endif // GmshMeshWriter_h
