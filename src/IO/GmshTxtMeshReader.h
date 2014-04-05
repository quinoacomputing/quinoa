//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.h
  \author    J. Bakosi
  \date      Sat 05 Apr 2014 12:03:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshTxtMeshReader_h
#define GmshTxtMeshReader_h

#include <vector>
#include <map>

#include <Reader.h>
#include <GmshMesh.h>

namespace quinoa {

//! GmshTxtMeshReader : Reader
class GmshTxtMeshReader : public tk::Reader {

  public:
    //! Constructor
    explicit GmshTxtMeshReader(const std::string filename, GmshMesh& mesh);

    //! Destructor, default compiler generated
    ~GmshTxtMeshReader() noexcept override = default;

    //! Read Gmsh mesh
    void read() override;

  private:
    //! Don't permit copy constructor
    GmshTxtMeshReader(const GmshTxtMeshReader&) = delete;
    //! Don't permit copy assigment
    GmshTxtMeshReader& operator=(const GmshTxtMeshReader&) = delete;
    //! Don't permit move constructor
    GmshTxtMeshReader(GmshTxtMeshReader&&) = delete;
    //! Don't permit move assigment
    GmshTxtMeshReader& operator=(GmshTxtMeshReader&&) = delete;

    //! Read mandatory "$MeshFormat--$EndMeshFormat" section
    void readMeshFormat();

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();

    //! Add new element
    void addElem(int type, const std::vector< int >& nodes);

    //! Add new element tags
    void addElemTags(int type, const std::vector< int >& tags);

    GmshMesh& m_mesh;              //!< Mesh object

    //! Element types and their number of nodes
    std::map<int,int> m_GmshElemNodes;
};

} // quinoa::

#endif // GmshTxtMeshReader_h
