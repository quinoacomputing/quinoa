//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.h
  \author    J. Bakosi
  \date      Thu 03 Oct 2013 08:32:57 PM MDT
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

namespace quinoa {

//! GmshTxtMeshReader : Reader
class GmshTxtMeshReader : public Reader {

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

    //! Count up elements, nodes, physicals
    void count();

    //! Read "$Nodes--$EndNodes" section and count nodes
    void countNodes();

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section and count elements
    void countElements();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section and count physicals
    void countPhysicalNames();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();

    //! Add new element
    void addElem(int type, std::vector<int>& nodes);

    //! Add new element tags
    void addElemTags(int type, std::vector<int>& tags);

    GmshMesh& m_mesh;              //!< Mesh object

    //! Element types and their number of nodes
    std::map<int,int> m_GmshElemNodes;

    int m_nnodes;                  //!< Number of nodes
    int m_nLins;                   //!< Number of line elements
    int m_nTris;                   //!< Number of triangle elements

    int m_nodeCnt;                 //!< Counter for nodes added
    int m_linCnt;                  //!< Counter for line elems added
    int m_triCnt;                  //!< Counter for triangle elems added
};

} // namespace quinoa

#endif // GmshTxtMeshReader_h
