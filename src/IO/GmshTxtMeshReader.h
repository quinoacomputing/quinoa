//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.h
  \author    J. Bakosi
  \date      Tue May  7 12:52:33 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshTxtMeshReader_h
#define GmshTxtMeshReader_h

#include <map>

using namespace std;

#include <MeshReader.h>

namespace Quinoa {

//! GmshTxtMeshReader : MeshReader
class GmshTxtMeshReader : public MeshReader {

  public:
    //! Constructor
    explicit GmshTxtMeshReader(const string filename,
                               UnsMesh* const mesh,
                               Memory* const memory);

    //! Destructor, default compiler generated
    virtual ~GmshTxtMeshReader() noexcept = default;

    //! Read Gmsh mesh
    void read();

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
    void addElem(int type, vector<int>& nodes);

    //! Add new element tags
    void addElemTags(int type, vector<int>& tags);

    map<int,int> m_GmshElemNodes; //!< Element types and their number of nodes

    int m_nnodes;                  //!< Number of nodes
    int m_nLins;                   //!< Number of line elements
    int m_nTris;                   //!< Number of triangle elements

    int m_nodeCnt;                 //!< Counter for nodes added
    int m_linCnt;                  //!< Counter for line elems added
    int m_triCnt;                  //!< Counter for triangle elems added
};

} // namespace Quinoa

#endif // GmshTxtMeshReader_h
