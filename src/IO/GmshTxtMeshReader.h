//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.h
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 09:13:18 PM MDT
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
class GmshTxtMeshReader : MeshReader {

  public:
    //! Constructor
    GmshTxtMeshReader(string filename, UnsMesh* mesh, Memory* memory);

    //! Destructor, default compiler generated
    ~GmshTxtMeshReader() = default;

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
    void addElem(Int type, vector<Int>& nodes);

    //! Add new element tags
    void addElemTags(Int type, vector<Int>& tags);

    map<Int, Int> m_GmshElemNodes; //!< Element types and their number of nodes

    Int m_nnodes;                  //!< Number of nodes
    Int m_nLins;                   //!< Number of line elements
    Int m_nTris;                   //!< Number of triangle elements

    Int m_nodeCnt;                 //!< Counter for nodes added
    Int m_linCnt;                  //!< Counter for line elems added
    Int m_triCnt;                  //!< Counter for triangle elems added
};

} // namespace Quinoa

#endif // GmshTxtMeshReader_h
