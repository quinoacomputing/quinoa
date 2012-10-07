//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.h
  \author    J. Bakosi
  \date      Sun 07 Oct 2012 09:51:35 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshTxtMeshReader_h
#define GmshTxtMeshReader_h

#include <map>

using namespace std;

#include <UnsMesh.h>
#include <MeshReader.h>
#include <MeshException.h>

namespace Quinoa {

//! Base names for various mesh memory entries
const string     NODEID_NAME = "nodeId";
const string      COORD_NAME = "coord";
const string      LINID_NAME = "linId";
const string      TRIID_NAME = "triId";
const string   LINNODES_NAME = "linNode";
const string   TRINODES_NAME = "triNode";
const string     LINTAG_NAME = "linTag";
const string     TRITAG_NAME = "triTag";

//! Gmsh element types and their number of nodes,
//! all Gmsh-supported listed, Quinoa-supported at this time uncommented
const map<Int, Int> GmshElemNodes = {
          { 1,  2},  //! 2-node line
          { 2,  3},  //! 3-node triangle
//        { 3,  4},  //! 4-node quadrangle
//        { 4,  4},  //! 4-node tetrahedron
//        { 5,  8},  //! 8-node hexahedron
//        { 6,  6},  //! 6-node prism
//        { 7,  5},  //! 5-node pyramid
//        { 8,  3},  //! 3-node second order line
//        { 9,  6},  //! 6-node second order triangle
//        {10,  9},  //! 9-node second order quadrangle
//        {11, 10},  //! 10-node second order tetrahedron
//        {12, 27},  //! 27-node second order hexahedron
//        {13, 18},  //! 18-node second order prism
//        {14, 14},  //! 14-node second order pyramid
//        {15,  1},  //! 1-node point
//        {16,  8},  //! 8-node second order quadrangle
//        {17, 20},  //! 20-node second order hexahedron
//        {18, 15},  //! 15-node second order prism
//        {19, 13},  //! 13-node second order pyramid
//        {20,  9},  //! 9-node third order incomplete triangle
//        {21, 10},  //! 10-node third order triangle
//        {22, 12},  //! 12-node fourth order incomplete triangle
//        {23, 15},  //! 15-node fourth order triangle
//        {24, 15},  //! 15-node fifth order incomplete triangle
//        {25, 21},  //! 21-node fifth order complete triangle
//        {26,  4},  //! 4-node third order edge
//        {27,  5},  //! 5-node fourth order edge
//        {28,  6},  //! 6-node fifth order edge
//        {29, 20},  //! 20-node third order tetrahedron
//        {30, 35},  //! 35-node fourth order tetrahedron
//        {31, 56},  //! 56-node fifth order tetrahedron
//        {92, 64},  //! 64-node third order hexahedron
//        {93,125}   //! 125-node fourth order hexahedron
};

//! GmshTxtMeshReader : MeshReader
class GmshTxtMeshReader : protected MeshReader {

  //! Memory entry keys holding Mesh data
  typedef unordered_set<MemoryEntry*> MeshSet;

  protected:
    //! Constructor
    GmshTxtMeshReader(string filename, UnsMesh* mesh, Memory* memory) :
      MeshReader(filename, mesh, memory) {}

    //! Destructor: free mesh entries
    ~GmshTxtMeshReader();

    //! Read Gmsh mesh
    void read();

    //! Add new MeshSet entry (e.g. list of nodes, elements, node ids, etc.)
    template<class V> V* newEntry(size_t number,
                                  ValType value,
                                  VarType variable,
                                  string name,
                                  Bool plot = false,
                                  Bool restart = false) {
      // Allocate new memory entry
      MemoryEntry* entry = m_memory->newEntry(number,
                                              value,
                                              variable,
                                              name,
                                              plot,
                                              restart);
      // Store new MeshSet entry
      pair<MeshSet::iterator,Bool> n = m_meshEntry.insert(entry);
      if (!n.second)
        throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_INSERT);
      // Get pointer to new entry 
      return m_memory->getPtr<V>(entry);
    }

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

    //! Allocate memory to read mesh in
    void alloc();

    //! Add new element
    void addElem(Int type, vector<Int>& nodes);

    //! Add new element tags
    void addElemTags(Int type, vector<Int>& tags);

    //! Echo element tags and connectivity in all element sets
    void echoElemSets();

    MeshSet m_meshEntry;         //!< Memory entry keys for Mesh data

    Int m_nnodes = 0;            //!< Number of nodes
    Int m_nLins = 0;             //!< Number of line elements
    Int m_nTris = 0;             //!< Number of triangle elements
    Int m_nLinNodes = 0;         //!< Number of line element nodes
    Int m_nTriNodes = 0;         //!< Number of triangle element nodes
    Int m_nLinTags = 0;          //!< Number of line of tags
    Int m_nTriTags = 0;          //!< Number of triangle tags

    Int m_nodeCnt = 0;           //!< Counter for nodes added
    Int m_linCnt = 0;            //!< Counter for line elems added
    Int m_triCnt = 0;            //!< Counter for triangles added
    Int m_linNodeCnt = 0;        //!< Counter for line nodes added
    Int m_triNodeCnt = 0;        //!< Counter for triangle nodes added
    Int m_linTagCnt = 0;         //!< Counter for line tags added
    Int m_triTagCnt = 0;         //!< Counter for triangle tags added

    Int* m_node = nullptr;       //!< Node Ids
    Real* m_coord = nullptr;     //!< Node coordinates
    Int* m_linId = nullptr;      //!< Line element Ids
    Int* m_triId = nullptr;      //!< Triangle element Ids
    Int* m_lin = nullptr;        //!< Line elements
    Int* m_tri = nullptr;        //!< Triangle elements
    Int* m_linTag = nullptr;     //!< Line element tags
    Int* m_triTag = nullptr;     //!< Triangle element tags
};

} // namespace Quinoa

#endif // GmshTxtMeshReader_h
