//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.h
  \author    J. Bakosi
  \date      Fri 21 Sep 2012 09:40:23 AM MDT
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
const string   NODEID_NAME = "nodeID";
const string    COORD_NAME = "coord";
const string   ELEMID_NAME = "elmID";
const string ELEMTYPE_NAME = "elmType";

//! Gmsh element types and their number of nodes,
//! all Gmsh-supported listed, Quinoa-supported uncommented
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
class GmshTxtMeshReader : MeshReader {

  //! Memory entry keys holding Mesh data
  typedef unordered_set<MemoryEntry*> MeshSet;

  protected:
    //! Constructor
    GmshTxtMeshReader(string filename, UnsMesh* mesh, Memory* memory) :
      MeshReader(filename, mesh, memory), m_nodesets(0), m_elemsets(0) {}

    //! Destructor: free mesh entries
    ~GmshTxtMeshReader();

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

    //! Read mandatory "$MeshFormat--$EndMeshFormat" section
    void readMeshFormat();

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();

    //! Add new element
    void addElem(vector<Int>& nodes);

    //! Add new element tags
    void addElemTags(vector<Int>& tags);

    //! Reserve element capacity
    void reserveElem(vector<vector<Int>>::size_type n);

    //! Echo element tags and connectivity in all element sets
    void echoElemSets();

    MeshSet m_meshEntry;         //!< Memory entry keys for Mesh data
    Int m_nodesets;              //!< Number of node sets
    Int m_elemsets;              //!< Number of element sets
    vector<vector<Int>> m_elem;  //!< Elements
    vector<vector<Int>> m_tag;   //!< Element tags
};

} // namespace Quinoa

#endif // GmshTxtMeshReader_h
