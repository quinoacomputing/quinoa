//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.h
  \author    J. Bakosi
  \date      Tue 18 Sep 2012 09:26:42 PM MDT
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

//! Gmsh element types and their number of nodes
const map<Int, Int> GmshElemNodes = {
          { 1,  2},  //! 2-node line
          { 2,  3},  //! 3-node triangle
          { 3,  4},  //! 4-node quadrangle
          { 4,  4},  //! 4-node tetrahedron
          { 5,  8},  //! 8-node hexahedron
          { 6,  6},  //! 6-node prism
          { 7,  5},  //! 5-node pyramid
          { 8,  3},  //! 3-node second order line
          { 9,  6},  //! 6-node second order triangle
          {10,  9},  //! 9-node second order quadrangle
          {11, 10},  //! 10-node second order tetrahedron
          {12, 27},  //! 27-node second order hexahedron
          {13, 18},  //! 18-node second order prism
          {14, 14},  //! 14-node second order pyramid
          {15,  1},  //! 1-node point
          {16,  8},  //! 8-node second order quadrangle
          {17, 20},  //! 20-node second order hexahedron
          {18, 15},  //! 15-node second order prism
          {19, 13},  //! 13-node second order pyramid
          {20,  9},  //! 9-node third order incomplete triangle
          {21, 10},  //! 10-node third order triangle
          {22, 12},  //! 12-node fourth order incomplete triangle
          {23, 15},  //! 15-node fourth order triangle
          {24, 15},  //! 15-node fifth order incomplete triangle
          {25, 21},  //! 21-node fifth order complete triangle
          {26,  4},  //! 4-node third order edge
          {27,  5},  //! 5-node fourth order edge
          {28,  6},  //! 6-node fifth order edge
          {29, 20},  //! 20-node third order tetrahedron
          {30, 35},  //! 35-node fourth order tetrahedron
          {31, 56},  //! 56-node fifth order tetrahedron
          {92, 64},  //! 64-node third order hexahedron
          {93,125}   //! 125-node fourth order hexahedron
};

//! GmshTxtMeshReader : MeshReader
class GmshTxtMeshReader : MeshReader {

  public:
    //! Constructor
    GmshTxtMeshReader(string filename, UnsMesh* mesh, Memory* memory) :
      MeshReader(filename, mesh, memory) {}

    //! Destructor, default compiler generated
    ~GmshTxtMeshReader() = default;

    //! Public interface for read Gmsh mesh
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

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();
};

} // namespace Quinoa

#endif // GmshTxtMeshReader_h
