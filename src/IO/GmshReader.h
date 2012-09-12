//******************************************************************************
/*!
  \file      src/IO/GmshReader.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 05:29:23 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshReader_h
#define GmshReader_h

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

//! GmshReader : MeshReader
class GmshReader : MeshReader {

  public:
    //! Constructor
    GmshReader(string filename, UnsMesh* mesh, Memory* memory) :
      MeshReader(filename, mesh, memory) {};

    //! Destructor
    ~GmshReader() {};

    //! Interface for read
    void read();

  private:
    //! Don't permit copy operator
    GmshReader(const GmshReader&);
    //! Don't permit assigment operator
    GmshReader& operator=(const GmshReader&);

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

#endif // GmshReader_h
