//******************************************************************************
/*!
  \file      src/Mesh/GmshTxtMeshReader2D.C
  \author    J. Bakosi
  \date      Fri 21 Sep 2012 01:02:44 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition
*/
//******************************************************************************

#include <GmshTxtMeshReader2D.h>

using namespace Quinoa;

void
GmshTxtMeshReader2D::read()
//******************************************************************************
//  Public interface for read 2D txt Gmsh mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Read in txt Gmsh mesh into STL containers
  GmshTxtMeshReader::read();
  // Compress STL containers into memory entries
  compress();
}

void
GmshTxtMeshReader2D::compress()
//******************************************************************************
//  Compress mesh data from STL to memory entries
//! \author J. Bakosi
//******************************************************************************
{

}
