//******************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.C
  \author    J. Bakosi
  \date      Fri 12 Jul 2013 09:13:03 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) reader class definition
  \details   ASCII STL (STereoLithography) reader class definition
*/
//******************************************************************************

#include <STLTxtMeshReader.h>

using namespace Quinoa;

void
STLTxtMeshReader::read()
//******************************************************************************
//  Read ASCII STL mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Read header
  std::string solid;
  m_inFile >> solid >> m_name;
  ErrChk(solid == "solid", ExceptType::FATAL,
         "ASCII STL header corrupt: First line should be 'solid <name>'");

  // Read and store facets until line 'endsolid'
  bool endsolid = false;
  while (!endsolid) {
    //inMesh
  }
}
