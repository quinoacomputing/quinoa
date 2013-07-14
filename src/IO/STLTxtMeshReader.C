//******************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.C
  \author    J. Bakosi
  \date      Sat 13 Jul 2013 10:31:33 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) reader class definition
  \details   ASCII STL (STereoLithography) reader class definition
*/
//******************************************************************************

#include <sstream>

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
         "Corrupt ASCII STL header: First line of input file '" +  m_filename +
         "' should be 'solid <model_name>'");

  // Read and store facets until line 'endsolid'
  bool endsolid = false;
  //while (!endsolid) {
{
    // possible keywords in ASCI STL file
    std::string facet, normal, outer_loop, vertex, end_loop, end_facet;
    // normal
    real nx, ny, nz;
    // three vertices of a triangle
    Triangle t;

    // Read in normal (and throw it away as it is redundant)
    m_inFile >> facet >> normal >> nx >> ny >> nz;
    ErrChk(facet == "facet", ExceptType::FATAL,
           "Corruption in ASCII STL file '" + m_filename +
           "' while parsing keyword '" + facet + "'");
    ErrChk(normal == "normal", ExceptType::FATAL,
           "Corruption in ASCII STL file '" + m_filename +
           "' while parsing keyword '" + normal + "'");

    // Read in triangle
  }
}
