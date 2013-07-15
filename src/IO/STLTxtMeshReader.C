//******************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.C
  \author    J. Bakosi
  \date      Sun 14 Jul 2013 08:41:43 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) reader class definition
  \details   ASCII STL (STereoLithography) reader class definition
*/
//******************************************************************************

#include <sstream>
#include <iostream>

#include <STLTxtMeshReader.h>

using namespace Quinoa;

void
STLTxtMeshReader::read()
//******************************************************************************
//  Read ASCII STL mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Define possible keywords in ASCI STL file: objects and their correct
  // string values (for error checking)
  STLKeyword solid("solid"), facet("facet"), normal("normal"), outer("outer"),
             loop("loop"), vertex("vertex"), endloop("endloop"),
             endfacet("endfacet");

  // Read in solids with their facets until eof
  while (!m_inFile.eof()) {
    // Start reading new solid
    std::string solidname;
    m_inFile >> solid >> solidname;

    // Read and store facets
    bool newfacet = true;
    while (newfacet) {
      real nx, ny, nz;               // Normal
      Triangle t;                    // Triangle

      // Read in normal (and throw it away as it is redundant)
      m_inFile >> facet >> normal >> nx >> ny >> nz;

      // Read in and store off triangle
      m_inFile >> outer >> loop;
      m_inFile >> vertex >> t.A.x >> t.A.y >> t.A.z;
      m_inFile >> vertex >> t.B.x >> t.B.y >> t.B.z;
      m_inFile >> vertex >> t.C.x >> t.C.y >> t.C.z;
      m_inFile >> endloop >> endfacet;
      m_Triangles.push_back(t);

      // Read in next keyword
      std::streampos back = m_inFile.tellg();     // save file position
      std::string kw;                // not an STLKeyword: no error checking
      m_inFile >> kw;
      if (kw == "facet") {
        m_inFile.seekg(back);        // seek back
        newfacet = true;             // there is more to this solid
      } else if (kw == "endsolid") {
        m_inFile >> solidname;       // read in solidname last time
        newfacet = false;            // solid finished, try to read next one
        std::streampos back = m_inFile.tellg();
        m_inFile >> kw;              // try to read on
        if (!m_inFile.eof()) m_inFile.seekg(back);
      } else {
        Throw(ExceptType::FATAL,
              "Corruption in ASCII STL file while parsing keyword '" + kw +
              "': keyword 'endfacet' must be followed by either 'facet' or "
              "'endsolid'");
      }
    }   // while (newfacet)
  }   // while (newsolid)

  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inFile.clear();
}
