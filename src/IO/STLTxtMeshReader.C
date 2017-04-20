// *****************************************************************************
/*!
  \file      src/IO/STLTxtMeshReader.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     ASCII STL (STereoLithography) reader class definition
  \details   ASCII STL (STereoLithography) reader class definition.
*/
// *****************************************************************************

#include "Macro.h"
#include "STLMesh.h"
#include "STLTxtMeshReader.h"

using tk::STLTxtMeshReader;

STLTxtMeshReader::STLTxtMeshReader( const std::string filename, STLMesh& mesh )
  : Reader( filename ), m_mesh( mesh )
// *****************************************************************************
// Constructor
//! \param[in] filename File to read STL data from
//! \param[inout] mesh STLMesh object to store STL mesh
//! \author J. Bakosi
// *****************************************************************************
{
  // Set mesh name as filename modulo extension
  mesh.setName( filename.substr( 0, filename.find_last_of(".") ) );
}

void
STLTxtMeshReader::readMesh()
// *****************************************************************************
//  Read ASCII STL mesh
//! \author J. Bakosi
// *****************************************************************************
{
  // Count up number of vertices in STL mesh
  auto nnodes = readFacets( COUNT );
  Assert( nnodes%3 == 0, "Number of nodes in STL file must be divisible by 3");

  // Allocate memory to store coordinates and face list
  m_mesh.alloc( nnodes );

  // Read and store mesh
  readFacets( STORE, m_mesh.getx(), m_mesh.gety(), m_mesh.getz() );
}

std::size_t
STLTxtMeshReader::readFacets( const bool store,
                              tk::real* const x,
                              tk::real* const y,
                              tk::real* const z )
// *****************************************************************************
//  Read ASCII STL mesh
//  \param[in]  store  Whether to store the facets or not (i.e., only count)
//  \param[in]  x      Vertex x coordinates
//  \param[in]  y      Vertex y coordinates
//  \param[in]  z      Vertex z coordinates
//  \return            Number of vertices counted
//! \author J. Bakosi
// *****************************************************************************
{
  #ifdef STRICT_GNUC
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Winline"
  #endif
  // Define possible keywords in ASCI STL file: objects and their correct
  // string values (for error checking)
  STLKeyword solid("solid"), facet("facet"), normal("normal"), outer("outer"),
             loop("loop"), vertex("vertex"), endloop("endloop"),
             endfacet("endfacet");
  #ifdef STRICT_GNUC
    #pragma GCC diagnostic pop
  #endif

  // Read in solids with their facets until eof
  std::size_t num = 0;
  while ( !m_inFile.eof() ) {
    // Start reading new solid
    std::string solidname;
    m_inFile >> solid >> solidname;

    // Read and store facets
    bool newfacet = true;
    while (newfacet) {
      // Coordinates of normal, triangle vertex A, B, and C
      tk::real nx, ny, nz, Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz;

      // Read in normal (and throw it away as it is redundant)
      m_inFile >> facet >> normal >> nx >> ny >> nz;

      // Read in and store off triangle
      m_inFile >> outer >> loop;
      m_inFile >> vertex >> Ax >> Ay >> Az;
      m_inFile >> vertex >> Bx >> By >> Bz;
      m_inFile >> vertex >> Cx >> Cy >> Cz;
      m_inFile >> endloop >> endfacet;

      // Store coordinates of facet if requested
      if (store) {
        x[num] = Ax;
        y[num] = Ay;
        z[num] = Az;
        ++num;
        x[num] = Bx;
        y[num] = By;
        z[num] = Bz;
        ++num;
        x[num] = Cx;
        y[num] = Cy;
        z[num] = Cz;
        ++num;
      } else {
        num += 3;       // only increase number of vertices
      }

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
        back = m_inFile.tellg();
        m_inFile >> kw;              // try to read on
        if (!m_inFile.eof()) m_inFile.seekg(back);
      } else {
        Throw( "Corruption in ASCII STL file while parsing keyword '" + kw +
               "': keyword 'endfacet' must be followed by either 'facet' or "
               "'endsolid'" );
      }
    }   // while (newfacet)
  }   // while (newsolid)

  // Seek to beginning of file
  m_inFile.seekg(0, std::ios::beg);

  // Return number of vertices
  return num;
}
