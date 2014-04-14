//******************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.C
  \author    J. Bakosi
  \date      Mon Apr 14 15:31:47 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ExodusII mesh writer
  \details   ExodusII mesh writer
*/
//******************************************************************************

#include <iostream>

#include <exodusII.h>
#include <ne_nemesisI.h>

#include <ExodusIIMeshWriter.h>
#include <Exception.h>

using quinoa::ExodusIIMeshWriter;

ExodusIIMeshWriter::ExodusIIMeshWriter( const std::string& filename,
                                        UnsMesh& mesh,
                                        int cpuwordsize,
                                        int iowordsize ) :
  m_filename(filename), m_mesh(mesh), m_outFile(0)
//******************************************************************************
//  Constructor: create Exodus II file
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile = ex_create( filename.c_str(),
                         EX_CLOBBER | EX_LARGE_MODEL,
                         &cpuwordsize,
                         &iowordsize );

  ErrChk( m_outFile > 0, tk::ExceptType::FATAL,
          "Failed to create file: " + filename );
}

ExodusIIMeshWriter::~ExodusIIMeshWriter()
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
  if ( ex_close(m_outFile) < 0 )
    std::cout << "WARNING: Failed to close file: " << m_filename << std::endl;
}

void
ExodusIIMeshWriter::write()
//******************************************************************************
//  Write ExodusII mesh file
//! \author J. Bakosi
//******************************************************************************
{
  writeHeader();
  writeNodes();
  writeElements();
}

void
ExodusIIMeshWriter::writeHeader()
//******************************************************************************
//  Write ExodusII header
//! \author J. Bakosi
//******************************************************************************
{
  ErrChk( ex_put_init( m_outFile, "Written by Quinoa", 3, m_mesh.nnode(),
                       m_mesh.tetinpoel().size(), 1, 0, 0 ) == 0,
          tk::ExceptType::FATAL,
          "Failed to write header to file: " + m_filename );
}

void
ExodusIIMeshWriter::writeNodes()
//******************************************************************************
//  Write node coordinates to ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  auto& c = m_mesh.coord();

  for (std::size_t n=0; n<m_mesh.nnode(); ++n) {
    ErrChk( ne_put_n_coord( m_outFile, n+1, 1, &c[n][0], &c[n][1],
                            &c[n][2] ) == 0,    // 1-based node ids
            tk::ExceptType::FATAL,
            "Failed to write coordinates to file: " + m_filename );
  }
}

void
ExodusIIMeshWriter::writeElements()
//******************************************************************************
//  Write element connectivity to ExodusII file
//! \author J. Bakosi
//******************************************************************************
{
  auto& c = m_mesh.tetinpoel();

  // Write element block information
  int elclass = 1;
  ErrChk( ex_put_elem_block( m_outFile, elclass, "TET", c.size(), 4, 0 ) == 0,
            tk::ExceptType::FATAL,
            "Failed to write element block to file: " + m_filename );

  // Write element connectivity
  for (std::size_t e=0; e<c.size(); ++e) {
    std::vector< int > n( c[e] );
    for (auto& i : n) ++i;      // 1-based node ids
    ErrChk( ne_put_n_elem_conn( m_outFile, elclass, e+1, 1, n.data() ) == 0,
            tk::ExceptType::FATAL,
            "Failed to write element connectivity to file: " + m_filename );
  }
}
