// *****************************************************************************
/*!
  \file      src/IO/RootMeshWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Root mesh-based data writer
  \details   Root mesh-based data writer class definition.
*/
// *****************************************************************************

#include <algorithm>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <cstdint>
#include <cstdio>

#include "RootMeshWriter.h"
#include "TFile.h"
#include "Exception.h"
#include "UnsMesh.h"

using tk::RootMeshWriter;

RootMeshWriter::RootMeshWriter( const std::string& filename,
                                        RootWriter mode,
                                        int cpuwordsize,
                                        int iowordsize ) :
  m_filename( filename ), m_outFile( 0 )
// *****************************************************************************
//  Constructor: create/open Root file
//! \param[in] filename File to open as Root file
//! \param[in] mode Root writer constructor mode: ExoWriter::CREATE for
//!   creating a new file, ExoWriter::OPEN for opening an existing file for
//!   appending
//! \param[in] cpuwordsize Set CPU word size, see Root documentation
//! \param[in] iowordsize Set I/O word size, see Root documentation
//! \author J. Bakosi
// *****************************************************************************
{
  ErrChk( m_outFile > 0, "Failed to create/open Root file: " + filename );
}

RootMeshWriter::~RootMeshWriter() noexcept
// *****************************************************************************
//  Destructor
//! \author J. Bakosi
// *****************************************************************************
{
}

void
RootMeshWriter::writeMesh( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write Root mesh file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  writeHeader( mesh );
}

void
RootMeshWriter::writeHeader( const UnsMesh& mesh ) const
// *****************************************************************************
//  Write Root header
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
}
