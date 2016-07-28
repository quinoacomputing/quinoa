// *****************************************************************************
/*!
  \file      src/IO/H5PartWriter.C
  \author    J. Bakosi
  \date      Thu 28 Jul 2016 09:01:24 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     H5Part particles data writer
  \details   H5Part particles data writer class definition, facilitating writing
    particle coordinates and associated particle fields into HDF5-based H5Part
    data files in parallel, using MPI-IO.
*/
// *****************************************************************************

#include <NoWarning/mpi.h>
#include <NoWarning/H5Part.h>
#include <NoWarning/H5Block.h>

#include "H5PartWriter.h"
#include "Exception.h"

using tk::H5PartWriter;

H5PartWriter::H5PartWriter( const std::string& filename ) :
  m_filename( filename ), m_outFile( nullptr )
// *****************************************************************************
//  Constructor: create/open H5Part file
//! \param[in] filename File to open as H5Part file
//! \author J. Bakosi
// *****************************************************************************
{
  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wold-style-cast"
  #endif

  m_outFile =
    H5PartOpenFileParallel( filename.c_str(), H5PART_WRITE, MPI_COMM_WORLD );

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #endif

  ErrChk( m_outFile != nullptr,
         "Failed to create/open H5Part file: " + filename );

  H5PartWriteFileAttribString( m_outFile, "Origin", "Written by Quinoa" );
}

H5PartWriter::~H5PartWriter() noexcept
// *****************************************************************************
//  Destructor
//! \author J. Bakosi
// *****************************************************************************
{
  H5PartCloseFile( m_outFile );
}
