// *****************************************************************************
/*!
  \file      src/IO/H5PartWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     H5Part particles data writer
  \details   H5Part particles data writer class definition, facilitating writing
    particle coordinates and associated particle fields into HDF5-based H5Part
    data files in parallel, using MPI-IO.
*/
// *****************************************************************************

#include "H5PartWriter.h"
#include "Exception.h"

using tk::H5PartWriter;

H5PartWriter::H5PartWriter( const std::string& filename ) :
  m_filename( filename )
// *****************************************************************************
//  Constructor: create/open H5Part file
//! \param[in] filename File to open as H5Part file
//! \details It is okay to call this constructor with empty filename. In that
//!   case no IO will be performed. This is basically a punt to enable skipping
//!   H5Part I/O. Particles are a highly experimental feature at this point.
//! \note If the file exists, it will be truncated.
//! \author J. Bakosi
// *****************************************************************************
{
  if (m_filename.empty()) return;

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wold-style-cast"
  #endif

  auto f =
    H5PartOpenFileParallel(filename.c_str(), H5PART_WRITE, MPI_COMM_WORLD);

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #endif

  ErrChk( f, "Failed to create/open H5Part file: " + filename );

  ErrChk( H5PartWriteFileAttribString( f, "Origin", "Written by Quinoa" ) ==
          H5PART_SUCCESS, "Failed to write file attribute to " + filename );

  ErrChk( H5PartCloseFile( f ) == H5PART_SUCCESS,
          "Failed to close file " + filename );
}

void
H5PartWriter::writeCoords( uint64_t it,
                           const std::vector< tk::real >& x,
                           const std::vector< tk::real >& y,
                           const std::vector< tk::real >& z ) const
// *****************************************************************************
//  Write particle coordinates to H5Part file
//! \param[in] it Iteration number
//! \param[in] x X coordinates of particles
//! \param[in] y Y coordinates of particles
//! \param[in] z Z coordinates of particles
//! \author J. Bakosi
// *****************************************************************************
{
  if (m_filename.empty()) return;

  Assert( x.size() == y.size() && y.size() == z.size(),
          "Particle coordinates array sizes mismatch" );

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wold-style-cast"
  #endif

  auto f =
    H5PartOpenFileParallel(m_filename.c_str(), H5PART_APPEND, MPI_COMM_WORLD);

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #endif

  ErrChk( f, "Failed to open H5Part file for appending: " + m_filename );

  ErrChk( H5PartSetStep( f, static_cast<h5part_int64_t>(it) ) == H5PART_SUCCESS,
          "Failed to set time step in file " + m_filename );

  ErrChk( H5PartSetNumParticles( f, static_cast<h5part_int64_t>(x.size()) ) ==
          H5PART_SUCCESS, "Failed to set number of particles in file " +
                             m_filename );

  ErrChk( H5PartWriteDataFloat64( f, "x", x.data() ) == H5PART_SUCCESS,
          "Failed to write x particle coordinates to file " + m_filename );
  ErrChk( H5PartWriteDataFloat64( f, "y", y.data() ) == H5PART_SUCCESS,
          "Failed to write y particle coordinates to file " + m_filename );
  ErrChk( H5PartWriteDataFloat64( f, "z", z.data() ) == H5PART_SUCCESS,
          "Failed to write z particle coordinates to file " + m_filename );

  ErrChk( H5PartCloseFile( f ) == H5PART_SUCCESS,
          "Failed to close file " + m_filename );
}
