// *****************************************************************************
/*!
  \file      src/IO/H5PartWriter.C
  \author    J. Bakosi
  \date      Thu 28 Jul 2016 10:16:28 AM MDT
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

void
H5PartWriter::writeTimeStamp( uint64_t it, uint64_t npar ) const
// *****************************************************************************
//  Write a new time stamp to H5Part file
//! \param[in] it Iteration number
//! \param[in] npar Number of particles we will write in this iteration
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( npar > 0, "Attempting to write 0 particles into H5Part file" );

  H5PartSetStep( m_outFile, static_cast< h5part_int64_t >( it ) );
  H5PartSetNumParticles( m_outFile, static_cast< h5part_int64_t >( npar ) );
}

void
H5PartWriter::writeCoords( const std::vector< tk::real >& x,
                           const std::vector< tk::real >& y,
                           const std::vector< tk::real >& z ) const
// *****************************************************************************
//  Write particle coordinates to H5Part file
//! \param[in] x X coordinates of particles
//! \param[in] y Y coordinates of particles
//! \param[in] z Z coordinates of particles
//! \author J. Bakosi
// *****************************************************************************
{
  H5PartWriteDataFloat64( m_outFile, "x", x.data() );
  H5PartWriteDataFloat64( m_outFile, "y", y.data() );
  H5PartWriteDataFloat64( m_outFile, "z", z.data() );
}
