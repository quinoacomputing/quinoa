// *****************************************************************************
/*!
  \file      src/IO/ParticleWriter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ group for outputing particle data to file via H5Part
  \details   Charm++ group for outputing particle data to file via H5Part in
     parallel using MPI-IO.
*/
// *****************************************************************************

#include "ParticleWriter.hpp"
#include "Exception.hpp"

void
tk::ParticleWriter::writeCoords( uint64_t it,
                                 const std::vector< tk::real >& x,
                                 const std::vector< tk::real >& y,
                                 const std::vector< tk::real >& z,
                                 CkCallback c )
// *****************************************************************************
//  Write particle coordinates to file
//! \param[in] it Output iteration count
//! \param[in] x X coordinates of particles
//! \param[in] y Y coordinates of particles
//! \param[in] z Z coordinates of particles
//! \param[in] c Function to continue with after the write
// *****************************************************************************
{
  Assert( x.size() == y.size() && y.size() == z.size(),
          "Particle coordinates array sizes mismatch" );

  // WARNING: This will not work with virtualization, because the positions to
  // be written must be buffered up for that to work.
  m_writer.writeCoords( it, x, y, z );

  c.send();
}

#include "NoWarning/particlewriter.def.h"
