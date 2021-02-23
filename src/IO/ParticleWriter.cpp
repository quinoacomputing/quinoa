// *****************************************************************************
/*!
  \file      src/IO/ParticleWriter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ group for outputing particle data to file via H5Part
  \details   Charm++ group for outputing particle data to file via H5Part in
     parallel using MPI-IO.
*/
// *****************************************************************************

#include "ParticleWriter.hpp"
#include "Exception.hpp"

using tk::ParticleWriter;

void
ParticleWriter::npar( std::size_t n, CkCallback c )
// *****************************************************************************
//  Chares contribute their number of particles they will output on my node
//! \param[in] n Number of particles will be contributed
//! \param[in] c Function to continue with
// *****************************************************************************
{
  m_npar += n;
  c.send();
}

void
ParticleWriter::writeCoords( uint64_t it,
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
//! \param[in] c Function to continue with after the write is complete
// *****************************************************************************
{
  Assert( x.size() == y.size() && y.size() == z.size(),
          "Particle coordinates array sizes mismatch" );

  // buffer up coordinates
  m_x.insert( end(m_x), begin(x), end(x) );
  m_y.insert( end(m_y), begin(y), end(y) );
  m_z.insert( end(m_z), begin(z), end(z) );

  // if received from all chares on my PE, write to file
  if (m_x.size() == m_npar) {
    m_writer.writeCoords( it, m_x, m_y, m_z );
    m_npar = 0;
    m_x.clear();
    m_y.clear();
    m_z.clear();
  }

  c.send();
}

#include "NoWarning/particlewriter.def.h"
