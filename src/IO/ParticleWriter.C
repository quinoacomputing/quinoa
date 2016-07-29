// *****************************************************************************
/*!
  \file      src/IO/ParticleWriter.C
  \author    J. Bakosi
  \date      Fri 29 Jul 2016 03:04:02 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ group for outputing particle data to file via H5Part
  \details   Charm++ group for outputing particle data to file via H5Part in
     parallel using MPI-IO.
*/
// *****************************************************************************

#include "Exception.h"
#include "ParticleWriter.h"

using tk::ParticleWriter;

void
ParticleWriter::writeTimeStamp( uint64_t it, uint64_t npar )
// *****************************************************************************
//  Write a new time stamp to H5Part file
//! \param[in] it Iteration number
//! \param[in] npar Number of particles we will write in this iteration
//! \note This function does not have to be declared as a Charm++ entre method
//!    since it is always called by chares on the same PE.
//! \author J. Bakosi
// *****************************************************************************
{
  m_npar = npar;

  // Since potentiall this function may be called mutliple times (by multiple
  // chare array elements on our PE), make sure the function call to output the
  // new time stamp is only called once per PE
  if (!m_timestamped) m_writer.writeTimeStamp( it, npar );
  m_timestamped = true;
}

void
ParticleWriter::writeCoords( const std::vector< tk::real >& x,
                             const std::vector< tk::real >& y,
                             const std::vector< tk::real >& z )
// *****************************************************************************
// Receive and write particle coordinates
//! \param[in] x X coordinates of particles
//! \param[in] y Y coordinates of particles
//! \param[in] z Z coordinates of particles
//! \note This function does not have to be declared as a Charm++ entre method
//!    since it is always called by chares on the same PE.
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( m_timestamped,
          "Outputing the time stamp must precede the particle coordinates" );
  Assert( x.size() == y.size(), "Particle coordinates array sizes mismatch" );
  Assert( y.size() == z.size(), "Particle coordinates array sizes mismatch" );

  m_x.insert( end(m_x), begin(x), end(x) );
  m_y.insert( end(m_y), begin(y), end(y) );
  m_z.insert( end(m_z), begin(z), end(z) );

  std::cout << CkMyPe() << ": npar: " << m_npar << ", " << m_x.size() << '\n';

  if (m_x.size() == m_npar) {
    m_writer.writeCoords( m_x, m_y, m_z );
    m_x.clear();        // prepare for next step
    m_y.clear();
    m_z.clear();
    m_timestamped = false;
  }
}

#include "NoWarning/particlewriter.def.h"
