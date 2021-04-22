// *****************************************************************************
/*!
  \file      src/IO/ParticleWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ group for outputing particle data to file via H5Part
  \details   Charm++ group for outputing particle data to file via H5Part in
     parallel using MPI-IO.
*/
// *****************************************************************************
#ifndef ParticleWriter_h
#define ParticleWriter_h

#include <string>
#include <vector>

#include "H5PartWriter.hpp"

#include "NoWarning/particlewriter.decl.h"

namespace tk {

//! \brief Charm++ group used to output particle data to file in parallel using
//!   H5Part and MPI-IO
class ParticleWriter : public CBase_ParticleWriter {

  public:
    //! Constructor
    //! \param[in] filename Filename of particle output file
    explicit ParticleWriter( const std::string& filename ) :
      m_writer( filename ),
      m_npar( 0 ),
      m_x(),
      m_y(),
      m_z() {}

    //! Chares contribute their number of particles they will output on my node
    void npar( std::size_t n, CkCallback c );

    //! Write particle coordinates to file
    void writeCoords( uint64_t it,
                      const std::vector< tk::real >& x,
                      const std::vector< tk::real >& y,
                      const std::vector< tk::real >& z,
                      CkCallback c );

  private:
    tk::H5PartWriter m_writer;     //!< Particle file format writer
    std::size_t m_npar;            //!< Number of particles to be written
    std::vector< tk::real > m_x;   //!< Buffer collecting x coordinates
    std::vector< tk::real > m_y;   //!< Buffer collecting y coordinates
    std::vector< tk::real > m_z;   //!< Buffer collecting z coordinates
};

} // tk::

#endif // ParticleWriter_h
