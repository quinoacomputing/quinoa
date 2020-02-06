// *****************************************************************************
/*!
  \file      src/IO/ParticleWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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
      m_writer( filename ) {}

    //! Write particle coordinates to file
    void writeCoords( uint64_t it,
                      const std::vector< tk::real >& x,
                      const std::vector< tk::real >& y,
                      const std::vector< tk::real >& z,
                      CkCallback c );

  private:
    tk::H5PartWriter m_writer;     //!< Particle file format writer
};

} // tk::

#endif // ParticleWriter_h
