// *****************************************************************************
/*!
  \file      src/IO/ParticleWriter.h
  \author    J. Bakosi
  \date      Fri 29 Jul 2016 02:53:49 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ group for outputing particle data to file via H5Part
  \details   Charm++ group for outputing particle data to file via H5Part in
     parallel using MPI-IO.
*/
// *****************************************************************************
#ifndef ParticleWriter_h
#define ParticleWriter_h

#include "H5PartWriter.h"

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
      m_timestamped( 0 ),
      m_nchare( 0 ),
      m_x(),
      m_y(),
      m_z() {}

    //! Chares register on my PE
    //! \note This function does not have to be declared as a Charm++ entre method
    //!    since it is always called by chares on the same PE.
    void checkin() { ++m_nchare; }

    //! Write a new time stamp to the particle file
    void writeTimeStamp( uint64_t it, uint64_t npar );

    //! Receive and write particle coordinates
    void writeCoords( const std::vector< tk::real >& x,
                      const std::vector< tk::real >& y,
                      const std::vector< tk::real >& z );

  private:
    tk::H5PartWriter m_writer;      //!< Particle file format writer
    bool m_timestamped;  //!< Whether time stamp has been written out on this PE
    int m_nchare;        //!< Nnumber of chares contributing to my PE
    uint64_t m_npar;               //!< Number of particles to be written
    std::vector< tk::real > m_x;   //!< Buffer collecting x coordinates
    std::vector< tk::real > m_y;   //!< Buffer collecting y coordinates
    std::vector< tk::real > m_z;   //!< Buffer collecting z coordinates
};

} // tk::

#endif // ParticleWriter_h
