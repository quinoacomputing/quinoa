//******************************************************************************
/*!
  \file      src/SDE/SDE.h
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 10:04:20 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE
  \details   SDE
*/
//******************************************************************************
#ifndef SDE_h
#define SDE_h

#include <cstdint>

#include <Model.h>
#include <Base.h>
#include <InitPolicy.h>
#include <LayoutPolicy.h>

namespace quinoa {

//! SDE
template< class Init, class Layout >
class SDE : public Model {

  protected:
    //! Constructor: protected, designed to be base-only
    explicit SDE( const Base& base,
                  tk::real* const particles,
                  int offset,
                  int ncomp ) :
      m_particles( particles ),
      m_npar( base.control.get< tag::component, tag::npar >() ),
      m_nprop( base.control.nprop() ),
      m_offset( offset ),
      m_ncomp( ncomp )
    {
      // Initialize particle properties (and throw away init policy
      Init initialize( m_particles, m_npar, m_nprop, m_offset, m_ncomp );
      // Instantiate RNG
      initRNG( base );
    }

    tk::real* const m_particles;    //!< Particle properties
    const uint64_t m_npar;          //!< Total number of particles
    const int m_nprop;              //!< Total number of particle properties
    const int m_offset;             //!< Offset SDE operates from
    const int m_ncomp;              //!< Number of components

  private:
    //! Don't permit copy constructor
    SDE(const SDE&) = delete;
    //! Don't permit copy assigment
    SDE& operator=(const SDE&) = delete;
    //! Don't permit move constructor
    SDE(SDE&&) = delete;
    //! Don't permit move assigment
    SDE& operator=(SDE&&) = delete;

    //! Instantiate random number genrator
    void initRNG( const Base& base ) {
      // Get vector of selected RNGs
//       const std::vector< tk::ctr::RNGType > rngs =
//         base.control.get< tag::selected, tk::tag::rng>();
//       // For now, only instantiate the first one of the RNGs
//       if (rngs[0] != tk::ctr::RNGType::NO_RNG) {
//         m_rng = std::unique_ptr< tk::RNG >( base.rng[rngs[0]]() );
//       }
//       ErrChk( m_rng, tk::ExceptType::FATAL, "No RNG requested");
    }

    std::unique_ptr< tk::RNG > m_rng;           //!< Random number generator
};

} // quinoa::

#endif // SDE_h
