//******************************************************************************
/*!
  \file      src/SDE/SDE.h
  \author    J. Bakosi
  \date      Tue 14 Jan 2014 08:37:24 PM MST
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

namespace quinoa {

//! SDE
template< class Init >
class SDE : public Model {

  protected:
    //! Constructor: protected, designed to be base-only
    explicit SDE( const Base& base,
                  tk::real* const particles,
                  int offset,
                  int ncomp ) :
      m_particles( particles ),
      m_npar( base.control.get< ctr::component, ctr::npar >() ),
      m_nprop( base.control.nprop() ),
      m_offset( offset ),
      m_ncomp( ncomp ),
      m_initialize( m_particles, m_npar, m_nprop, m_offset, m_ncomp ) {}

  private:
    //! Don't permit copy constructor
    SDE(const SDE&) = delete;
    //! Don't permit copy assigment
    SDE& operator=(const SDE&) = delete;
    //! Don't permit move constructor
    SDE(SDE&&) = delete;
    //! Don't permit move assigment
    SDE& operator=(SDE&&) = delete;

    tk::real* const m_particles;    //!< Particle properties
    const uint64_t m_npar;          //!< Total number of particles
    const int m_nprop;              //!< Total number of particle properties
    const int m_offset;             //!< Offset SDE operates from
    const int m_ncomp;              //!< Number of components

    Init m_initialize;              //!< Initialization policy
};

} // quinoa::

#endif // SDE_h
