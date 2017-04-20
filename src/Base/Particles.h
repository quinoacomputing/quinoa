// *****************************************************************************
/*!
  \file      src/Base/Particles.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Particles used to store particle data.
  \details   Particles used to store data at particles as a specialization of
    tk::Data. See also Base/Data.h and the rationale discussed in the
    [design](layout.html) document.
*/
// *****************************************************************************
#ifndef Particles_h
#define Particles_h

#include "QuinoaConfig.h"
#include "Data.h"

namespace tk {

//! Select data layout policy for particle data at compile-time
#if   defined PARTICLE_DATA_LAYOUT_AS_PARTICLE_MAJOR
using Particles = Data< UnkEqComp >;
#elif defined PARTICLE_DATA_LAYOUT_AS_EQUATION_MAJOR
using Particles = Data< EqCompUnk >;
#endif

} // tk::

#endif // Particles_h
