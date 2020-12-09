// *****************************************************************************
/*!
  \file      src/Base/Particles.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Particles used to store particle data.
  \details   Particles used to store data at particles as a specialization of
    tk::Data. See also Base/Data.h and the rationale discussed in the
    [design](layout.html) document.
*/
// *****************************************************************************
#ifndef Particles_h
#define Particles_h

#include "QuinoaBuildConfig.hpp"
#include "Data.hpp"

namespace tk {

//! Select data layout policy for particle data at compile-time
#if   defined PARTICLE_DATA_LAYOUT_AS_PARTICLE_MAJOR
using Particles = Data< UnkEqComp >;
#elif defined PARTICLE_DATA_LAYOUT_AS_EQUATION_MAJOR
using Particles = Data< EqCompUnk >;
#endif

} // tk::

#endif // Particles_h
