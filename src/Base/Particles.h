//******************************************************************************
/*!
  \file      src/Base/Particles.h
  \author    J. Bakosi
  \date      Sun 31 Jan 2016 07:16:07 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Particles used to store particle data.
  \details   Particles used to store data at particles as a specialization of
    DataLayout. See also Base/DataLayout.h and the rationale discussed in the
    [design](layout.html) document.
*/
//******************************************************************************
#ifndef Particles_h
#define Particles_h

#include "Config.h"
#include "DataLayout.h"

namespace tk {

//! Select data layout policy for particle data at compile-time
#if   defined PARTICLE_DATA_LAYOUT_AS_PARTICLE_MAJOR
using Particles = DataLayout< UnkEqComp >;
#elif defined PARTICLE_DATA_LAYOUT_AS_EQUATION_MAJOR
using Particles = DataLayout< EqCompUnk >;
#endif

} // tk::

#endif // Particles_h
