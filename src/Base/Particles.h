//******************************************************************************
/*!
  \file      src/Base/Particles.h
  \author    J. Bakosi
  \date      Sun 03 Apr 2016 10:07:21 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Particles used to store particle data.
  \details   Particles used to store data at particles as a specialization of
    DataLayout. See also Base/DataLayout.h and the rationale discussed in the
    [design](layout.html) document.
*/
//******************************************************************************
#ifndef Particles_h
#define Particles_h

#include "QuinoaConfig.h"
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
