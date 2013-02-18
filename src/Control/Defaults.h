//******************************************************************************
/*!
  \file      src/Control/Defaults.h
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 03:54:07 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Defaults for control
  \details   Defaults for control
*/
//******************************************************************************
#ifndef Defaults_h
#define Defaults_h

#include <limits>

#include <ControlTypes.h>

using namespace std;

namespace Quinoa {

namespace control {

//! Default bundle for parsed data
const Bundle Defaults(
  "",                         //!< 0: Title
  NO_PHYSICS,                 //!< 1: Physics
  NO_HYDRO,                   //!< 2: Hydrodynamics model
  NO_MIX,                     //!< 3: Material mix model
  numeric_limits<int>::max(), //!< 4: Number of time steps to take
  1.0,                        //!< 5: Time to terminate time stepping
  0.5,                        //!< 6: Size of time step
  1,                          //!< 7: Number of mixing scalars
  1,                          //!< 8: Total number of particles
  1                           //!< 9: One-liner info every few time steps
);

} // namespace control

} // namespace Quinoa

#endif // Defaults_h
