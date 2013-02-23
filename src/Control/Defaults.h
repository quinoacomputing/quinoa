//******************************************************************************
/*!
  \file      src/Control/Defaults.h
  \author    J. Bakosi
  \date      Sat 23 Feb 2013 11:46:01 AM MST
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

// ICC: This could be inserted below once initializer lists are properly
// supported
const vector<real> ZEROVEC;

//! Default bundle for parsed data
const Bundle DEFAULTS(
  "",                         //!<  0: Title
  NO_PHYSICS,                 //!<  1: Physics
  NO_HYDRO,                   //!<  2: Hydrodynamics model
  NO_MIX,                     //!<  3: Material mix model
  numeric_limits<int>::max(), //!<  4: Number of time steps to take
  1.0,                        //!<  5: Time to terminate time stepping
  0.5,                        //!<  6: Size of time step
  1,                          //!<  7: Number of mixing scalars
  1,                          //!<  8: Total number of particles
  1,                          //!<  9: TTY output interval
  0,                          //!< 10: Dump output interval
  0,                          //!< 11: Plot output interval
  1,                          //!< 12: PDF output interval
  ZEROVEC,                    //!< 13: Vector of parameters 'b'
  ZEROVEC,                    //!< 14: Vector of paramaters 'S'
  ZEROVEC,                    //!< 15: Vector of parameters 'kappa'
  ZEROVEC                     //!< 16: Vector of parameters 'c_ij'
);

//! Default basefilenames for output
const string JPDF_FILENAME_BASE = "jpdf";

} // namespace control

} // namespace Quinoa

#endif // Defaults_h
