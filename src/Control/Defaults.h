//******************************************************************************
/*!
  \file      src/Control/Defaults.h
  \author    J. Bakosi
  \date      Wed 13 Mar 2013 07:55:39 PM MDT
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
  1,                          //!< 13: Glob output interval
  "jpdf",                     //!< 14: Default jpdf base filename
  "glob",                     //!< 15: Default glob filename
  "plot",                     //!< 16: Default plot base filename
  vector<real>(),             //!< 17: Parameters 'b'
  vector<real>(),             //!< 18: Paramaters 'S'
  vector<real>(),             //!< 19: Parameters 'kappa'
  vector<real>(),             //!< 20: Parameters 'c_ij'
  vector<Product>()           //!< 21: Statistics
);

} // namespace control

} // namespace Quinoa

#endif // Defaults_h
