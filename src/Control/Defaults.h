//******************************************************************************
/*!
  \file      src/Control/Defaults.h
  \author    J. Bakosi
  \date      Tue 30 Jul 2013 07:56:28 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Defaults for control
  \details   Defaults for control
*/
//******************************************************************************
#ifndef Defaults_h
#define Defaults_h

#include <limits>

#include <ControlTypes.h>

namespace Quinoa {

namespace control {

//! Default bundle for parsed data
const Bundle DEFAULTS(
  "",                                  //!< Title
  select::GeometryTypes::NO_GEOMETRY,  //!< Geometry definition
  select::PhysicsTypes::NO_PHYSICS,    //!< Physics
  select::PositionTypes::NO_POSITION,  //!< Position model
  select::MassTypes::NO_MASS,          //!< Mass model
  select::HydroTypes::NO_HYDRO,        //!< Hydrodynamics model
  select::EnergyTypes::NO_ENERGY,      //!< Internal energy model
  select::MixTypes::NO_MIX,            //!< Material mix model
  select::FrequencyTypes::NO_FREQUENCY,//!< Turbulence frequency model
  select::MixRateTypes::NO_MIXRATE,    //!< Material mix rate model
  select::RNGTestTypes::NO_RNGTEST,    //!< RNG test suite
  std::numeric_limits<uint64_t>::max(),//!< Number of time steps to take
  1.0,                                 //!< Time to terminate time stepping
  0.5,                                 //!< Size of time step
  0,                                   //!< Number of position components
  0,                                   //!< Number of density components
  0,                                   //!< Number of velocity components
  0,                                   //!< Number of scalar components
  0,                                   //!< Number of frequency components
  1,                                   //!< Total number of particles
  1,                                   //!< TTY output interval
  0,                                   //!< Dump output interval
  0,                                   //!< Plot output interval
  1,                                   //!< PDF output interval
  1,                                   //!< Glob output interval
  "",                                  //!< Input filename
  "",                                  //!< Output filename
  "jpdf",                              //!< Default jpdf filename
  "glob",                              //!< Default glob filename
  "stat",                              //!< Default statistics filename
  std::vector<real>(),                 //!< Parameters 'b'
  std::vector<real>(),                 //!< Paramaters 'S'
  std::vector<real>(),                 //!< Parameters 'kappa'
  std::vector<real>(),                 //!< Parameters 'c_ij'
  2.1,                                 //!< Parameter C0
  0.5,                                 //!< Parameter Atwood number
  0.5,                                 //!< Parameter C1 in gamma freq. model
  0.73,                                //!< Parameter C2 in gamma freq. model
  5.0,                                 //!< Parameter C3 in gamma freq. model
  0.25,                                //!< Parameter C4 in gamma freq. model
  std::vector<real>(),                 //!< Sextets for boxes for anal. geom.
  std::vector<Product>()               //!< Statistics
);

} // namespace control

} // namespace Quinoa

#endif // Defaults_h
