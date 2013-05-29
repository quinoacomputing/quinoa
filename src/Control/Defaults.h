//******************************************************************************
/*!
  \file      src/Control/Defaults.h
  \author    J. Bakosi
  \date      Wed May 29 07:19:30 2013
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
  "",                                  //!< Title
  select::PhysicsTypes::NO_PHYSICS,    //!< Physics
  select::PositionTypes::NO_POSITION,  //!< Position model
  select::MassTypes::NO_MASS,          //!< Mass model
  select::HydroTypes::NO_HYDRO,        //!< Hydrodynamics model
  select::EnergyTypes::NO_ENERGY,      //!< Internal energy model
  select::MixTypes::NO_MIX,            //!< Material mix model
  select::FrequencyTypes::NO_FREQUENCY,//!< Material mix model
  numeric_limits<int>::max(),          //!< Number of time steps to take
  1.0,                                 //!< Time to terminate time stepping
  0.5,                                 //!< Size of time step
  0,                                   //!< Number of position components
  0,                                   //!< Number of density components
  0,                                   //!< Number of velocity components
  0,                                   //!< Number of scalar components
  1,                                   //!< Total number of particles
  1,                                   //!< TTY output interval
  0,                                   //!< Dump output interval
  0,                                   //!< Plot output interval
  1,                                   //!< PDF output interval
  1,                                   //!< Glob output interval
  "jpdf",                              //!< Default jpdf base filename
  "glob",                              //!< Default glob filename
  "plot",                              //!< Default plot base filename
  vector<real>(),                      //!< Parameters 'b'
  vector<real>(),                      //!< Paramaters 'S'
  vector<real>(),                      //!< Parameters 'kappa'
  vector<real>(),                      //!< Parameters 'c_ij'
  2.1,                                 //!< Parameter C0
  0.5,                                 //!< Parameter Atwood number
  vector<Product>()                    //!< Statistics
);

} // namespace control

} // namespace Quinoa

#endif // Defaults_h
