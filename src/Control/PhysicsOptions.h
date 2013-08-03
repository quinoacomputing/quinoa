//******************************************************************************
/*!
  \file      src/Control/PhysicsOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:42:29 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics options and associations
  \details   Physics options and associations
*/
//******************************************************************************
#ifndef PhysicsOptions_h
#define PhysicsOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Physics types
enum class PhysicsType : uint8_t { NO_PHYSICS=0,
                                   HOMOGENEOUS_MIX,
                                   HOMOGENEOUS_HYDRO,
                                   HOMOGENEOUS_RAYLEIGH_TAYLOR,
                                   SPINSFLOW,
                                   RNGTEST };

//! Class with base templated on the above enum class with associations
class Physics : public Toggle<PhysicsType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    explicit Physics() : Toggle<PhysicsType>(names, values) {
      //! Enums -> names
      names[PhysicsType::NO_PHYSICS] = "No physics";
      names[PhysicsType::HOMOGENEOUS_MIX] = "Homogeneous material mixing";
      names[PhysicsType::HOMOGENEOUS_HYDRO] = "Homogeneous hydrodynamics";
      names[PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR] =
        "Homogeneous Rayleigh-Taylor";
      names[PhysicsType::SPINSFLOW] =
        "Standalone-Particle Incompressible Navier-Stokes Flow";
      names[PhysicsType::RNGTEST] = "Random number generator tests";
      //! keywords -> Enums
      values["no_physics"] = PhysicsType::NO_PHYSICS;
      values["hommix"] = PhysicsType::HOMOGENEOUS_MIX;
      values["homhydro"] = PhysicsType::HOMOGENEOUS_HYDRO;
      values["homrt"] = PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR;
      values["spinsflow"] = PhysicsType::SPINSFLOW;
      values["rngtest"] = PhysicsType::RNGTEST;
    }

  private:
    //! Don't permit copy constructor
    Physics(const Physics&) = delete;
    //! Don't permit copy assigment
    Physics& operator=(const Physics&) = delete;
    //! Don't permit move constructor
    Physics(Physics&&) = delete;
    //! Don't permit move assigment
    Physics& operator=(Physics&&) = delete;

    std::map<PhysicsType, std::string> names;
    std::map<std::string, PhysicsType> values;
};

} // namespace select

} // namespace Quinoa

#endif // PhysicsOptions_h
