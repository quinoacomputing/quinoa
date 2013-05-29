//******************************************************************************
/*!
  \file      src/Control/PhysicsOptions.h
  \author    J. Bakosi
  \date      Wed May 29 07:33:05 2013
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
enum class PhysicsTypes : uint8_t { NO_PHYSICS=0,
                                    HOMOGENEOUS_MIX,
                                    HOMOGENEOUS_HYDRO,
                                    HOMOGENEOUS_RAYLEIGH_TAYLOR,
                                    SPINSFLOW };

//! Class with base templated on the above enum class with associations
class Physics : public Toggle<PhysicsTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    explicit Physics() : Toggle<PhysicsTypes>(names, values) {
      //! Enums -> names
      names[PhysicsTypes::NO_PHYSICS] = "No physics";
      names[PhysicsTypes::HOMOGENEOUS_MIX] = "Homogeneous material mixing";
      names[PhysicsTypes::HOMOGENEOUS_HYDRO] = "Homogeneous hydrodynamics";
      names[PhysicsTypes::HOMOGENEOUS_RAYLEIGH_TAYLOR] =
        "Homogeneous Rayleigh-Taylor";
      names[PhysicsTypes::SPINSFLOW] =
        "Standalone-Particle Incompressible Navier-Stokes Flow";
      //! keywords -> Enums
      values["no_physics"] = PhysicsTypes::NO_PHYSICS;
      values["hommix"] = PhysicsTypes::HOMOGENEOUS_MIX;
      values["homhydro"] = PhysicsTypes::HOMOGENEOUS_HYDRO;
      values["homrt"] = PhysicsTypes::HOMOGENEOUS_RAYLEIGH_TAYLOR;
      values["spinsflow"] = PhysicsTypes::SPINSFLOW;
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

    std::map<PhysicsTypes, std::string> names;
    std::map<std::string, PhysicsTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // PhysicsOptions_h
