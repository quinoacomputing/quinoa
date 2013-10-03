//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Physics.h
  \author    J. Bakosi
  \date      Thu Oct  3 17:40:44 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics options and associations
  \details   Physics options and associations
*/
//******************************************************************************
#ifndef QuinoaPhysicsOptions_h
#define QuinoaPhysicsOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace sel {

//! Physics types
enum class PhysicsType : uint8_t { NO_PHYSICS=0,
                                   HOMOGENEOUS_MIX,
                                   HOMOGENEOUS_HYDRO,
                                   HOMOGENEOUS_RAYLEIGH_TAYLOR,
                                   SPINSFLOW };

//! Class with base templated on the above enum class with associations
class Physics : public Toggle<PhysicsType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Physics() : Toggle<PhysicsType>("Physics", names, values) {}

  private:
    //! Don't permit copy constructor
    Physics(const Physics&) = delete;
    //! Don't permit copy assigment
    Physics& operator=(const Physics&) = delete;
    //! Don't permit move constructor
    Physics(Physics&&) = delete;
    //! Don't permit move assigment
    Physics& operator=(Physics&&) = delete;

    //! Get access to physics keywords
    const kw::hommix hommix {};
    const kw::homhydro homhydro {};
    const kw::homrt homrt {};
    const kw::spinsflow spinsflow {};

    //! Enums -> names
    const std::map<PhysicsType, std::string> names {
      { PhysicsType::NO_PHYSICS, "n/a" },
      { PhysicsType::HOMOGENEOUS_MIX, hommix.name() },
      { PhysicsType::HOMOGENEOUS_HYDRO, homhydro.name() },
      { PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR, homrt.name() },
      { PhysicsType::SPINSFLOW, spinsflow.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, PhysicsType> values {
      { "no_physics", PhysicsType::NO_PHYSICS },
      { hommix.string(), PhysicsType::HOMOGENEOUS_MIX },
      { homhydro.string(), PhysicsType::HOMOGENEOUS_HYDRO },
      { homrt.string(), PhysicsType::HOMOGENEOUS_RAYLEIGH_TAYLOR },
      { spinsflow.string(), PhysicsType::SPINSFLOW }
    };
};

} // sel::
} // quinoa::

#endif // QuinoaPhysicsOptions_h
