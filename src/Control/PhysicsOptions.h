//******************************************************************************
/*!
  \file      src/Control/PhysicsOptions.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 02:42:38 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics options and associations
  \details   Physics options and associations
*/
//******************************************************************************
#ifndef PhysicsOptions_h
#define PhysicsOptions_h

#include <map>
#include <sstream>

#include <Exception.h>

namespace Quinoa {

namespace select {

class Physics {

  public:
    //! Physics types
    enum class Enum { NO_PHYSICS=0,
                      HOMOGENEOUS_MIX,
                      HOMOGENEOUS_HYDRO,
                      HOMOGENEOUS_RAYLEIGH_TAYLOR,
                      SPINSFLOW };

    //! Operator << for writing Enum to output streams
    friend std::ostream& operator<< (std::ostream& os, const Enum& e) {
      os << static_cast<std::underlying_type<Enum>::type>(e);
      return os;
    }

    //! Operator + for adding Enum to a std::string
    friend std::string operator+ (const std::string& lhs, Enum e) {
      std::stringstream ss;
      ss << lhs << e ;
      std::string rhs = ss.str();
      return rhs;
    }

    //! Constructor initializing associations
    Physics() :
      //! Enums -> names
      names{ { Enum::NO_PHYSICS, "No physics" },
             { Enum::HOMOGENEOUS_MIX, "Homogeneous material mixing"},
             { Enum::HOMOGENEOUS_HYDRO, "Homogeneous hydrodynamics"},
             { Enum::HOMOGENEOUS_RAYLEIGH_TAYLOR,
               "Homogeneous Rayleigh-Taylor"},
             { Enum::SPINSFLOW,
               "Standalone-Particle Incompressible Navier-Stokes Flow"} },
      //! keywords -> Enums
      values{ { "no_physics", Enum::NO_PHYSICS },
              { "hommix",     Enum::HOMOGENEOUS_MIX },
              { "homhydro",   Enum::HOMOGENEOUS_HYDRO },
              { "homrt",      Enum::HOMOGENEOUS_RAYLEIGH_TAYLOR },
              { "spinsflow",  Enum::SPINSFLOW } } {}

    //! Lookup Enum value based on keyword
    Enum value(const std::string keyword) const {
      auto it = values.find(keyword);
      Assert(it != values.end(), FATAL,
            "Cannot find value for keyword \"" + keyword + "\"");
      return it->second;
    }

    //! Lookup option name based on Enum
    const std::string& name(Enum value) const {
      auto it = names.find(value);
      Assert(it != names.end(), FATAL,
             std::string("Cannot find name for value \"") + value + "\"");
      return it->second;
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

    std::map<Enum, std::string> names;
    std::map<std::string, Enum> values;
};

} // namespace select

} // namespace Quinoa

#endif // PhysicsOptions_h
