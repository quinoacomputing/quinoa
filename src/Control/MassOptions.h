//******************************************************************************
/*!
  \file      src/Control/MassOptions.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:03:05 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mass model options and associations
  \details   Mass model options and associations
*/
//******************************************************************************
#ifndef MassOptions_h
#define MassOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace quinoa {

namespace select {

//! Mass model types
enum class MassType : uint8_t { NO_MASS=0,
                                BETA };

//! Class with base templated on the above enum class with associations
class Mass : public Toggle<MassType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Mass() : Toggle<MassType>(names, values) {}

  private:
    //! Don't permit copy constructor
    Mass(const Mass&) = delete;
    //! Don't permit copy assigment
    Mass& operator=(const Mass&) = delete;
    //! Don't permit move constructor
    Mass(Mass&&) = delete;
    //! Don't permit move assigment
    Mass& operator=(Mass&&) = delete;

    //! Enums -> names
    const std::map<MassType, std::string> names {
      { MassType::NO_MASS, "No mass" },
      { MassType::BETA, "Beta" }
    };

    //! keywords -> Enums
    const std::map<std::string, MassType> values {
      { "no_mass", MassType::NO_MASS },
      { "mass_beta", MassType::BETA }
    };
};

} // namespace select

} // namespace quinoa

#endif // MassOptions_h
