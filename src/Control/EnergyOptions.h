//******************************************************************************
/*!
  \file      src/Control/EnergyOptions.h
  \author    J. Bakosi
  \date      Fri Sep 27 08:52:53 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Energy model options and associations
  \details   Energy model options and associations
*/
//******************************************************************************
#ifndef EnergyOptions_h
#define EnergyOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace quinoa {
namespace sel {

//! Energy model types
enum class EnergyType : uint8_t { NO_ENERGY=0 };

//! Class with base templated on the above enum class with associations
class Energy : public Toggle<EnergyType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Energy() : Toggle<EnergyType>("Energy", names, values) {}

  private:
    //! Don't permit copy constructor
    Energy(const Energy&) = delete;
    //! Don't permit copy assigment
    Energy& operator=(const Energy&) = delete;
    //! Don't permit move constructor
    Energy(Energy&&) = delete;
    //! Don't permit move assigment
    Energy& operator=(Energy&&) = delete;

    //! Get access to energy keywords
    //const grm::kw::hommix hommix {};

    //! Enums -> names
    const std::map<EnergyType, std::string> names {
      { EnergyType::NO_ENERGY, "n/a" }
    };

    //! keywords -> Enums
    const std::map<std::string, EnergyType> values {
      { "no_energy", EnergyType::NO_ENERGY }
    };
};

} // sel::
} // quinoa::

#endif // EnergyOptions_h
