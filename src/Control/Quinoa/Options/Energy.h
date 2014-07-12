//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Energy.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 04:33:29 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Energy model options and associations
  \details   Energy model options and associations
*/
//******************************************************************************
#ifndef QuinoaEnergyOptions_h
#define QuinoaEnergyOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Energy model types
enum class EnergyType : uint8_t { NO_ENERGY=0 };

//! Class with base templated on the above enum class with associations
class Energy : public tk::Toggle<EnergyType> {

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
    //const kw::hommix hommix {};

    //! Enums -> names
    const std::map<EnergyType, std::string> names {
      { EnergyType::NO_ENERGY, "n/a" }
    };

    //! keywords -> Enums
    const std::map<std::string, EnergyType> values {
      { "no_energy", EnergyType::NO_ENERGY }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaEnergyOptions_h
